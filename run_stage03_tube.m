% =========================================================================
% run_stage03_tube
% =========================================================================
clear; clc; close all;
if ~exist('OCTAVE_VERSION','builtin'), error('This simulation project targets GNU Octave.'); end

%% 1. Resolve project path & Diary
script_file = mfilename('fullpath');
if isempty(script_file), project_dir = pwd; else, project_dir = fileparts(script_file); end
addpath(project_dir); rehash();

out_dir = fullfile(project_dir,'results');
log_file = fullfile(out_dir, 'stage03_log.txt');
if exist(log_file, 'file'), delete(log_file); end
diary(log_file);

fprintf('=== Stage 3: uncertainty and invariant Tube ===\n');

%% 2. External numerical dependency
try, pkg('load','control'); catch ME, error('Stage 3 requires the Octave control package: %s',ME.message); end
if exist('dlqr','file') == 0 || exist('dlyap','file') == 0, error('Stage 3: dlqr or dlyap is unavailable.'); end
if exist('glpk','file') == 0, error('Stage 3: GLPK is unavailable.'); end

%% 3. Check local dependencies
required = {'apply_preprocessing', 'project_scores', 'collect_model_residuals', 'estimate_uncertainty', ...
    'compute_state_error_set', 'design_tube_feedback', 'compute_invariant_ellipsoid', 'compute_tightening', 'check_tightened_feasibility'};
for j = 1:numel(required)
    expected = fullfile(project_dir,[required{j},'.m']);
    if exist(expected,'file') ~= 2, error('Missing required file: %s',expected); end
    resolved = which(required{j});
    [resolved_dir,~,~] = fileparts(resolved);
    if ~strcmp(resolved_dir,project_dir), error('Function %s resolves outside this project: %s', required{j},resolved); end
end

%% 4. Load identification results
stage02_file = fullfile(out_dir,'stage02_identification.mat');
if exist(stage02_file,'file') ~= 2, error('Run run_stage02_identification.m before Stage 3.'); end
S = load(stage02_file);
Params = S.Params; Reference = S.Reference; Model = S.Model; Prep = S.Prep;

if ~strcmp(Params.code_version,'cryst_modular_v1'), error('Stage 3: incompatible artifact version.'); end
if Params.A_lv ~= Model.A_lv || Params.n_arx ~= Model.n_arx, error('Stage 3: inconsistent effective model dimensions.'); end

%% 5. Recompute residuals
Training_residuals = collect_model_residuals(S.Train_data,Model,Prep,Params);
Selection_residuals = collect_model_residuals(S.Val_data,Model,Prep,Params);
Calibration_residuals = collect_model_residuals(S.Cal_data,Model,Prep,Params);
UncSets = estimate_uncertainty(Calibration_residuals,Training_residuals,Selection_residuals,Model,Params);

%% 6. Initial error and full MIMO feedback
StateError = compute_state_error_set(Model,Params);
Feedback = design_tube_feedback(Model,Params);

%% 7. Compute invariant ellipsoid and physical tightening
Invariant = compute_invariant_ellipsoid(Model,Feedback,UncSets,StateError,Params);
Tube = compute_tightening(Model,Feedback,Invariant,UncSets,StateError,Params);

%% 8. Check sampled initial errors against the Tube
L = chol(Tube.P,'lower');
sets = {Training_residuals, Selection_residuals, Calibration_residuals};
labels = {'training','selection','calibration'};
InitialCoverage = struct();

fprintf('\nObserved initial-error coverage in the Tube\n');
for j = 1:numel(sets)
    E0 = sets{j}.e_initial; whitened = L\E0; quadratic = sum(whitened.^2,1);
    report.coverage = mean(quadratic <= 1+1e-7);
    report.max_quadratic = max(quadratic);
    report.n_batches = size(E0,2);
    InitialCoverage.(labels{j}) = report;
    fprintf('  %-11s: coverage=%6.2f%%, max quadratic=%.6g.\n', labels{j},100*report.coverage,report.max_quadratic);
    if report.max_quadratic > 1+1e-6
        error('Stage 3: sampled initial errors exceed the constructed initial Tube.');
    end
end

%% 9. Complete nominal trajectory feasibility
Feasibility = check_tightened_feasibility(Model,Tube,Params);
Tube.nominal_feasible = Feasibility.feasible;
if Feasibility.feasible, Tube.feasible_seed = Feasibility.t_feasible; else, Tube.feasible_seed = []; end

%% 10. Save stage artifact
Artifact = struct('stage',3,'code_version',Params.code_version,'octave_version',OCTAVE_VERSION, ...
    'created_at',datestr(now,31),'state_convention','observed_score_history','ready_for_nominal_optimization',Feasibility.feasible);

CV = S.CV; Validation = S.Validation; Calibration = S.Calibration;
stage03_file = fullfile(out_dir,'stage03_tube.mat');
save('-mat7-binary',stage03_file, 'Params','Reference','Prep','Model', 'CV','Validation','Calibration', ...
    'Training_residuals','Selection_residuals','Calibration_residuals', ...
    'UncSets','StateError','Feedback','Invariant','Tube','InitialCoverage','Feasibility','Artifact');
fprintf('\nSaved: %s\n',stage03_file);

if ~Feasibility.feasible, error('Stage 3 diagnostics were saved, but tightened nominal trajectory was not certified feasible.'); end

%% 11. Diagnostic plots
if Params.plot_results
    figure('Name','Stage03 uncertainty bounds', 'Color','w','Position',[100,100,1100,500]);
    subplot(1,2,1); bar(UncSets.eta_bound); grid on; xlabel('Latent component'); ylabel('Innovation halfwidth'); title('Calibration innovation bounds');
    subplot(1,2,2); bar(UncSets.r_bound); grid on; xlabel('Physical output channel'); ylabel('Standardized residual halfwidth'); title('Output reconstruction bounds');

    figure('Name','Stage03 input tightening', 'Color','w','Position',[130,130,1100,500]);
    ni = Params.n_in_opt; upper_margin = Tube.delta_v(1:ni); lower_margin = Tube.delta_v(ni+1:end);
    physical_width = Params.u_phys_max(Model.input_idx) - Params.u_phys_min(Model.input_idx);
    consumed_fraction = (upper_margin+lower_margin)./physical_width;
    bar(100*consumed_fraction); grid on; xlabel('Optimized input channel'); ylabel('Input interval removed (%)'); title('Fraction of physical input interval used by Tube tightening');
end

fprintf('\nStage 3 completed.\n');
fprintf('  Closed-loop spectral radius: %.8f.\n',Tube.spectral_radius);
fprintf('  Invariant containment ratio: %.10f.\n',Tube.invariant_ratio);
fprintf('  Complete tightened nominal trajectory: feasible.\n');
diary off;
