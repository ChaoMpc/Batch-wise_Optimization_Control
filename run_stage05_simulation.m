% =========================================================================
% run_stage05_simulation
% =========================================================================
clear; clc; close all;
if ~exist('OCTAVE_VERSION','builtin'), error('This simulation project targets GNU Octave.'); end

%% 1. Project path & Diary
script_file = mfilename('fullpath');
if isempty(script_file), project_dir = pwd; else, project_dir = fileparts(script_file); end
addpath(project_dir); rehash();

out_dir = fullfile(project_dir,'results');
log_file = fullfile(out_dir, 'stage05_log.txt');
if exist(log_file, 'file'), delete(log_file); end
diary(log_file);

fprintf('=== Stage 5: nonlinear online simulation ===\n');

%% 2. Check dependencies
required = {'cryst_plant', 'cryst_measure', 'generate_simulation_scenario', 'assemble_latent_state', ...
    'initialize_score_history', 'update_score_history', 'simulate_control_batch', 'run_ilmpc', ...
    'run_xue_ilmpc_v2', 'run_oh_ilmpc', 'compute_metrics', 'report_simulation_comparison', ...
    'solve_batch_update', 'solve_qp_checked', 'evaluate_nominal_trajectory'};
for j = 1:numel(required)
    expected = fullfile(project_dir,[required{j},'.m']);
    if exist(expected,'file') ~= 2, error('Missing required file: %s',expected); end
end

%% 3. Load artifacts
nominal_file = fullfile(out_dir,'stage04_nominal.mat');
tube_file = fullfile(out_dir,'stage03_tube.mat');
if exist(nominal_file,'file') ~= 2 || exist(tube_file,'file') ~= 2, error('Run Stages 3 and 4 before Stage 5.'); end

S4 = load(nominal_file);
S3 = load(tube_file,'Model','Tube','UncSets');

Params = S4.Params; Reference = S4.Reference; Prep = S4.Prep;
Model = S4.Model; Tube = S4.Tube;
Problem = S4.Problem; Initial = S4.Initial; Optimal = S4.Optimal;
UncSets = S3.UncSets;

%% 4. Shared scenario
Scenario = generate_simulation_scenario(Params);

%% 5. Proposed method
fprintf('\nRunning Tube method...\n');
timer = tic();
Data_prop = run_ilmpc(Problem,Initial,Optimal,Scenario, Model,Tube,UncSets,Params);
total_time_prop = toc(timer);

if isfield(S4,'Nominal')
    expected = S4.Nominal.T_vectors(:,2:end);
    actual = reshape(Data_prop.T_hist,Problem.nt,Params.N_batches);
    discrepancy = norm(actual-expected,'fro') / max(1,norm(expected,'fro'));
    fprintf('Stage 4 / online nominal-update discrepancy: %.3e.\n', discrepancy);
    if discrepancy > 1e-7, error('Stage 5: nominal update differs from Stage 4.'); end
end

checkpoint_file = fullfile(out_dir,'stage05_tube_checkpoint.mat');
save('-mat7-binary',checkpoint_file, 'Data_prop','Scenario','Params','Model','Tube','Optimal');

%% 6. Comparison method (Xue)
fprintf('\nRunning adapted Xue comparison...\n');
timer = tic();
Data_xue = run_xue_ilmpc_v2(Problem,Initial,Optimal,Scenario, Model,Tube,UncSets,Params);
total_time_xue = toc(timer);

%% 7. Comparison method (Oh)
fprintf('\nRunning adapted Oh comparison [Ref 12]...\n');
timer = tic();
Data_oh = run_oh_ilmpc(Problem,Initial,Optimal,Scenario, Model,Tube,UncSets,Params);
total_time_oh = toc(timer);

%% 8. Metrics
Metrics_prop = compute_metrics(Data_prop,Params);
Metrics_xue = compute_metrics(Data_xue,Params);

if exist('Data_oh', 'var') && ~isempty(Data_oh)
    Metrics_oh = compute_metrics(Data_oh,Params);
    Metrics_oh.total_simulation_time = total_time_oh;
else
    Metrics_oh = [];
end

Metrics_prop.total_simulation_time = total_time_prop;
Metrics_xue.total_simulation_time = total_time_xue;

report_simulation_comparison(Metrics_prop,Metrics_xue,Metrics_oh);

%% 9. Save complete results
Artifact = struct('stage',5,'code_version',Params.code_version,'octave_version',OCTAVE_VERSION, ...
    'created_at',datestr(now,31),'state_convention','observed_score_history','shared_scenario',true,'ready_for_reporting',true);

Validation = []; Val_data = []; CV = [];
identification_file = fullfile(out_dir,'stage02_identification.mat');
if exist(identification_file,'file') == 2
    S2 = load(identification_file,'Validation','Val_data','CV');
    if isfield(S2,'Validation'), Validation = S2.Validation; end
    if isfield(S2,'Val_data'), Val_data = S2.Val_data; end
    if isfield(S2,'CV'), CV = S2.CV; end
end

output_file = fullfile(out_dir,'stage05_simulation.mat');
save('-mat7-binary',output_file, 'Params','Reference','Prep','Model','Tube','UncSets', ...
    'Problem','Initial','Optimal','Scenario', 'Data_prop','Data_xue','Data_oh', ...
    'Metrics_prop','Metrics_xue','Metrics_oh', 'Validation','Val_data','CV','Artifact');

fprintf('\nSaved: %s\n',output_file);
fprintf('Stage 5 completed.\n');
diary off;
