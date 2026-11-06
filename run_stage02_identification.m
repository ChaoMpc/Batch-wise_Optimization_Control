% =========================================================================
% run_stage02_identification
% =========================================================================
clear; clc; close all;

if ~exist('OCTAVE_VERSION','builtin')
    error('This simulation project targets GNU Octave.');
end

%% 1. Project directory and Setup Diary
script_file = mfilename('fullpath');
if isempty(script_file), project_dir = pwd; else, project_dir = fileparts(script_file); end
addpath(project_dir); rehash();

out_dir = fullfile(project_dir,'results');
if exist(out_dir,'dir') ~= 7
    [ok,msg] = mkdir(out_dir);
    if ~ok, error('Cannot create result directory: %s',msg); end
end

log_file = fullfile(out_dir, 'stage02_log.txt');
if exist(log_file, 'file'), delete(log_file); end
diary(log_file);

%% 2. Dependencies
required = { ...
    'init_params', 'cryst_plant', 'cryst_measure', 'generate_reference', ...
    'generate_experiment_data', 'generate_excitation', 'fit_preprocessing', ...
    'apply_preprocessing', 'fit_pls', 'project_scores', 'fit_latent_arx', ...
    'realize_latent_model', 'build_identified_model', 'validate_model', ...
    'select_hyperparameters', 'report_identification'};

fprintf('=== Stage 2: identification and cross-validation ===\n');
for j = 1:numel(required)
    local_file = fullfile(project_dir,[required{j},'.m']);
    if exist(local_file,'file') ~= 2, error('Missing required file: %s',local_file); end
    resolved = which(required{j});
    [resolved_dir,~,~] = fileparts(resolved);
    if ~strcmp(resolved_dir,project_dir), error('Function %s resolves outside this project: %s', required{j},resolved); end
end

%% 3. Load Stage 1
stage01_file = fullfile(out_dir,'stage01_setup.mat');
if exist(stage01_file,'file') ~= 2, error('Run run_stage01_setup.m before Stage 2.'); end
S = load(stage01_file);
Params = S.Params; Reference = S.Reference;
if ~isfield(Params,'code_version') || ~strcmp(Params.code_version,'cryst_modular_v1'), error('Stage 1 artifact has an incompatible code version.'); end

%% 4. Generate independent datasets
Train_data = generate_experiment_data(Params,Params.N_train_batches,Params.rng_seed,'training');
Val_data = generate_experiment_data(Params,Params.N_val_batches,Params.rng_seed_val,'selection');
Cal_data = generate_experiment_data(Params,Params.N_unc_batches,Params.rng_seed_unc,'calibration');

%% 5. Fit and apply training preprocessing
Prep = fit_preprocessing(Train_data,Params);
Standard_train = apply_preprocessing(Train_data,Prep,Params);
Standard_val = apply_preprocessing(Val_data,Prep,Params);
Standard_cal = apply_preprocessing(Cal_data,Prep,Params);

%% 6. Cross-validation and final selected model
[Model,CV] = select_hyperparameters(Prep,Standard_train,Standard_val,Val_data,Params);
Params.A_lv = Model.A_lv; Params.n_arx = Model.n_arx;

Validation = validate_model(Model,Standard_val,Val_data,Params);
Calibration = validate_model(Model,Standard_cal,Cal_data,Params);
report_identification(Model,Prep,Validation,Params);

%% 7. Realization consistency check
n = Model.n_arx; A = Model.A_lv; block = Model.n_state_per_lv;
Scores_train = project_scores(Standard_train,Model);
j = Params.ctrl_start+n;
Tb = Scores_train.T(:,:,1); Ub = Scores_train.U(:,:,1);

x_check = zeros(Model.nx,1); phi_check = zeros(2*A*n,1);
for a = 1:A
    idx = (a-1)*block+(1:block);
    x_check(idx(1:n)) = Ub(j:-1:j-n+1,a);
    if n > 1, x_check(idx(n+1:end)) = Tb(j-1:-1:j-n+1,a); end
    iu = (a-1)*n+(1:n); it = A*n+iu;
    phi_check(iu) = Ub(j:-1:j-n+1,a); phi_check(it) = Tb(j:-1:j-n+1,a);
end
input_check = Tb(j,:)';

u_arx = Model.Theta_matrix'*phi_check;
u_ss = Model.barC * (Model.barA*x_check+Model.barB*input_check);
realization_error = norm(u_arx-u_ss,inf);

if realization_error > 1e-10*max(1,norm(u_arx,inf))
    error('Stage 2: ARX/state-space realization mismatch %.3e.', realization_error);
end
fprintf('State-space realization consistency error: %.3e\n', realization_error);

%% 8. Save artifact
Artifact = struct('stage',2,'code_version',Params.code_version,'octave_version',OCTAVE_VERSION, ...
    'created_at',datestr(now,31),'realization_error',realization_error);

stage02_file = fullfile(out_dir,'stage02_identification.mat');
save('-mat7-binary',stage02_file, 'Params','Reference','Train_data','Val_data','Cal_data', ...
    'Prep','Model','CV','Validation','Calibration','Artifact');
fprintf('\nSaved: %s\n',stage02_file);

%% 9. Plot
if Params.plot_results
    columns = Validation.output_columns; time = Validation.time;
    figure('Name','Stage02 independent free-run validation', 'Color','w','Position',[100,100,1200,700]);
    for out = 1:Params.n_out
        subplot(2,3,out); hold on; grid on;
        measured = reshape(Val_data.Y(out,columns,1),1,[]);
        predicted = reshape(Validation.Y_free(out,:,1),1,[]);
        plot(time,measured,'k-','LineWidth',1.2);
        plot(time,predicted,'r--','LineWidth',1.4);
        xlabel('Time (min)'); ylabel(Params.output_names{out},'Interpreter','none');
        title(sprintf('%s: NRMSE %.3f, VAF %.1f%%', Params.output_names{out}, Validation.NRMSE(out),Validation.VAF(out)), 'Interpreter','none');
        if out == 1, legend('Measured','Free-run','Location','northeast'); end
        xlim([Params.ctrl_start*Params.Ts, Params.Ble*Params.Ts]);
    end
end
fprintf('\nStage 2 completed.\n');
diary off;
