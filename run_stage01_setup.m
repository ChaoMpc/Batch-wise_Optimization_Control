% =========================================================================
% run_stage01_setup
%
% Stage 1:
%   initialize parameters;
%   generate nonlinear baseline;
%   check measurement/integration conventions;
%   save results/stage01_setup.mat;
%   plot baseline inputs and outputs.
% =========================================================================

clear;
clc;
close all;

if ~exist('OCTAVE_VERSION','builtin')
    error('This simulation project targets GNU Octave.');
end

%% 1. Resolve project directory & Setup Diary
script_file = mfilename('fullpath');
if isempty(script_file)
    project_dir = pwd;
else
    project_dir = fileparts(script_file);
end
addpath(project_dir);
rehash();

out_dir = fullfile(project_dir,'results');
if exist(out_dir,'dir') ~= 7
    [ok,msg] = mkdir(out_dir);
    if ~ok, error('Cannot create result directory: %s',msg); end
end

log_file = fullfile(out_dir, 'stage01_log.txt');
if exist(log_file, 'file'), delete(log_file); end
diary(log_file);

%% 2. Check all dependencies in this batch
required = { ...
    'init_params', ...
    'cryst_plant', ...
    'cryst_measure', ...
    'generate_reference' ...
};

fprintf('=== Stage 1: nonlinear plant and baseline ===\n');
fprintf('Octave version: %s\n',OCTAVE_VERSION);
fprintf('Project directory: %s\n',project_dir);

for j = 1:numel(required)
    filename = fullfile(project_dir,[required{j},'.m']);
    if exist(filename,'file') ~= 2
        error('Missing required file: %s',filename);
    end
    resolved_file = which(required{j});
    fprintf('  %s -> %s\n',required{j},resolved_file);
    [resolved_dir,~,~] = fileparts(resolved_file);
    if ~strcmp(resolved_dir,project_dir)
        error(['Function %s resolves outside this project. ', ...
               'Remove conflicting files from the Octave path.'], ...
              required{j});
    end
end

%% 3. Parameters and baseline
Params = init_params();
[Reference,Params] = generate_reference(Params);

%% 4. Check the measurement-only interface
[y0,x_same] = cryst_plant(Params.x0_nom,[],Params,false);

if norm(x_same-Params.x0_nom,inf) > 1e-14
    error('Stage 1: measurement-only mode advanced the state.');
end
if norm(y0-Reference.Y(:,1),inf) > 1e-12
    error('Stage 1: initial measurement is inconsistent.');
end

%% 5. Check one deterministic integration step
[y1,x1] = cryst_plant(Params.x0_nom,Reference.V(:,1),Params,false);

if norm(x1-Reference.X(:,2),inf) > 1e-12 || ...
   norm(y1-Reference.Y(:,2),inf) > 1e-12
    error('Stage 1: baseline first-step integration is inconsistent.');
end

%% 6. Check explicit noise injection without random-number dependence
noise_example = 0.5*Params.measurement_noise_halfwidth;
[y_noise,y_true] = cryst_measure(Params.x0_nom,Params,true,noise_example);
expected_noisy = max(0,y_true+noise_example);

if norm(y_noise-expected_noisy,inf) > 1e-12
    error('Stage 1: explicit noise injection is inconsistent.');
end

%% 7. Print baseline diagnostics
fprintf('\nBaseline output ranges:\n');
for j = 1:Params.n_out
    fprintf('  %-6s: min=%12.6g, max=%12.6g\n', ...
        Params.output_names{j}, Reference.output_min(j), Reference.output_max(j));
end
fprintf('Maximum baseline output violation: %.6e\n', Reference.max_output_violation);

fprintf('\nTime indexing:\n');
fprintf('  V(:,61) acts from 60 to 61 min.\n');
fprintf('  Y(:,61) is the output at 60 min.\n');
fprintf('  Y(:,62) is the output at 61 min.\n');

%% 8. Save stage artifact
Artifact = struct();
Artifact.stage = 1;
Artifact.code_version = Params.code_version;
Artifact.octave_version = OCTAVE_VERSION;
Artifact.created_at = datestr(now,31);

out_file = fullfile(out_dir,'stage01_setup.mat');
save('-mat7-binary',out_file, 'Params','Reference','Artifact');
fprintf('\nSaved: %s\n',out_file);

%% 9. Plot baseline
if Params.plot_results
    figure('Name','Stage01 baseline inputs', 'Color','w','Position',[100,100,1200,650]);
    for j = 1:Params.n_in
        subplot(2,3,j); hold on; grid on;
        stairs(Reference.time_v,Reference.V(j,:), 'b-','LineWidth',1.4);
        plot([0,Params.Ble*Params.Ts], Params.u_phys_min(j)*[1,1],'k--');
        plot([0,Params.Ble*Params.Ts], Params.u_phys_max(j)*[1,1],'k--');
        xlabel('Time (min)'); ylabel(Params.input_names{j},'Interpreter','none');
        title(Params.input_names{j},'Interpreter','none'); xlim([0,Params.Ble*Params.Ts]);
    end

    figure('Name','Stage01 baseline outputs', 'Color','w','Position',[130,130,1200,650]);
    for j = 1:Params.n_out
        subplot(2,3,j); hold on; grid on;
        plot(Reference.time_y,Reference.Y(j,:), 'b-','LineWidth',1.5);
        if j == Params.track_idx
            plot([0,Params.Ble*Params.Ts], Params.y_phys_min(j)*[1,1],'k--');
            plot([0,Params.Ble*Params.Ts], Params.y_phys_max(j)*[1,1],'k--');
            plot(Reference.tracking_times, Reference.tracking_targets, 'ro','MarkerFaceColor','r','MarkerSize',7);
        end
        yl = ylim();
        plot(Params.ctrl_start*Params.Ts*[1,1], yl,'m:','LineWidth',1.2);
        ylim(yl);
        xlabel('Time (min)'); ylabel(Params.output_names{j},'Interpreter','none');
        title(Params.output_names{j},'Interpreter','none'); xlim([0,Params.Ble*Params.Ts]);
    end
end

fprintf('\nStage 1 completed.\n');
diary off;
