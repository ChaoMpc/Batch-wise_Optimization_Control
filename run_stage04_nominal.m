% =========================================================================
% run_stage04_nominal
% =========================================================================
clear; clc; close all;
if ~exist('OCTAVE_VERSION','builtin'), error('This simulation project targets GNU Octave.'); end

%% 1. Project path & Diary
script_file = mfilename('fullpath');
if isempty(script_file), project_dir = pwd; else, project_dir = fileparts(script_file); end
addpath(project_dir); rehash();

out_dir = fullfile(project_dir,'results');
log_file = fullfile(out_dir, 'stage04_log.txt');
if exist(log_file, 'file'), delete(log_file); end
diary(log_file);

fprintf('=== Stage 4: common nominal optimization ===\n');

%% 2. Check dependencies
required = {'build_lifting_matrices','build_nominal_problem','solve_qp_checked','compute_initial_feasible', ...
    'compute_optimal_reachable','solve_batch_update','evaluate_nominal_trajectory'};
for j = 1:numel(required)
    expected = fullfile(project_dir,[required{j},'.m']);
    if exist(expected,'file') ~= 2, error('Missing required file: %s',expected); end
    resolved = which(required{j});
    [resolved_dir,~,~] = fileparts(resolved);
    if ~strcmp(resolved_dir,project_dir), error('Function %s resolves outside this project: %s', required{j},resolved); end
end
if exist('qp','file') == 0, error('Stage 4: Octave qp is unavailable.'); end

%% 3. Load Stage 3 artifact
input_file = fullfile(out_dir,'stage03_tube.mat');
if exist(input_file,'file') ~= 2, error('Run run_stage03_tube.m before Stage 4.'); end
S = load(input_file);

required_fields = {'Params','Reference','Prep','Model','Tube','Artifact'};
for j = 1:numel(required_fields)
    if ~isfield(S,required_fields{j}), error('Stage 3 artifact is missing %s.',required_fields{j}); end
end
if ~isfield(S.Artifact,'ready_for_nominal_optimization') || ~S.Artifact.ready_for_nominal_optimization
    error('Stage 3 did not certify nominal feasibility.');
end

Params = S.Params; Reference = S.Reference; Prep = S.Prep; Model = S.Model; Tube = S.Tube;
if ~strcmp(Params.code_version,'cryst_modular_v1'), error('Stage 4: incompatible artifact version.'); end

%% 4. Build the common problem
Problem = build_nominal_problem(Model,Tube,Params);
fprintf('Horizon: %d; decision variables: %d; constraints: %d.\n', Problem.H,Problem.nt,size(Problem.A,1));

%% 5. Check lifting against direct state recursion
probe = 0.01*sin((1:Problem.nt)');
ProbeTrajectory = evaluate_nominal_trajectory(Problem,probe);
x = Model.x0_ctrl; max_lifting_error = 0;
Mv = diag(Model.X_std)*Model.P; My = diag(Model.Y_std)*Model.Q*Model.barC;

for k = 1:Problem.H
    j = Params.ctrl_start+k; t_k = ProbeTrajectory.T(:,k);
    max_lifting_error = max(max_lifting_error, norm(x-ProbeTrajectory.X_pre(:,k),inf));
    v = Params.u_base(Model.input_idx,j)+Model.X_mean+Mv*t_k;
    x = Model.barA*x+Model.barB*t_k;
    y = Params.y_base(:,j+1)+Model.Y_mean+My*x;
    max_lifting_error = max([max_lifting_error, norm(x-ProbeTrajectory.X_post(:,k),inf), norm(v-ProbeTrajectory.V(:,k),inf), norm(y-ProbeTrajectory.Y(:,k),inf)]);
end
probe_scale = max([1; abs(ProbeTrajectory.X_post(:)); abs(ProbeTrajectory.V(:)); abs(ProbeTrajectory.Y(:))]);
if max_lifting_error > 1e-9*probe_scale, error('Stage 4: lifting consistency check failed %.3e.', max_lifting_error); end
fprintf('Lifting/direct recursion discrepancy: %.3e.\n', max_lifting_error);

%% 6. Feasible initialization and limiting optimum
Initial = compute_initial_feasible(Problem,Tube,Params);
Optimal = compute_optimal_reachable(Problem,Initial.t,Params);

%% 7. Finite nominal iterations
B = Params.N_batches; nselected = numel(Problem.selected_rows);
Nominal = struct();
Nominal.T_vectors = zeros(Problem.nt,B+1);
Nominal.selected_outputs = zeros(nselected,B+1);
Nominal.tracking_cost = zeros(B+1,1);
Nominal.step_norm = zeros(B,1); Nominal.proximal_cost = zeros(B,1); Nominal.descent_residual = zeros(B,1);
Nominal.QP_times = zeros(B,1); Nominal.max_scaled_violation = zeros(B,1);

t = Initial.t;
Nominal.T_vectors(:,1) = t;
Nominal.selected_outputs(:,1) = Initial.Trajectory.selected_outputs;
Nominal.tracking_cost(1) = Initial.Trajectory.tracking_cost;

for b = 1:B
    Update = solve_batch_update(Problem,t,Params);
    t = Update.t;
    Nominal.T_vectors(:,b+1) = t;
    Nominal.selected_outputs(:,b+1) = Update.Trajectory.selected_outputs;
    Nominal.tracking_cost(b+1) = Update.new_cost;
    Nominal.step_norm(b) = Update.step_norm;
    Nominal.proximal_cost(b) = Update.proximal_cost;
    Nominal.descent_residual(b) = Update.descent_residual;
    Nominal.QP_times(b) = Update.QP.elapsed;
    Nominal.max_scaled_violation(b) = Update.QP.max_scaled_violation;

    output_distance = norm(Update.Trajectory.selected_outputs-Optimal.selected_outputs);
    fprintf('Batch %2d/%2d: J=%.8e; ||d||=%.4e; selected-output distance=%.4e; QP=%.4fs.\n', b,B,Update.new_cost,Update.step_norm, output_distance,Update.QP.elapsed);
end

Nominal.optimal_cost = Optimal.cost;
Nominal.cost_difference = Nominal.tracking_cost-Optimal.cost;
Nominal.output_distance = sqrt(sum(bsxfun(@minus,Nominal.selected_outputs, Optimal.selected_outputs).^2,1))';
Nominal.final_trajectory = evaluate_nominal_trajectory(Problem,t);

comparison_tol = 1e-7*max([Initial.Trajectory.tracking_cost,Optimal.cost,1e-12]);
if any(Nominal.cost_difference < -comparison_tol), error('Stage 4: an iterative feasible point has lower cost than the reported limiting optimum.'); end

%% 8. Save stage artifact
Artifact = struct('stage',4,'code_version',Params.code_version,'octave_version',OCTAVE_VERSION, ...
    'created_at',datestr(now,31),'state_convention','observed_score_history','max_lifting_error',max_lifting_error,'ready_for_online_simulation',true);

output_file = fullfile(out_dir,'stage04_nominal.mat');
save('-mat7-binary',output_file, 'Params','Reference','Prep','Model','Tube','Problem','Initial','Optimal','Nominal','Artifact');
fprintf('\nSaved: %s\n',output_file);

%% 9. Plots
if Params.plot_results
    figure('Name','Stage04 nominal objective', 'Color','w','Position',[100,100,900,450]);
    hold on; grid on;
    plot(0:B,Nominal.tracking_cost,'b-o','LineWidth',1.5);
    plot([0,B],Optimal.cost*[1,1],'k--','LineWidth',1.2);
    xlabel('Number of batch updates'); ylabel('Weighted nominal tracking cost');
    legend('Proximal iteration','Limiting optimum', 'Location','northeast');

    figure('Name','Stage04 nominal selected outputs', 'Color','w','Position',[130,130,1000,500]);
    for j = 1:nselected
        subplot(nselected,1,j); hold on; grid on;
        plot(0:B,Nominal.selected_outputs(j,:), 'b-o','LineWidth',1.3);
        plot([0,B],Optimal.selected_outputs(j)*[1,1], 'g--','LineWidth',1.2);
        plot([0,B],Problem.selected_refs(j)*[1,1], 'r:','LineWidth',1.2);
        xlabel('Number of batch updates'); ylabel(sprintf('Selected output %d',j));
        if j == 1, legend('Nominal iteration','Optimal reachable','True target', 'Location','northeast'); end
    end
end

fprintf('\nStage 4 completed.\n');
fprintf('Initial tracking cost: %.8e.\n',Nominal.tracking_cost(1));
fprintf('Final tracking cost:   %.8e.\n',Nominal.tracking_cost(end));
fprintf('Optimal tracking cost: %.8e.\n',Optimal.cost);
fprintf('Final output distance: %.8e.\n',Nominal.output_distance(end));
diary off;
