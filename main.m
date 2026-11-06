% ========================================================================
% main.m -- DyPLS-IL-MPC closed-loop simulation (Wood-Berry)
%
% Updated for the N-step lifting:
%   * y(1..N) is predicted (N blocks); v(0..N-1) is predicted (N blocks).
%   * gen_validation_data uses the training statistics (Ident_data) so
%     validation/test X0/Y0 are in the SAME standardized coordinates as
%     the training data.
%   * free_run_validation now returns 5 outputs (incl. 1-step metrics).
%   * All closed-loop data are DEVIATION quantities; plotting adds
%     y_ss/u_ss back.
% ========================================================================

clear; clc; close all;

if exist('OCTAVE_VERSION', 'builtin')
    pkg load control optim;
end

fprintf('=== DyPLS-IL-MPC Simulation (Wood-Berry) ===\n');

% ---------- Step 1: Parameters ----------
Params = init_params();
if isfield(Params, 'rng_seed') && ~isempty(Params.rng_seed)
    if exist('rng', 'builtin'), rng(Params.rng_seed);
    else, rand('seed', Params.rng_seed); randn('seed', Params.rng_seed); end
end

% ---------- Step 2: Identification data ----------
[Ident_data, Params] = gen_prbs_excitation(Params);
fprintf('main: identification data generated (N=%d).\n', size(Ident_data.X0,1));

% ---------- Step 3: Validation data (uses TRAINING statistics) ----------
[Val_data, Params] = gen_validation_data(Params, 'val', Ident_data);
fprintf('main: validation data generated (N_val=%d).\n', size(Val_data.X0,1));

% ---------- Step 4: Modeling ----------
Model = build_dypls(Ident_data, Params);
fprintf('main: DyPLS-ARX model built (A=%d, n_arx=%d).\n', Params.A_lv, Params.n_arx);

diag_arx_transient(Ident_data, Val_data, Params);

% Model 是用 Ident_data 训练好的最终模型
results_test = evaluate_on_test_set(Model, Params, Ident_data, 10);

% ---------- Step 5: Free-run validation ----------
fprintf('\n--- Free-run validation ---\n');
[Y0_sim, nrmses, vafs, nrmses_1step, vafs_1step] = ...
    free_run_validation(Model, Val_data, Params);
fprintf('  Global Validation NRMSE (free-run): [%s]\n', num2str(nrmses, '%.4f '));
fprintf('  Global Validation VAF   (free-run): [%s] %%\n', num2str(vafs, '%.4f '));
fprintf('  Global Validation NRMSE (1-step)  : [%s]\n', num2str(nrmses_1step, '%.4f '));
fprintf('  Global Validation VAF   (1-step)  : [%s] %%\n', num2str(vafs_1step, '%.4f '));

start_idx = 1; window_len = 60;
[Y0_sim_win, nrmses_win, vafs_win, nrmses_win_1step, vafs_win_1step] = ...
    free_run_validation(Model, Val_data, Params, start_idx, window_len);
ValResults.Y0_sim         = Y0_sim_win;        % windowed, standardized deviation
ValResults.nrmses         = nrmses_win;
ValResults.vafs           = vafs_win;
ValResults.nrmses_1step   = nrmses_win_1step;
ValResults.vafs_1step     = vafs_win_1step;
ValResults.Val_data       = Val_data;
ValResults.start_idx      = start_idx;
ValResults.window_len     = window_len;

% ---------- Step 6: Uncertainty sets ----------
UncSets = estimate_uncertainty(Model, Ident_data, Params);
fprintf('main: uncertainty sets estimated.\n');

% ---------- Step 7: Tube parameters ----------
Tube = compute_tube(Ident_data, Model, UncSets, Params);
fprintf('main: Tube-MPC parameters computed. Spectral radius = %.4f\n', Tube.spectral_radius);

diag_tube_limits(Ident_data, Model, Tube, Params);

opts = struct();
opts.alpha_list = [0.10, 0.1077, 0.15];       % 3 个 alpha
opts.beta_grid  = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0];
opts.latex_dir  = 'figures';
opts.tag        = '';
scan_tube_table(Ident_data, Model, UncSets, Params, opts);


% ---------- Step 7.5: Report ----------
fprintf('\n=== Uncertainty and Invariant Set Parameters ===\n');
fprintf('[Uncertainty Sets]\n');
fprintf('  eta_bound (ARX innovation): [%s]\n', num2str(UncSets.eta_bound(:)', '%.4e '));
if norm(UncSets.W_r, 'fro') == 0
    fprintf('  r_bound   (PLS residual)  : 0 (A_lv == n_out, no truncation)\n');
else
    fprintf('  r_bound   (PLS residual)  : [%s]\n', num2str(UncSets.r_bound(:)', '%.4e '));
end

fprintf('\n[Ellipsoid Shape Matrices (Diagonal Elements)]\n');
fprintf('  W_eta (diag) : [%s]\n', num2str((diag(Tube.W_eta))', '%.4e '));
if norm(Tube.W_r, 'fro') > 0
    fprintf('  W_r   (diag) : [%s]\n', num2str((diag(Tube.W_r))', '%.4e '));
else
    fprintf('  W_r   (diag) : 0 (A_lv == n_out, no truncation)\n');
end

eig_P = eig(Tube.P_e);
fprintf('\n[Invariant Ellipsoid Matrix P_e]\n');
fprintf('  Dimension         : %d x %d\n', size(Tube.P_e,1), size(Tube.P_e,2));
fprintf('  min(eig(P_e))     : %.4e\n', min(eig_P));
fprintf('  max(eig(P_e))     : %.4e\n', max(eig_P));
fprintf('  Condition Number  : %.4e\n', max(eig_P)/min(eig_P));
fprintf('  trace(P_e)        : %.4e\n', trace(Tube.P_e));

fprintf('\n[Tightening Margins]\n');
fprintf('  delta_v : [%s]\n', num2str(Tube.delta_v(:)', '%.4e '));
fprintf('  delta_y : [%s]\n', num2str(Tube.delta_y(:)', '%.4e '));

fprintf('\n[Tuning Parameters]\n');
fprintf('  alpha (Tube) : %.4f (adaptive)\n', Tube.alpha);
if isfield(Tube, 'beta')
    fprintf('  beta  (Ellip): %.4f\n', Tube.beta);
else
    fprintf('  beta  (Ellip): 1.0\n');
end
fprintf('================================================\n\n');

if ~exist('figures','dir'), mkdir('figures'); end
fid = fopen(fullfile('figures', 'uncertainty_params.txt'), 'w');
fprintf(fid, 'eta_bound = %s\n', num2str(UncSets.eta_bound(:)', '%.6e '));
if norm(UncSets.W_r, 'fro') == 0
    fprintf(fid, 'r_bound   = 0 (A_lv == n_out, no truncation)\n');
else
    fprintf(fid, 'r_bound   = %s\n', num2str(UncSets.r_bound(:)', '%.6e '));
end
fprintf(fid, 'W_eta_diag = %s\n', num2str((diag(Tube.W_eta))', '%.6e '));
if norm(Tube.W_r, 'fro') > 0
    fprintf(fid, 'W_r_diag   = %s\n', num2str((diag(Tube.W_r))', '%.6e '));
else
    fprintf(fid, 'W_r_diag   = 0 (A_lv==n_out)\n');
end
fprintf(fid, 'P_e_eigs   = %s\n', num2str(eig_P', '%.6e '));
fprintf(fid, 'delta_v    = %s\n', num2str(Tube.delta_v(:)', '%.6e '));
fprintf(fid, 'delta_y    = %s\n', num2str(Tube.delta_y(:)', '%.6e '));
fclose(fid);

% ---------- Step 8: Unreachability verification ----------
fprintf('\n--- Direct verification of unreachable reference (Problem 55) ---\n');
VerifyResult = verify_unreachability(Ident_data, Model, Tube, Params);
fprintf('  y*_s(k_i) at Kp = [%s]  (wt%% dev)\n', ...
    num2str(VerifyResult.y_s_at_Kp(:,1)', '%.4f '));
fprintf('  Reference     at Kp = [%s]  (wt%% dev)\n', ...
    num2str(Params.Ref_points(:)', '%.4f '));
fprintf('  Nominal error at y*_s : [%s]\n', ...
    num2str(VerifyResult.nom_error_at_Kp(:)', '%.4f '));
fprintf('  Active constraints: %d\n', numel(VerifyResult.active));

% ---------- Step 9: Closed-loop simulation ----------
Data = run_ilmpc(Ident_data, Model, Tube, UncSets, Params);

% ---------- Step 10: Metrics ----------
metrics = compute_metrics(Data, Params, Model, UncSets);

fprintf('\n=== Closed-loop metrics ===\n');
fprintf('  Constraint violations (input): %d (max %+.2e)\n', ...
    metrics.n_violations_v, metrics.max_violation_v);
fprintf('  Constraint violations (output): %d (max %+.2e)\n', ...
    metrics.n_violations_y, metrics.max_violation_y);
if isnan(metrics.coverage_eta)
    fprintf('  eta coverage: N/A (no samples)\n');
else
    fprintf('  eta coverage: %.2f%% (%d samples)\n', ...
        100*metrics.coverage_eta, metrics.n_eta_samples);
end
if isnan(metrics.coverage_r)
    fprintf('  r   coverage: N/A (A_lv == n_out, r == 0)\n');
else
    fprintf('  r   coverage: %.2f%% (%d samples)\n', ...
        100*metrics.coverage_r, metrics.n_r_samples);
end
fprintf('  Mean QP solve time: %.4f s\n', mean(metrics.qp_times));

% --- Report J_batch (matches Theorem 2: J^{b,*} = V^b + sum||dT||_R^2) ---
if isfield(Data, 'J_batch') && ~isempty(Data.J_batch)
    fprintf('\n  J_batch (first, final) = (%.4e, %.4e)\n', ...
        Data.J_batch(1), Data.J_batch(end));
    % Monotonicity check for V^b (paper: V^{b+1} + sum||dT^{b+1}||_R^2 <= V^b)
    dJ = diff(Data.J_batch);
    if all(dJ <= 1e-8 * max(1, abs(Data.J_batch(1))))
        fprintf('  J_batch is non-increasing (consistent with Theorem 2).\n');
    else
        fprintf('  WARNING: J_batch not strictly non-increasing; check feasibility/tube.\n');
    end
end

% ---------- Step 11: Initial-state test ----------
if isfield(Params, 'run_init_state_test') && Params.run_init_state_test
    fprintf('\n=== Test: Multiple initial-state perturbations within Z_e,inf ===\n');

    % 1. 先计算标称状态下的数据
    seed_test = Params.rng_seed + 1000;
    if exist('rng', 'builtin'), rng(seed_test); else, randn('seed', seed_test); end
    Params_temp_nom = Params;
    Params_temp_nom.init_state_variation = 0.0;
    Data_nom = run_ilmpc(Ident_data, Model, Tube, UncSets, Params_temp_nom);

    % 2. 准备边界扰动向量 (利用 Cholesky 分解 P_e = L*L')
    % 这样生成的扰动天然落在椭球 E(P_e) 的边界上
    L_Pe = chol(Tube.P_e, 'lower');
    z_rand = randn(size(Tube.P_e, 1), 1);
    z_rand = z_rand / sqrt(z_rand' * z_rand); % 归一化，使其范数为1
    e_boundary = L_Pe * z_rand; % e^T P_e^-1 e = 1

    % 3. 计算不同缩放比例下的扰动向量
    e_small = 0.3 * e_boundary; % 内部扰动 (二次型 = 0.09)
    e_mid   = 0.7 * e_boundary; % 中部扰动 (二次型 = 0.49)
    e_edge  = 0.99 * e_boundary;% 边界扰动 (二次型 = 0.9801)

    % 打印给审稿人看的数据
    fprintf('  Perturbation sizes (e(0)^T * P_e^(-1) * e(0)):\n');
    fprintf('    Small  : %.4f (e_norm = %.4e)\n', e_small'*(Tube.P_e\e_small), norm(e_small));
    fprintf('    Mid    : %.4f (e_norm = %.4e)\n', e_mid'*(Tube.P_e\e_mid), norm(e_mid));
    fprintf('    Edge   : %.4f (e_norm = %.4e)\n', e_edge'*(Tube.P_e\e_edge), norm(e_edge));

    % 4. 跑测试
    Params_temp_pert = Params;
    Params_temp_pert.init_state_vector = e_small;
    Data_pert_small = run_ilmpc(Ident_data, Model, Tube, UncSets, Params_temp_pert);

    Params_temp_pert.init_state_vector = e_mid;
    Data_pert_mid = run_ilmpc(Ident_data, Model, Tube, UncSets, Params_temp_pert);

    Params_temp_pert.init_state_vector = e_edge;
    Data_pert_edge = run_ilmpc(Ident_data, Model, Tube, UncSets, Params_temp_pert);

    Data_init_test.Nominal = Data_nom;
    Data_init_test.Small   = Data_pert_small;
    Data_init_test.Mid     = Data_pert_mid;
    Data_init_test.Edge    = Data_pert_edge;
    Data.init_test = Data_init_test;

    summary_metrics.pert_small = compute_metrics(Data_pert_small, Params, Model, UncSets);
    summary_metrics.pert_mid   = compute_metrics(Data_pert_mid, Params, Model, UncSets);
    summary_metrics.pert_edge  = compute_metrics(Data_pert_edge, Params, Model, UncSets);
end

fprintf('\n--- Summary for Validation Table ---\n');
summary_metrics.main = metrics;
if isfield(Data, 'init_test')
    % 适配新的多组扰动测试结构
    summary_metrics.pert_small = compute_metrics(Data.init_test.Small, Params, Model, UncSets);
    summary_metrics.pert_mid   = compute_metrics(Data.init_test.Mid, Params, Model, UncSets);
    summary_metrics.pert_edge  = compute_metrics(Data.init_test.Edge, Params, Model, UncSets);
end

% ---------- Step 11.5: Monte Carlo Statistics for Reviewer ----------
%if isfield(Params, 'run_mc_test') && Params.run_mc_test
%    run_mc_initial_states(Ident_data, Model, Tube, UncSets, Params);
%end
% ---------- Step 12: Plot ----------
if isfield(Params, 'plot_results') && Params.plot_results
    plot_results(Data, Params, Model, Tube, VerifyResult, metrics, ValResults, summary_metrics);
    fprintf('main: plotting completed.\n');
end

% ---------- Step 13: Export ----------
if isfield(Params, 'plot_results') && Params.plot_results
    try
        export_all_figures(Params, 'figures');
        fprintf('main: figures exported to figures/\n');
    catch ME
        fprintf(2, '[export_all_figures] failed: %s\n', ME.message);
    end
end

fprintf('\n=== Simulation completed ===\n');
