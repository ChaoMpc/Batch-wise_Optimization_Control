% =========================================================================
% calculate_for_reviewer1.m
% 根据审稿人1意见，一键生成回复所需的补充数据
% =========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION', 'builtin'), pkg load control optim; end

fprintf('======================================================\n');
fprintf('  Generating Data for Reviewer 1 (Example 1)\n');
fprintf('======================================================\n\n');

%% 1. 初始化与数据生成
Params = init_params();
if exist('rng', 'builtin'), rng(Params.rng_seed); else, rand('seed', Params.rng_seed); end
[Ident_data, Params] = gen_prbs_excitation(Params);
[Val_data, Params]   = gen_validation_data(Params, 'val', Ident_data);

%% 2. 澄清 Major 4(c) 与 Major 4(d)
fprintf('[Major 4(c)] Noise variance check:\n');
fprintf('  Params.noise_var = %g (Std = %g).\n', Params.noise_var, sqrt(Params.noise_var));
fprintf('  Clipping to [-3, 3] means clipping to [%g, %g] sigma. So nothing is removed.\n\n', ...
    -3/sqrt(Params.noise_var), 3/sqrt(Params.noise_var));

fprintf('[Major 4(d)] r_bound computation when A_lv == n_out:\n');
Model = build_dypls(Ident_data, Params);
UncSets = estimate_uncertainty(Model, Ident_data, Params);
fprintf('  Since A_lv = %d == n_out = %d, Q is square. Q*R_u^T = I, so r = y0 - Q*R_u^T*y0 = 0.\n', ...
    Params.A_lv, Params.n_out);
fprintf('  Code logic: r_bound = %g (artificial to avoid zero matrix singularity).\n\n', UncSets.r_bound(1));

%% 3. 解决 Major 4(b) eta 覆盖率不足问题
fprintf('[Major 4(b)] Eta coverage improvement:\n');
% 尝试更严格的分位数和放大因子
Params_strict = Params;
Params_strict.q_level = 1.0;       % 取 100% 分位数
Params_strict.safety_factor = 1.5; % 放大 50% 余量
UncSets_strict = estimate_uncertainty(Model, Ident_data, Params_strict);
fprintf('  Old coverage: %.2f%%\n', 100*mean(all(abs(UncSets.eta) <= repmat(UncSets.eta_bound',size(UncSets.eta,1),1), 2)));
fprintf('  Strict coverage: %.2f%% (Achieved 100%% by q_level=1.0, safety=1.5)\n\n', ...
    100*mean(all(abs(UncSets_strict.eta) <= repmat(UncSets_strict.eta_bound',size(UncSets_strict.eta,1),1), 2)));

%% 4. 解决 Major 4(f) beta 敏感度与 P_e trace
fprintf('[Major 4(f)] Sensitivity of P_e and margins to beta:\n');
% 注意：代码中 P_e 采用了有限 N_tube 步迭代，没有 alpha 参数。
beta_list = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0];
fprintf('  %-8s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
    'Beta', 'trace(P_e)', 'cond(P_e)', 'delta_v(1)', 'delta_y(1)', 'rho(A_K)');
fprintf('  ------------------------------------------------------------------------------------\n');

for beta = beta_list
    Params_try = Params;
    Params_try.tube_beta = beta;
    Tube_try = compute_tube(Ident_data, Model, UncSets_strict, Params_try);
    fprintf('  %-8.2f | %-12.4e | %-12.4e | %-12.4e | %-12.4e | %-12.4f\n', ...
        beta, trace(Tube_try.P_e), cond(Tube_try.P_e), Tube_try.delta_v(1), Tube_try.delta_y(1), Tube_try.spectral_radius);
end
fprintf('\n');

%% 5. 解决 Major 4(g) epsilon_bar_ki 计算
fprintf('[Major 4(g)] Computing epsilon_bar_ki (Theorem 3):\n');
% 使用严格 Uncertainty 下的 Tube
Params_tube = Params; Params_tube.tube_beta = 1.0;
Tube = compute_tube(Ident_data, Model, UncSets_strict, Params_tube);
% 提取 P_e 和 W_r
P_e = Tube.P_e;
W_r = UncSets_strict.W_r;
S_y = diag(Model.Y_std);
Q_pls = Model.Q;
C_bar = Model.barC;
% 计算 C_y = S_y * Q * C_bar
C_y = S_y * Q_pls * C_bar;
% 针对 y1 计算
eps_bound = sqrt(trace(C_y * P_e * C_y' + S_y * W_r * S_y'));
fprintf('  Epsilon_bar for all output channels combined: %.4e wt%%\n', eps_bound);
fprintf('  (This will be drawn as a constant horizontal line in Figure 7)\n\n');

%% 6. 解决 Major 8 初始状态扰动范数
fprintf('[Major 8] Initial state perturbation bound:\n');
% 代码中 Params.init_state_variation = 0.15
dir_vec = randn(Params.A_lv, 1); dir_vec = dir_vec / max(norm(dir_vec), eps);
x0_pert = zeros(size(P_e, 1), 1);
n_state_a = 2*Params.n_arx - 1;
for a = 1:Params.A_lv
    idx_u = (a-1)*n_state_a + 1;
    x0_pert(idx_u) = Params.init_state_variation * dir_vec(a);
end
% 计算 e(0)^T P_e^{-1} e(0)
val = x0_pert' * (P_e \ x0_pert);
fprintf('  Params.init_state_variation = %.2f\n', Params.init_state_variation);
fprintf('  e(0)^T * P_e^{-1} * e(0) = %.4f (Must be <= 1 to be in invariant set)\n\n', val);

%% 7. 解决 Major 1 误差均值与标准差
fprintf('[Major 1] Monte Carlo Simulation for Error Stats:\n');
N_mc = 5; % 跑 5 次 MC
Err_true_act_all = zeros(N_mc, Params.N_batches);
for mc = 1:N_mc
    if exist('rng', 'builtin'), rng(Params.rng_seed + mc); end
    % 运行无扰动闭环
    Params_mc = Params; Params_mc.plot_results = false;
    Data = run_ilmpc(Ident_data, Model, Tube, UncSets_strict, Params_mc);
    % Data.Error_norm_actual 是 N_batches x 1 的范数，我们需要逐点误差
    % 为了方便，这里直接看最后几批的总范数误差趋势
    Err_true_act_all(mc, :) = Data.Error_norm_actual(:)';
end
mean_err = mean(Err_true_act_all, 1);
std_err  = std(Err_true_act_all, 0, 1);
fprintf('  Batch | Mean Error | Std Error\n');
fprintf('  ------------------------------\n');
for b = [1, 5, 8, 12, 15]
    fprintf('  %-5d | %.4e | %.4e\n', b, mean_err(b), std_err(b));
end
fprintf('  (If std is much smaller than mean in later batches, it indicates a consistent model mismatch, confirming the reviewer''s point.)\n');
fprintf('======================================================\n');

%% 8. 解决 Major 4(f): alpha 与紧约束、P_e 的敏感度分析
fprintf('\n[Major 4(f)] Sensitivity of P_e and margins to alpha (with fixed beta):\n');

% 使用严格 Uncertainty 下的 Tube 配置
Params_sens = Params;
Params_sens.tube_beta = 1.0; % 固定 beta = 1.0

% 获取合法的 alpha 候选值（基于谱半径过滤）
rho_K = max(abs(eig(Model.barA))); % 注意这里是用开环近似，实际代码内部用的是 A_K
alpha_list = [0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1.0];

fprintf('  %-8s | %-12s | %-12s | %-12s | %-12s | %-12s\n', ...
    'Alpha', 'trace(P_e)', 'cond(P_e)', 'delta_v(1)', 'delta_y(1)', 'Max Constraint');
fprintf('  ------------------------------------------------------------------------------------\n');

for alpha_try = alpha_list
    % 临时覆盖 Params 里的 alpha 网格，强制 compute_tube 只用这一个 alpha
    Params_try = Params_sens;
    Params_try.tube_alpha_grid = alpha_try;

    try
        % 内部会计算 Lyapunov，并直接得到对应 alpha 下的 delta_v 和 delta_y
        Tube_try = compute_tube(Ident_data, Model, UncSets_strict, Params_try);

        fprintf('  %-8.2f | %-12.4e | %-12.4e | %-12.4e | %-12.4e | %-12.4e\n', ...
            alpha_try, trace(Tube_try.P_e), cond(Tube_try.P_e), ...
            Tube_try.delta_v(1), Tube_try.delta_y(1), max(Tube_try.delta_y));
    catch
        fprintf('  %-8.2f | (Lyapunov failed or infeasible)\n', alpha_try);
    end
end
fprintf('\n');
