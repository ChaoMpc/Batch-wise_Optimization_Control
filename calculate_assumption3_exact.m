% =========================================================================
% calculate_assumption3_exact.m
% 一键计算 Example 1 中 beta_v = 0.15 时的最大二次型
% =========================================================================
clear; clc;
if exist('OCTAVE_VERSION','builtin'), pkg load control optim; end

fprintf('=== Rebuilding Example 1 model and Tube ===\n');
% 1. 加载参数并生成数据
Params = init_params();
if exist('rng','builtin'), rng(Params.rng_seed);
else, rand('seed', Params.rng_seed); randn('seed', Params.rng_seed); end

[Ident_data, Params] = gen_prbs_excitation(Params);
Model = build_dypls(Ident_data, Params);

% 2. 估算不确定性并计算 Tube
UncSets = estimate_uncertainty(Model, Ident_data, Params);
Tube = compute_tube(Ident_data, Model, UncSets, Params);

% 3. 计算 Assumption 3 的精确最大二次型
fprintf('\n=== Computing Maximum Quadratic Value ===\n');
P_e = Tube.P_e;
A_lv = Params.A_lv;
n_state_a = 2*Params.n_arx - 1;
beta_v = 0.15;

% 提取非零元素的索引（每个通道块的第1个位置）
idx_pert = zeros(A_lv, 1);
for a = 1:A_lv
    idx_pert(a) = (a-1)*n_state_a + 1;
end

% 提取 P_e^-1 的对应子矩阵 (A_lv x A_lv)
P_inv = P_e \ eye(size(P_e));
M_sub = P_inv(idx_pert, idx_pert);

% 最大特征值即为最大二次型除以 beta_v^2
lambda_max = max(eig(M_sub));
max_quad_val = (beta_v^2) * lambda_max;

fprintf('\n======================================================\n');
fprintf('  Params.init_state_variation (beta_v) = %.2f\n', beta_v);
fprintf('  Max eigenvalue of restricted P_e^-1 = %.4f\n', lambda_max);
fprintf('  Maximum e(0)^T P_e^{-1} e(0) = %.6f\n', max_quad_val);
fprintf('  Is it strictly < 1?  %s\n', mat2str(max_quad_val < 1));
fprintf('======================================================\n');
