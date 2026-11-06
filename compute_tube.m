function Tube = compute_tube(Identification_data, Model, UncSets, Params)
% =========================================================================
% compute_tube
%
% 使用离散 Lyapunov 方程计算 RPI 椭球。
% 包含自适应 alpha 搜索，建立 alpha 与紧约束之间的关联。
% =========================================================================

if exist('OCTAVE_VERSION', 'builtin')
    pkg load control;
end

%% 1. 模型矩阵与维度
A_bar = Model.barA;
B_bar = Model.barB;
C_bar = Model.barC;
E_bar = Model.barE;
P_pls = Model.P;
Q_pls = Model.Q;

[nx, ~] = size(A_bar);
[~, nt] = size(B_bar);

A_lv = Params.A_lv;
n_arx = Params.n_arx;
n_state = 2*n_arx - 1;

%% 2. 不确定性边界提取
eta_bound = UncSets.eta_bound(:);
W_r = UncSets.W_r;
W_x = UncSets.W_x;

if isfield(Model, 'X_std')
    S_v = diag(Model.X_std(:));
    S_y = diag(Model.Y_std(:));
else
    S_v = diag(Identification_data.X_std(:));
    S_y = diag(Identification_data.Y_std(:));
end

H_v = Params.H_v;
H_y = Params.H_y;
h_v_phys = Params.h_v_phys(:);
h_y_phys = Params.h_y_phys(:);
n_v_constr = size(H_v,1);
n_y_constr = size(H_y,1);

%% 3. 结构化 LQR 反馈增益
use_lqr = true;
if isfield(Params, 'tube_use_lqr') && ~isempty(Params.tube_use_lqr)
    use_lqr = logical(Params.tube_use_lqr);
end

if use_lqr
    if isfield(Params, 'lqr_Q_scale'), q_scale = Params.lqr_Q_scale; else, q_scale = 1.0; end
    if isfield(Params, 'lqr_R_scale'), r_scale = Params.lqr_R_scale; else, r_scale = 1.0; end

    % 构建结构化的 Q 矩阵，只惩罚真实的物理状态
    Q_lqr = zeros(n_state, n_state);
    Q_lqr(1, 1) = 1.00;
    if n_arx >= 2
        Q_lqr(2:n_arx, 2:n_arx) = 0.30 * eye(n_arx - 1);
    end
    if n_state >= n_arx + 1
        Q_lqr(n_arx + 1, n_arx + 1) = 0.05;
    end
    Q_lqr = q_scale * Q_lqr;

    K_blk = zeros(nt, nx);
    for a = 1:A_lv
        idx = (a-1)*n_state+1 : a*n_state;
        try
            K_lqr_i = dlqr(A_bar(idx, idx), B_bar(idx, a), Q_lqr, r_scale);
        catch ME
            error('compute_tube: dlqr failed at LV %d: %s', a, ME.message);
        end
        K_blk(a, idx) = -K_lqr_i;
    end
    K_bar = K_blk;
    A_K   = A_bar + B_bar * K_bar;
else
    K_bar = zeros(nt, nx);
    A_K   = A_bar;
end

rho_K = max(abs(eig(A_K)));
if rho_K >= 1
    error('compute_tube: rho(A_K) = %.6f >= 1; cannot build tube.', rho_K);
end

%% 4. 聚合扰动集 W_e
beta = Params.tube_beta;
W_eta = A_lv * diag(eta_bound.^2);
W_e = (1 + beta) * (B_bar * K_bar) * W_x * (B_bar * K_bar)' ...
    + (1 + 1/beta) * E_bar * W_eta * E_bar';
W_e = 0.5 * (W_e + W_e');

%% 5. 自适应 Alpha 搜索与 Lyapunov 方程求解 P_e
rho2 = rho_K^2;
alpha_limit = 1/rho2 - 1;

% 读取候选网格，并附加一些边界自适应点
if isfield(Params, 'tube_alpha_grid') && ~isempty(Params.tube_alpha_grid)
    base_candidates = Params.tube_alpha_grid(:)';
else
    base_candidates = [0.01, 0.03, 0.1, 0.3, 1.0];
end
adaptive = alpha_limit * [0.02, 0.1, 0.3, 0.6];
candidates = unique([base_candidates, adaptive]);
% 过滤：必须大于 0 且满足稳定性条件 sqrt(1+alpha)*rho_K < 1
candidates = candidates(isfinite(candidates) & candidates > 0 & (1+candidates)*rho2 < 1-1e-10);

if isempty(candidates)
    error('compute_tube: no valid alpha candidate satisfying stability.');
end

best_P = [];
best_alpha = NaN;
best_score = Inf;

fprintf('\nLyapunov RPI ellipsoid construction (alpha search):\n');
fprintf('  %-10s | %-14s | %-12s | %-10s\n', 'alpha', 'trace(P_e)', 'rcond(P_e)', 'status');

for j = 1:numel(candidates)
    alpha = candidates(j);
    As = sqrt(1+alpha) * A_K;
    % 注意 Qshape 依赖于 alpha
    Qshape = (1+1/alpha) * W_e + 1e-12 * eye(nx);

    try
        P_test = dlyap(As, Qshape);
        P_test = 0.5 * (P_test + P_test');

        [L, flag] = chol(P_test, 'lower');
        if flag ~= 0
            fprintf('  %-10.4f | (not pos-def)\n', alpha);
            continue;
        end

        % 校验不变性 (对 alpha 和 beta 组合进行验证)
        image_shape = (1+alpha)*A_K*P_test*A_K' + (1+1/alpha)*W_e;
        inv_ratio = max(eig(symmetric_local(L \ image_shape / L')));

        if inv_ratio > 1 + 1e-7
            fprintf('  %-10.4f | (invariant fail)\n', alpha);
            continue;
        end

        score = trace(P_test);
        fprintf('  %-10.4f | %-14.4e | %-12.4e | OK\n', alpha, score, rcond(P_test));

        if score < best_score
            best_score = score;
            best_alpha = alpha;
            best_P = P_test;
        end
    catch
        % ignore failed alpha
    end
end

if isempty(best_P)
    error('compute_tube: Lyapunov solver failed for all alpha candidates.');
end

P_e = best_P;
fprintf('  --> Selected alpha: %.4f, trace(P_e): %.4e\n', best_alpha, best_score);

%% 6. 约束收紧计算
M_v = S_v * P_pls * K_bar;
delta_v = zeros(n_v_constr, 1);
if norm(M_v(:)) > 0
    for i = 1:n_v_constr
        h_i = H_v(i,:);
        t1 = h_i * M_v * P_e * M_v' * h_i';
        t2 = h_i * M_v * W_x * M_v' * h_i';
        if t1 < 0, t1 = 0; end
        if t2 < 0, t2 = 0; end
        delta_v(i) = sqrt(t1) + sqrt(t2);
    end
end

C_y = S_y * Q_pls * C_bar;
delta_y = zeros(n_y_constr, 1);
for i = 1:n_y_constr
    h_i = H_y(i,:);
    t1 = h_i * C_y * P_e * C_y' * h_i';
    t2 = h_i * S_y * W_r * S_y' * h_i';
    if t1 < 0, t1 = 0; end
    if t2 < 0, t2 = 0; end
    delta_y(i) = sqrt(t1) + sqrt(t2);
end

%% 7. 打包结果
h_v_tight = h_v_phys - delta_v;
h_y_tight = h_y_phys - delta_y;

infeas_v = find(h_v_tight < 0);
infeas_y = find(h_y_tight < 0);
if isempty(infeas_v) && isempty(infeas_y)
    fprintf('  *** TUBE FEASIBLE ***\n');
else
    fprintf('\n  *** INFEASIBLE: v=%s, y=%s ***\n', mat2str(infeas_v'), mat2str(infeas_y'));
end

Tube.P_e = P_e;
Tube.K = K_bar;
Tube.A_K = A_K;
Tube.spectral_radius = rho_K;
Tube.alpha = best_alpha;            % 记录当前选定的最优 alpha
Tube.beta = beta;
Tube.W_eta = W_eta;
Tube.W_x = W_x;
Tube.W_r = W_r;
Tube.W_e = W_e;
Tube.delta_v = delta_v;
Tube.delta_y = delta_y;
Tube.h_v_tight = h_v_tight;
Tube.h_y_tight = h_y_tight;
Tube.h_v_phys = h_v_phys;
Tube.h_y_phys = h_y_phys;
Tube.H_v = H_v;
Tube.H_y = H_y;
Tube.tube_use_lqr = use_lqr;
if use_lqr
    Tube.lqr_Q_scale = q_scale;
    Tube.lqr_R_scale = r_scale;
    Tube.Q_lqr = Q_lqr;
end

fprintf('\ncompute_tube: done.\n');
end

% =========================================================================
% 辅助函数
% =========================================================================
function S = symmetric_local(M)
    S = real((M+M')/2);
end
