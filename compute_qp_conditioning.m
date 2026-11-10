function R = compute_qp_conditioning(Problem, Model, Params)
% 同时给出潜空间与"等价原始物理空间"的 QP Hessian 条件数。

    H = Problem.H;

    % ---- 潜空间 ----
    L_lat     = Problem.L;
    H_lat     = 2*(L_lat'*L_lat + Problem.Pprox);
    H_lat     = 0.5*(H_lat + H_lat');
    R.lat_dim = size(H_lat,1);
    R.lat_cond = cond(H_lat);
    eig_lat   = eig(H_lat);
    R.lat_eigmin = min(eig_lat);
    R.lat_eigmax = max(eig_lat);

    % ---- 原始物理空间 ----
    % 物理决策变量 v_vec = [v(k_1); v(k_2); ...]，每个 k 是 n_opt 维。
    % 它到 latent 决策 t_vec 的映射：
    %   t_vec = kron(I_H, R_x'*S_v^{-1}) * v_vec + const
    Rx    = Model.R_x;
    Sv_inv = diag(1 ./ Model.X_std);
    M_map = kron(eye(H), Rx' * Sv_inv);     % (n_lv*H) x (n_opt*H)

    L_raw = L_lat * M_map;                   % tracking-error 对 raw v 的线性映射
    P_raw = M_map' * Problem.Pprox * M_map;  % proximal 在 raw 空间的拉回
    H_raw = 2*(L_raw'*L_raw + P_raw);
    H_raw = 0.5*(H_raw + H_raw');
    R.raw_dim  = size(H_raw,1);
    R.raw_cond = cond(H_raw);
    eig_raw    = eig(H_raw);

    % ================= 修复数值精度造成的微小负值 =================
    tol_eig = R.raw_dim * eps(norm(H_raw, 2)); % 机器精度的动态阈值
    eig_raw(eig_raw < 0 & eig_raw > -tol_eig) = 0; % 将极小负值拉平为 0

    R.raw_eigmin = min(eig_raw);
    if R.raw_eigmin == 0
        R.raw_eigmin_str = '0 (numerical truncation)'; % 用于打印显示
    else
        R.raw_eigmin_str = sprintf('%.3e', R.raw_eigmin);
    end
    R.raw_eigmax = max(eig_raw);

    % ---- 决策变量缩减 ----
    R.n_opt  = Params.n_in_opt;
    R.n_lv   = Model.A_lv;
    R.reduce = 100 * (1 - R.n_lv/R.n_opt);

    fprintf('\n[QP conditioning]\n');
    fprintf('  dim: raw %d vs latent %d (%.1f%% reduction)\n', ...
        R.raw_dim, R.lat_dim, R.reduce);
    fprintf('  cond(H): raw %.4e | latent %.4e\n', R.raw_cond, R.lat_cond);
    fprintf('  eig range (raw)   : [%.3e, %.3e]\n', R.raw_eigmin, R.raw_eigmax);
    fprintf('  eig range (latent): [%.3e, %.3e]\n', R.lat_eigmin, R.lat_eigmax);
end
