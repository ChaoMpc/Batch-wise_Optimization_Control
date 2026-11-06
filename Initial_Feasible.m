function [t_feasible, u_feasible, y_feasible, v_feasible, x_feasible] = ...
        Initial_Feasible(x_nominal_init, Identification_data, Model, Tube, Params)
% =========================================================================
% Initial_Feasible -- Find a feasible nominal latent-input trajectory.
%
% N-STEP version. All signals DEVIATION quantities.
%
% IMPORTANT:
%   * Only PHYSICAL input/output constraints are imposed here.
%     The batch-increment bound |Δt| <= Δt_max is a BATCH UPDATE
%     constraint, not an absolute bound on t. It is therefore NOT
%     applied in Initial_Feasible (there is no previous batch).
%   * Output/input lifting uses N steps:
%         y(1..N)  (N  blocks),  v(0..N-1) (N blocks)
% =========================================================================

    N     = Params.Ble;
    A_lv  = Params.A_lv;
    n_in  = Params.n_in;
    n_out = Params.n_out;
    nx    = length(x_nominal_init);

    barA = Model.barA; barB = Model.barB; barC = Model.barC;
    P_pls = Model.P; Q = Model.Q;

    if isfield(Model, 'X_mean_dev')
        X_mean_dev = Model.X_mean_dev(:);
        Y_mean_dev = Model.Y_mean_dev(:);
    else
        X_mean_dev = Identification_data.X_mean(:);
        Y_mean_dev = Identification_data.Y_mean(:);
    end
    if isfield(Model, 'X_std')
        S_v = diag(Model.X_std(:));
        S_y = diag(Model.Y_std(:));
    else
        S_v = diag(Identification_data.X_std(:));
        S_y = diag(Identification_data.Y_std(:));
    end

    h_v_tight = Tube.h_v_tight; h_y_tight = Tube.h_y_tight;
    H_v = Params.H_v; H_y = Params.H_y;

R_weight = Params.R_weight;
if isscalar(R_weight), R_weight = R_weight * eye(A_lv); end
% 让 batch 1 的初始解更保守（靠近稳态），便于展示 ILC 收敛过程
if isfield(Params, 'init_R_weight_scale') && ~isempty(Params.init_R_weight_scale)
    R_weight = Params.init_R_weight_scale * R_weight;
end

    % A_pow{m} = barA^{m-1}
    A_powers = cell(N + 1, 1); A_powers{1} = eye(nx);
    for i = 2:N+1, A_powers{i} = barA * A_powers{i - 1}; end

    n_t   = A_lv * N;
    Phi_x = zeros(nx * N, n_t);
    for k = 0:N-1
        r = k*nx+1 : (k+1)*nx;
        for j = 0:k
            c = j*A_lv+1 : (j+1)*A_lv;
            Phi_x(r, c) = A_powers{k - j + 1} * barB;
        end
    end

    Gamma_x = zeros(nx * N, nx);
    for k = 0:N-1
        r = k*nx+1 : (k+1)*nx;
        Gamma_x(r, :) = A_powers{k + 2};
    end
    x0_contrib = Gamma_x * x_nominal_init;

    C_y      = S_y * Q * barC;
    Phi_y    = zeros(n_out * N, n_t);
    y0_const = zeros(n_out * N, 1);
    for k = 1:N
        r  = (k-1)*n_out+1 : k*n_out;
        sr = (k-1)*nx+1    : k*nx;
        Phi_y(r, :)  = C_y * Phi_x(sr, :);
        y0_const(r)  = Y_mean_dev + C_y * x0_contrib(sr);
    end

    M_v      = S_v * P_pls;
    Phi_v    = zeros(n_in * N, n_t);
    v0_const = zeros(n_in * N, 1);
    for k = 0:N-1
        r = k*n_in+1 : (k+1)*n_in;
        c = k*A_lv+1 : (k+1)*A_lv;
        Phi_v(r, c) = M_v;
        v0_const(r) = X_mean_dev;
    end

    n_constr_v = size(H_v, 1); n_constr_y = size(H_y, 1);
    n_ineq     = N * (n_constr_v + n_constr_y);

    Aineq = zeros(n_ineq, n_t); bineq = zeros(n_ineq, 1);
    for k = 1:N
        row_base = (k-1) * (n_constr_v + n_constr_y);
        % v(k-1) 对应的物理输入约束
        rows_v = row_base + 1 : row_base + n_constr_v;
        Aineq(rows_v, :) = H_v * Phi_v((k-1)*n_in+1 : k*n_in, :);
        bineq(rows_v)    = h_v_tight - H_v * v0_const((k-1)*n_in+1 : k*n_in);
        % y(k) 对应的物理输出约束
        rows_y = row_base + n_constr_v + 1 : row_base + n_constr_v + n_constr_y;
        Aineq(rows_y, :) = H_y * Phi_y((k-1)*n_out+1 : k*n_out, :);
        bineq(rows_y)    = h_y_tight - H_y * y0_const((k-1)*n_out+1 : k*n_out);
    end

    % NOTE: |Δt| <= Δt_max and |δΔt| <= δΔt_max are batch-update bounds
    % enforced in run_ilmpc, NOT here (no previous batch in Initial_Feasible).

    H_obj = kron(eye(N), R_weight);
    f_obj = zeros(n_t, 1);
    H_obj = 0.5 * (H_obj + H_obj');

    x0   = zeros(n_t, 1);
    lb   = -inf(n_t, 1);
    ub   =  inf(n_t, 1);
    m    = size(Aineq, 1);

    [x_opt, ~, info, ~] = qp(x0, H_obj, f_obj, [], [], lb, ub, ...
                             -inf(m,1), Aineq, bineq);

    if isstruct(info) && isfield(info, 'info'), info_code = info.info; else, info_code = info; end
    if info_code ~= 0
        warning('Initial_Feasible: QP returned info = %d.', info_code);
    end
    if isempty(x_opt) || numel(x_opt) ~= n_t
        error('Initial_Feasible: QP did not return a valid solution.');
    end

    t_vec      = x_opt;
    t_feasible = reshape(t_vec, A_lv, N);                 % t(0..N-1)
    x_vec      = Phi_x * t_vec + x0_contrib;              % x(1..N)
    x_feasible = [x_nominal_init, reshape(x_vec, nx, N)]; % x(0..N)

    u_feasible = zeros(A_lv, N);
    for j = 1:N, u_feasible(:, j) = barC * x_feasible(:, j+1); end   % u(1..N)

    y_vec      = Phi_y * t_vec + y0_const;                % y(1..N)
    y_feasible = reshape(y_vec, n_out, N);

    v_vec      = Phi_v * t_vec + v0_const;                % v(0..N-1)
    v_feasible = reshape(v_vec, n_in, N);
end
