function Data = run_ilmpc(Identification_data, Model, Tube, UncSets, Params)
% =========================================================================
% run_ilmpc -- Algorithm 2 of the paper (N-step lifting).
%
% COORDINATE CONVENTION (IMPORTANT):
%   Every signal stored in Data is a DEVIATION quantity:
%       Y_hist  = y_dev = y_abs - y_ss,  indexed by k = 1..N
%       V_hist  = v_dev = v_abs - u_ss,  indexed by k = 1..N (v(0..N-1))
%       ref     = Ref_points (= y_ref_phys - y_ss)
%   plot_results.m adds y_ss / u_ss for display.
%
%   The physical plant is the DEVIATION system  G: v_dev -> y_dev.
%   Controller acts on v_dev only. X_actual is the latent state of
%   this deviation system (x_actual = 0 <=> y_dev ≈ Y_mean_dev ≈ 0).
%
% DECISION VARIABLE: δΔt(0..N-1), length A_lv*N.
% =========================================================================

if exist('OCTAVE_VERSION', 'builtin'), pkg load control optim; end

N = Params.Ble; A_lv = Params.A_lv; n_in = Params.n_in; n_out = Params.n_out;
n_arx = Params.n_arx; N_batches = Params.N_batches;

barA = Model.barA; barB = Model.barB; barC = Model.barC;
nx = size(barA, 1); n_state_a = 2*n_arx - 1;
Q_pls = Model.Q; P_pls = Model.P; R_u = Model.R_u;

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

K_bar     = Tube.K;
h_v_tight = Tube.h_v_tight;
h_y_tight = Tube.h_y_tight;
H_v = Params.H_v; H_y = Params.H_y;
Q_w = Params.Q_weight; R_w = Params.R_weight;

mask = zeros(n_out, N); ref_traj = zeros(n_out, N);
for i = 1:numel(Params.Kp)
    k_idx = Params.Kp(i);
    mask(1, k_idx) = 1;
    ref_traj(1, k_idx) = Params.Ref_points(i);
end

System_d   = ss(Params.sys_d);
System_d_d = ss(Params.sys_d_d);
Noise_var  = Params.noise_var;

[Phi_x, x0_contrib, Phi_y, y0_const, Phi_v, v0_const] = ...
    build_lifting_matrices(N, nx, n_out, n_in, A_lv, barA, barB, barC, ...
        S_y, Q_pls, S_v, P_pls, X_mean_dev, Y_mean_dev, zeros(nx,1));

n_t   = A_lv * N;
K_cum = kron(tril(ones(N, N)), eye(A_lv));
Phi_x_K = Phi_x * K_cum; Phi_y_K = Phi_y * K_cum; Phi_v_K = Phi_v * K_cum;

% ---- Reference/mask blocks: block i <=> y(i),  i = 1..N ----
Mask_vec  = zeros(n_out * N, 1);
Y_ref_vec = zeros(n_out * N, 1);
for i = 1:N
    row_idx = (i-1)*n_out+1 : i*n_out;
    Y_ref_vec(row_idx) = ref_traj(:, i);
    Mask_vec(row_idx)  = mask(:, i)';
end
Mask_mat = diag(Mask_vec);
Q_blk    = kron(eye(N), Q_w);
R_blk    = kron(eye(N), R_w);

% ---- QP cost ----
% J(ddt) = ||Mask*(e_prev + Phi_y_K*ddt)||_Q^2 + ||ddt||_R^2,
% e_prev = Y_prev - Y_ref_vec.
% Octave form 0.5 x'Hx + f'x:
%   H = 2*(A_K + R_blk),  f = 2*G_qp_base*e_prev,
% where G_qp_base = Phi_y_K' * Mask * Q * Mask,  A_K = G_qp_base * Phi_y_K.
% (Minimizer identical to H=A_K+R_blk, f=G_qp_base*e_prev; J value halved.)
G_qp_base = Phi_y_K' * Mask_mat * Q_blk * Mask_mat;
A_K       = G_qp_base * Phi_y_K;
A_K       = 0.5 * (A_K + A_K');
H_qp_base = 2 * (A_K + R_blk);
H_qp_base = 0.5 * (H_qp_base + H_qp_base') + 1e-10 * eye(n_t);

H_v_blk  = kron(eye(N), H_v); H_y_blk = kron(eye(N), H_y);
A_v_full = H_v_blk * Phi_v_K; A_y_full = H_y_blk * Phi_y_K;

dt_max  = Params.delta_t_max; ddt_max = Params.ddt_max;
A_dt_K  = [K_cum; -K_cum];       b_dt_K  = dt_max  * ones(2*n_t, 1);
A_ddt   = [eye(n_t); -eye(n_t)]; b_ddt   = ddt_max * ones(2*n_t, 1);

A_ineq_base = [A_v_full; A_y_full; A_dt_K; A_ddt];
m_ineq      = size(A_ineq_base, 1);

x_nominal_init = zeros(nx, 1);
[t_feas, ~, ~, ~, ~] = Initial_Feasible(x_nominal_init, ...
                                        Identification_data, Model, Tube, Params);
t_vec_prev = t_feas(:);
V_prev = Phi_v * t_vec_prev + v0_const;   % v(0..N-1), n_in*N
Y_prev = Phi_y * t_vec_prev + y0_const;   % y(1..N),   n_out*N

Y_hist     = zeros(n_out, N, N_batches);
V_hist     = zeros(n_in,  N, N_batches);
Y_nom_hist = zeros(n_out, N, N_batches+1);
V_nom_hist = zeros(n_in,  N, N_batches+1);
Err_nom    = zeros(N_batches, 1);
Err_act    = zeros(N_batches, 1);
QP_times   = zeros(N_batches, 1);
J_batch    = zeros(N_batches, 1);

eta_in_flags = -ones(N_batches, N, A_lv);
r_in_flags   = -ones(N_batches, N, n_out);
r_evaluated  = (norm(UncSets.W_r, 'fro') > 0) && any(UncSets.r_bound > 0);

% Initial (pre-batch-1) nominal trajectories: y(1..N), v(0..N-1)
Y_nom_hist(:, :, 1) = reshape(Y_prev, n_out, N);
V_nom_hist(:, :, 1) = reshape(V_prev, n_in,  N);

for b = 1:N_batches
    % --- Initial latent state for this batch ---
    if isfield(Params, 'init_state_vector') && ~isempty(Params.init_state_vector)
        % 直接使用 main.m 中计算好的精确扰动向量
        x_init_batch = x_nominal_init + Params.init_state_vector(:);
    elseif isfield(Params, 'init_state_variation') && Params.init_state_variation > 0
        % 保留原来的随机标量方式作为备选
        x_init_batch = x_nominal_init;
        dir_vec = randn(A_lv, 1); dir_vec = dir_vec / max(norm(dir_vec), eps);
        for a = 1:A_lv
            idx_u = (a-1)*n_state_a + 1;
            x_init_batch(idx_u) = x_init_batch(idx_u) + Params.init_state_variation * dir_vec(a);
        end
    else
        x_init_batch = x_nominal_init;
    end

    % --- QP: solve for δΔt(0..N-1) ---
    b_ineq = [repmat(h_v_tight, N, 1) - H_v_blk * V_prev;
              repmat(h_y_tight, N, 1) - H_y_blk * Y_prev;
              b_dt_K; b_ddt];
    f_qp   = 2 * G_qp_base * (Y_prev - Y_ref_vec);

    tic;
    [ddt_opt, ~, info, ~] = qp(zeros(n_t,1), H_qp_base, f_qp, [], [], ...
        -inf(n_t,1), inf(n_t,1), -inf(m_ineq,1), A_ineq_base, b_ineq);
    QP_times(b) = toc;

    if isstruct(info) && isfield(info, 'info'), info_code = info.info; else, info_code = info; end
    if info_code ~= 0
        warning('run_ilmpc: QP failed (batch %d, info=%d); using previous trajectory.', b, info_code);
        ddt_opt = zeros(n_t, 1);
    end

    % --- Reconstruct nominal trajectory ---
    t_vec_new = t_vec_prev + K_cum * ddt_opt;
    V_new     = Phi_v * t_vec_new + v0_const;
    Y_new     = Phi_y * t_vec_new + y0_const;
    X_new     = Phi_x * t_vec_new + x0_contrib;

    T_nom = reshape(t_vec_new, A_lv, N);      % t(0..N-1)
    V_nom = reshape(V_new, n_in, N);          % v(0..N-1)
    Y_nom = reshape(Y_new, n_out, N);         % y(1..N)

    % --- Record J^{b,*} (matches Theorem 2: V^b + sum||δΔt||_R^2) ---
    Y_next_vec = Y_prev + Phi_y_K * ddt_opt;
    e_next_vec = Mask_mat * (Y_next_vec - Y_ref_vec);
    J_batch(b) = e_next_vec' * Q_blk * e_next_vec + ddt_opt' * R_blk * ddt_opt;

    % --- Batch execution with tube feedback ---
    Y_batch = zeros(n_out, N);   % y(1..N)
    V_batch = zeros(n_in,  N);   % v(0..N-1)
    U_batch = zeros(A_lv,  N);   % u(1..N)
    T_batch = zeros(A_lv,  N);   % t(0..N-1)

    X_nominal = x_nominal_init;
    X_actual  = x_init_batch;

    % Initialize physical deviation plant at y(0) = Y_mean_dev + S_y Q barC x_init
    U_act_0 = barC * X_actual;
    Y_0     = Y_mean_dev + S_y * Q_pls * U_act_0;
    DC_gain = System_d.c * ((eye(size(System_d.a)) - System_d.a) \ System_d.b);
    v_init_actual = DC_gain \ Y_0;
    Sys_state     = (eye(size(System_d.a)) - System_d.a) \ (System_d.b * v_init_actual);
    Sys_state_dis = zeros(size(System_d_d.a, 1), 1);

    for k = 1:N
        % --- Tube feedback (deviation form): apply v(k-1) ---
        e_k      = X_actual - X_nominal;
        t_actual = T_nom(:, k) + K_bar * e_k;                 % t(k-1) + correction
        v_actual = X_mean_dev + S_v * P_pls * t_actual;       % v(k-1)
        V_batch(:, k) = v_actual;
        T_batch(:, k) = t_actual;

        % Nominal state advance: x(k) -> x(k+1)
        X_nominal = barA * X_nominal + barB * T_nom(:, k);

        % --- Physical plant simulation (deviation system) ---
        d_raw  = randn(1,1);
        d_raw  = max(min(d_raw, 3.0), -3.0);
        d_dist = sqrt(Noise_var) * d_raw;

        Sys_state     = System_d.a   * Sys_state     + System_d.b   * v_actual;
        y1            = System_d.c   * Sys_state;
        Sys_state_dis = System_d_d.a * Sys_state_dis + System_d_d.b * d_dist;
        y2            = System_d_d.c * Sys_state_dis;
        y_actual_k    = y1 + y2;                               % y(k)
        Y_batch(:, k) = y_actual_k;

        % --- Output-score estimate via PLS left inverse R_u ---
        u_actual_k = R_u' * (S_y \ (y_actual_k - Y_mean_dev)); % u(k)
        U_batch(:, k) = u_actual_k;

% --- Coverage flags: eta(k) uses ALL channels' u(k-1..k-n_arx), t(k-1..k-n_arx) ---
if k >= n_arx + 1
    for a = 1:A_lv
        theta = Model.Theta{a}(:);
        u_pred = 0;
        for bb = 1:A_lv                       % ← 关键：用 bb，不覆盖 batch b
            idx_u = (bb-1)*n_arx + (1:n_arx);
            idx_t = A_lv*n_arx + idx_u;
            u_lag = U_batch(bb, k-1 : -1 : k-n_arx);
            t_lag = T_batch(bb, k   : -1 : k-n_arx+1);
            u_pred = u_pred + theta(idx_u)'*u_lag' + theta(idx_t)'*t_lag';
        end
        eta_a = u_actual_k(a) - u_pred;
        eta_in_flags(b, k, a) = double(abs(eta_a) <= UncSets.eta_bound(a));  % ← b 是 batch 号
    end
end

        if r_evaluated
            Y0_meas = S_y \ (y_actual_k - Y_mean_dev);
            r_k     = Y0_meas - Q_pls * u_actual_k;
            for i = 1:n_out
                r_in_flags(b, k, i) = double(abs(r_k(i)) <= UncSets.r_bound(i));
            end
        end

        % --- Latent state shift x_actual(k) -> x_actual(k+1) ---
        X_tt = zeros(nx, 1);
        for a = 1:A_lv
            row     = (a-1)*n_state_a+1 : a*n_state_a;
            X_temp1 = X_actual(row);
            u_seq   = X_temp1(1:n_arx-1);
            t_seq   = X_temp1(n_arx+1:2*n_arx-2);
            X_tt(row) = [u_actual_k(a); u_seq; t_actual(a); t_seq];
        end
        X_actual = X_tt;
    end

    Y_hist(:, :, b) = Y_batch;
    V_hist(:, :, b) = V_batch;


    % --- Tracking errors at prescribed instants (Weighted by Q_trk) ---
    ea = 0; en = 0;
    q_trk_1 = Params.Q_weight(1, 1); % 提取 y1 的追踪权重
    for i = 1:numel(Params.Kp)
        k_i = Params.Kp(i);
        ea = ea + q_trk_1 * (Y_batch(1, k_i) - Params.Ref_points(i))^2;
        en = en + q_trk_1 * (Y_nom(1,   k_i) - Params.Ref_points(i))^2;
    end
    Err_act(b) = sqrt(ea);
    Err_nom(b) = sqrt(en);

    % --- Prepare next batch ---
    t_vec_prev = t_vec_new;
    V_prev     = V_new;
    Y_prev     = Y_new;
    Y_nom_hist(:, :, b+1) = Y_nom;
    V_nom_hist(:, :, b+1) = V_nom;

    fprintf('Batch %d: nom err=%.4e, act err=%.4e, J=%.4e, QP time=%.4fs\n', ...
        b, Err_nom(b), Err_act(b), J_batch(b), QP_times(b));
end

% ---- Package (all deviation-form) ----
Data.Y_hist             = Y_hist;              % y(1..N)
Data.V_hist             = V_hist;              % v(0..N-1)
Data.Y_nominal_hist     = Y_nom_hist;          % y(1..N), incl. pre-batch-1
Data.V_nominal_hist     = V_nom_hist;          % v(0..N-1)
Data.Error_norm_actual  = Err_act;
Data.Error_norm_nominal = Err_nom;
Data.QP_times           = QP_times;
Data.J_batch            = J_batch;
Data.eta_in_flags       = eta_in_flags;
Data.r_in_flags         = r_in_flags;

Data.u_ss = Identification_data.u_ss;
Data.y_ss = Identification_data.y_ss;

if isfield(Model, 'X_mean_dev')
    Data.X_mean_dev = Model.X_mean_dev;
    Data.Y_mean_dev = Model.Y_mean_dev;
else
    Data.X_mean_dev = Identification_data.X_mean;
    Data.Y_mean_dev = Identification_data.Y_mean;
end
end
