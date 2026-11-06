function Result = verify_unreachability(Identification_data, Model, Tube, Params)
% =========================================================================
% verify_unreachability -- Direct solution of Problem (55), N-STEP version.
%
%   Solves  min_{t_vec} ||Mask*(y_dev(1..N) - y_ref(1..N))||_Q^2
%   s.t.    PHYSICAL tightened bounds on v(0..N-1), y(1..N).
%
% All quantities are DEVIATION-form. Kp reference values are provided
% via Params.Ref_points (= y_ref_phys - y_ss).
% =========================================================================

N = Params.Ble; A_lv = Params.A_lv; n_out = Params.n_out; n_in = Params.n_in;
barA = Model.barA; barB = Model.barB; barC = Model.barC;
nx = size(barA, 1);

Q_pls = Model.Q; P_pls = Model.P;

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

[~, ~, Phi_y, y0_const, Phi_v, v0_const] = ...
    build_lifting_matrices(N, nx, n_out, n_in, A_lv, barA, barB, barC, ...
        S_y, Q_pls, S_v, P_pls, X_mean_dev, Y_mean_dev, zeros(nx,1));

% mask / ref blocks: block i <=> y(i),  i = 1..N
mask = zeros(n_out, N); ref_traj = zeros(n_out, N);
for i = 1:numel(Params.Kp)
    k_idx = Params.Kp(i);
    mask(1, k_idx) = 1;
    ref_traj(1, k_idx) = Params.Ref_points(i);
end
Mask_vec  = zeros(n_out * N, 1);
Y_ref_vec = zeros(n_out * N, 1);
for i = 1:N
    row_idx = (i-1)*n_out+1 : i*n_out;
    Y_ref_vec(row_idx) = ref_traj(:, i);
    Mask_vec(row_idx)  = mask(:, i)';
end
Mask_mat = diag(Mask_vec);
Q_blk    = kron(eye(N), Params.Q_weight);

% ---- QP: min 0.5 x'Hx + f'x, H = 2*G_qp*Phi_y, f = 2*G_qp*(y0-Y_ref) ----
G_qp = Phi_y' * Mask_mat * Q_blk * Mask_mat;   % n_t × n_out*N
A_K  = G_qp * Phi_y;
A_K  = 0.5 * (A_K + A_K');
n_t  = A_lv * N;
H_qp = 2 * A_K + 1e-10 * eye(n_t);
H_qp = 0.5 * (H_qp + H_qp');
f_qp = 2 * G_qp * (y0_const - Y_ref_vec);

H_v_blk = kron(eye(N), Params.H_v);
H_y_blk = kron(eye(N), Params.H_y);
A_ineq  = [H_v_blk * Phi_v;
           H_y_blk * Phi_y];
b_ineq  = [repmat(Tube.h_v_tight, N, 1) - H_v_blk * v0_const;
           repmat(Tube.h_y_tight, N, 1) - H_y_blk * y0_const];

m = size(A_ineq, 1);

[t_vec_opt, ~, info, ~] = qp(zeros(n_t,1), H_qp, f_qp, [], [], ...
    -inf(n_t,1), inf(n_t,1), -inf(m,1), A_ineq, b_ineq);

if isstruct(info) && isfield(info, 'info'), info_code = info.info; else, info_code = info; end
if info_code ~= 0
    warning('verify_unreachability: QP returned info = %d.', info_code);
end

Y_stack = Phi_y * t_vec_opt + y0_const;   % y(1..N), n_out*N
Y_nom   = reshape(Y_stack, n_out, N);
V_nom   = reshape(Phi_v * t_vec_opt + v0_const, n_in, N);

Result.t_vec_opt = t_vec_opt;
Result.Y_nom     = Y_nom;         % y(1..N)
Result.V_nom     = V_nom;         % v(0..N-1)
Result.J_opt     = 0.5*t_vec_opt'*H_qp*t_vec_opt + f_qp'*t_vec_opt;

% y_s^* at prescribed instants (block k_idx is y(k_idx))
Result.y_s_at_Kp       = zeros(numel(Params.Kp), n_out);
Result.nom_error_at_Kp = zeros(numel(Params.Kp), 1);
for i = 1:numel(Params.Kp)
    k_idx = Params.Kp(i);
    Result.y_s_at_Kp(i, :)    = Y_nom(:, k_idx)';
    Result.nom_error_at_Kp(i) = Y_nom(1, k_idx) - Params.Ref_points(i);
end

tol = 1e-6; viol = A_ineq * t_vec_opt - b_ineq;
Result.active   = find(viol > -tol);
Result.n_active = numel(Result.active);
end
