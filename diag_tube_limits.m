function diag_tube_limits(Ident_data, Model, Tube, Params)
% Read-only. Confirms whether the physical constraints are compatible
% with the tube size (no code modification, no Params change).

N = Params.Ble; A_lv = Params.A_lv; n_out = Params.n_out; n_in = Params.n_in;
barA = Model.barA; barB = Model.barB; barC = Model.barC;
nx = size(barA,1);

% --- Y_std information ---
if isfield(Model,'Y_std'), Y_std = Model.Y_std(:); else, Y_std = Model.Y_std; end
if isfield(Model,'Y_mean_dev'), Y_mean = Model.Y_mean_dev(:); else, Y_mean = Ident_data.Y_mean(:); end
fprintf('\n===== Tube Feasibility Diagnostic =====\n');
fprintf('Y_std   = [%s]\n', num2str(Y_std', '%.4e '));
fprintf('Y_mean_dev = [%s]\n', num2str(Y_mean', '%.4e '));

% --- Physical constraint bounds in "physical" units ---
h_y = Params.h_y_phys(:);  h_v = Params.h_v_phys(:);
fprintf('h_y_phys (dev bounds, physical units) = [%s]\n', num2str(h_y', '%.4e '));
fprintf('h_v_phys (dev bounds, physical units) = [%s]\n', num2str(h_v', '%.4e '));

% --- Tube tightening (dev form) ---
d_y = Tube.delta_y(:); d_v = Tube.delta_v(:);
fprintf('delta_y = [%s]\n', num2str(d_y', '%.4e '));
fprintf('delta_v = [%s]\n', num2str(d_v', '%.4e '));

% --- Ratio: how many σ is each constraint? ---
%   The constraint / the corresponding tightening.
fprintf('\n[Constraint vs. tightening ratio (must be > 1 for feasibility)]\n');
fprintf('  y: h_y / delta_y = [%s]\n', num2str((h_y./max(d_y,eps))', '%.3f '));
fprintf('  v: h_v / delta_v = [%s]\n', num2str((h_v./max(d_v,eps))', '%.3f '));

% --- Bound expressed in training-sigma ---
fprintf('\n[Physical constraint expressed in training-sigma (y only)]\n');
%   h_y first half corresponds to + bound, second half - bound; same magnitude
%   (since [max; -min]), so h_y(1:n_out) is the + deviation limit.
fprintf('  +y bound in sigma: [%s]\n', num2str((h_y(1:n_out) ./ Y_std)', '%.3f '));

% --- Nominal-only QP feasibility with h_phys (no tightening) ---
%   Build lifting matrices without any tube, use Params.h_v_phys / h_y_phys.
[Phi_x, x0c, Phi_y, y0c, Phi_v, v0c] = build_lifting_matrices( ...
    N, nx, n_out, n_in, A_lv, barA, barB, barC, ...
    diag(Model.Y_std), Model.Q, diag(Model.X_std), Model.P, ...
    Model.X_mean_dev, Model.Y_mean_dev, zeros(nx,1));
n_t = A_lv*N;
H_v_blk = kron(eye(N), Params.H_v);
H_y_blk = kron(eye(N), Params.H_y);

% --- (i) Nominal QP with PHYSICAL bounds (no tightening) ---
A_ineq_phys = [H_v_blk * Phi_v; H_y_blk * Phi_y];
b_ineq_phys = [repmat(h_v, N, 1) - H_v_blk * v0c;
               repmat(h_y, N, 1) - H_y_blk * y0c];
m_ineq = size(A_ineq_phys, 1);

% Use tiny quadratic cost for feasibility test
H_qp = 1e-3 * eye(n_t); f_qp = zeros(n_t,1);
[t_phys, ~, info_phys, ~] = qp(zeros(n_t,1), H_qp, f_qp, [], [], ...
    -inf(n_t,1), inf(n_t,1), -inf(m_ineq,1), A_ineq_phys, b_ineq_phys);
if isstruct(info_phys) && isfield(info_phys,'info'), ic = info_phys.info; else, ic = info_phys; end
fprintf('\n[Feasibility with PHYSICAL bounds, no tightening]: info=%d (%s)\n', ...
    ic, ternary_str(ic == 0, 'FEASIBLE', 'INFEASIBLE'));

% --- (ii) Nominal QP with TUBE-tightened bounds ---
b_ineq_tight = [repmat(Tube.h_v_tight, N, 1) - H_v_blk * v0c;
                repmat(Tube.h_y_tight, N, 1) - H_y_blk * y0c];
[t_tight, ~, info_tight, ~] = qp(zeros(n_t,1), H_qp, f_qp, [], [], ...
    -inf(n_t,1), inf(n_t,1), -inf(m_ineq,1), A_ineq_phys, b_ineq_tight);
if isstruct(info_tight) && isfield(info_tight,'info'), ic2 = info_tight.info; else, ic2 = info_tight; end
fprintf('[Feasibility with TIGHTENED bounds]: info=%d (%s)\n', ...
    ic2, ternary_str(ic2 == 0, 'FEASIBLE', 'INFEASIBLE'));

% --- (iii) What would be the "barely feasible" tightening? ---
%   Compute max allowable delta such that h - delta > 0 (per component)
delta_max_allowable_y = h_y;   % trivially h_y (delta = h_y is the limit)
delta_max_allowable_v = h_v;
fprintf('\n[Maximum allowable tightening for non-emptiness]:\n');
fprintf('  y: need delta_y(i) < h_y(i)  =>  delta_y(1:2) < [%s]\n', ...
    num2str(h_y(1:2)', '%.4e '));
fprintf('  v: need delta_v(i) < h_v(i)  =>  delta_v(1:2) < [%s]\n', ...
    num2str(h_v(1:2)', '%.4e '));

fprintf('=======================================\n\n');
end

function s = ternary_str(c,a,b), if c, s=a; else, s=b; end, end
