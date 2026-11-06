function M = compute_metrics(Data, Params, Model, UncSets)
% =========================================================================
% compute_metrics -- Aggregate closed-loop metrics:
%   constraint violations, empirical coverage of eta and r, QP solve times.
%
% Coverage flags < 0 mean "not evaluated" and are excluded. When no
% samples are available (e.g., r is not evaluated because A_lv == n_out),
% coverage is reported as NaN and n_samples = 0.
% =========================================================================

N_batches = Params.N_batches;
N = Params.Ble;
H_v = Params.H_v; H_y = Params.H_y;
h_v_phys = Params.h_v_phys; h_y_phys = Params.h_y_phys;

M.max_violation_v = 0; M.max_violation_y = 0;
M.n_violations_v  = 0; M.n_violations_y  = 0;

for b = 1:N_batches
    for k = 1:N                                    % <-- 由 1:N-1 改为 1:N
        v = Data.V_hist(:, k, b);
        viol_v = H_v*v - h_v_phys;
        M.n_violations_v  = M.n_violations_v + sum(viol_v > 1e-6);
        M.max_violation_v = max(M.max_violation_v, max(viol_v));
    end
    for k = 1:N
        y = Data.Y_hist(:, k, b);
        viol_y = H_y*y - h_y_phys;
        M.n_violations_y  = M.n_violations_y + sum(viol_y > 1e-6);
        M.max_violation_y = max(M.max_violation_y, max(viol_y));
    end
end

% ---- eta coverage ----
flags = Data.eta_in_flags(:); flags = flags(flags>=0);
M.n_eta_samples = numel(flags);
if M.n_eta_samples == 0
    M.coverage_eta = NaN;
else
    M.coverage_eta = mean(flags);
end

% ---- r coverage (may be N/A when A_lv == n_out) ----
flags = Data.r_in_flags(:); flags = flags(flags>=0);
M.n_r_samples = numel(flags);
if M.n_r_samples == 0
    M.coverage_r = NaN;
else
    M.coverage_r = mean(flags);
end

M.qp_times = Data.QP_times(:);

fprintf('compute_metrics: %d/%d input violations, %d/%d output violations\n', ...
    M.n_violations_v, N_batches*(N-1), M.n_violations_y, N_batches*(N-1));

if isnan(M.coverage_eta)
    fprintf('  eta coverage: N/A (no samples)\n');
else
    fprintf('  eta coverage: %.2f%% (%d samples)\n', ...
        100*M.coverage_eta, M.n_eta_samples);
end
if isnan(M.coverage_r)
    fprintf('  r   coverage: N/A (A_lv == n_out, r == 0)\n');
else
    fprintf('  r   coverage: %.2f%% (%d samples)\n', ...
        100*M.coverage_r, M.n_r_samples);
end
end
