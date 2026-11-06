function [Y0_sim, nrmses, vafs, nrmses_1step, vafs_1step] = ...
    free_run_validation(Model, Val_data, Params, start_idx, window_len)
if nargin < 5, start_idx = []; window_len = []; end

A_lv = Model.A_lv; n_arx = Model.n_arx; n_out = Model.n_out;
Q_pls = Model.Q; R_x = Model.R_x; R_u = Model.R_u;

N_val = size(Val_data.X0, 1);
T_val = Val_data.X0 * R_x;
U_obs = Val_data.Y0 * R_u;

% --- Global free-run (MIMO ARX) ---
U_sim = zeros(N_val, A_lv);
U_sim(1:n_arx, :) = U_obs(1:n_arx, :);
for k = (n_arx+1):N_val
    for a = 1:A_lv
        theta_a = Model.Theta{a}(:);
        u_pred = 0;
        for b = 1:A_lv
            idx_u = (b-1)*n_arx + (1:n_arx);
            idx_t = A_lv*n_arx + idx_u;
            u_lag = U_sim(k-1:-1:k-n_arx, b);
            t_lag = T_val(k-1:-1:k-n_arx, b);
            u_pred = u_pred + theta_a(idx_u)'*u_lag + theta_a(idx_t)'*t_lag;
        end
        U_sim(k, a) = u_pred;
    end
end
Y0_sim_full = U_sim * Q_pls';

nrmses = zeros(1, n_out); vafs = zeros(1, n_out);
for i = 1:n_out
    err = Val_data.Y0(:,i) - Y0_sim_full(:,i);
    denom = norm(Val_data.Y0(:,i) - mean(Val_data.Y0(:,i)));
    nrmses(i) = norm(err) / max(denom, eps);
    vafs(i) = 100 * (1 - var(err) / max(var(Val_data.Y0(:,i)), eps));
end

% --- 1-step ahead (MIMO ARX) ---
U_pred_1step = zeros(N_val, A_lv);
for k = (n_arx+1):N_val
    for a = 1:A_lv
        theta_a = Model.Theta{a}(:);
        u_pred = 0;
        for b = 1:A_lv
            idx_u = (b-1)*n_arx + (1:n_arx);
            idx_t = A_lv*n_arx + idx_u;
            u_lag = U_obs(k-1:-1:k-n_arx, b);
            t_lag = T_val(k-1:-1:k-n_arx, b);
            u_pred = u_pred + theta_a(idx_u)'*u_lag + theta_a(idx_t)'*t_lag;
        end
        U_pred_1step(k, a) = u_pred;
    end
end
Y0_pred_1step = U_pred_1step * Q_pls';

nrmses_1step = zeros(1, n_out); vafs_1step = zeros(1, n_out);
for i = 1:n_out
    err_1 = Val_data.Y0(n_arx+1:end, i) - Y0_pred_1step(n_arx+1:end, i);
    y_true_1 = Val_data.Y0(n_arx+1:end, i);
    denom_1 = norm(y_true_1 - mean(y_true_1));
    nrmses_1step(i) = norm(err_1) / max(denom_1, eps);
    vafs_1step(i)   = 100 * (1 - var(err_1) / max(var(y_true_1), eps));
end

fprintf('\n  [Model Validation Metrics on Validation Data (MIMO ARX)]\n');
fprintf('  Free-Run  NRMSE: [%s], VAF: [%s] %%\n', num2str(nrmses,'%.4f '), num2str(vafs,'%.2f '));
fprintf('  1-Step    NRMSE: [%s], VAF: [%s] %%\n', num2str(nrmses_1step,'%.4f '), num2str(vafs_1step,'%.2f '));

if ~isempty(start_idx) && ~isempty(window_len)
    end_idx = min(start_idx + window_len - 1, N_val);
    Y0_sim = Y0_sim_full(start_idx:end_idx, :);
else
    Y0_sim = Y0_sim_full;
end
end
