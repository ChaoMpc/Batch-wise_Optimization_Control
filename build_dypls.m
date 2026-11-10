function Model = build_dypls(Identification_data, Params)
% =========================================================================
% build_dypls -- Build a decoupled Dynamic PLS-ARX (DyPLS-ARX) model.
%
% MIMO ARX REVISION:
%   The previous per-channel ARX (u_a(k) = f(u_a hist, t_a hist)) ignored
%   cross-channel coupling (G12, G21), which is the dominant source of the
%   1-step model residual eta. The new regressor is
%       u_a(k) = sum_b [sum_j alpha_ab,j u_b(k-j) + beta_ab,j t_b(k-j)]
%   i.e., each output channel depends on ALL channels' u and t histories.
%   State-space structure (per channel block, 2n-1 dims):
%       x_a = [u_a(k-1..k-n), t_a(k-2..k-n)]  (t_a(k-1) is the input)
%   Total state dim A*(2n-1) and input dim A are UNCHANGED, so downstream
%   lifting / tube / MPC files need no modification.
% =========================================================================

    X0 = Identification_data.X0;
    Y0 = Identification_data.Y0;

    A_lv  = Params.A_lv;
    n_arx = Params.n_arx;
    N     = Params.N;
    n_out = Params.n_out;

    if isfield(Params, 'eps_Q') && ~isempty(Params.eps_Q), eps_Q = Params.eps_Q; else, eps_Q = 1e-6; end
    if isfield(Params, 'plot_model') && ~isempty(Params.plot_model), plot_model = logical(Params.plot_model); else, plot_model = false; end

    %% 1. Consistency checks
    if isempty(X0) || isempty(Y0), error('build_dypls: X0 or Y0 empty.'); end
    if size(X0,1) ~= size(Y0,1), error('build_dypls: X0/Y0 row count mismatch.'); end
    if size(X0,1) ~= N
        warning(['build_dypls: N mismatch; using actual length %d.'], size(X0,1));
        N = size(X0,1);
    end
    if size(Y0,2) ~= n_out, error('build_dypls: Y0 columns != n_out.'); end
    if A_lv <= 0 || floor(A_lv) ~= A_lv, error('build_dypls: A_lv must be pos int.'); end
    if n_arx <= 0 || floor(n_arx) ~= n_arx, error('build_dypls: n_arx must be pos int.'); end
    if N <= n_arx, error('build_dypls: N <= n_arx.'); end
    if A_lv > n_out, error('build_dypls: A_lv > n_out.'); end

    n_x = size(X0, 2);

    %% 2. NIPALS with U-deflation (unchanged)
    T = zeros(N, A_lv); U = zeros(N, A_lv);
    W = zeros(n_x, A_lv); P = zeros(n_x, A_lv);
    C = zeros(n_out, A_lv); Q = zeros(n_out, A_lv);
    X_res = X0; Y_res = Y0;
    max_iter = 3000; tol = 1e-6; eps_val = 1e-12;

    for a = 1:A_lv
        if norm(X_res,'fro') < eps_val || norm(Y_res,'fro') < eps_val
            warning('build_dypls: residual ~0 at LV %d.', a);
        end
        col_norms = sqrt(sum(Y_res.^2,1));
        [max_col_norm, init_col] = max(col_norms);
        if isempty(max_col_norm) || max_col_norm < eps_val
            u = zeros(N,1); u(1)=1;
        else
            u = Y_res(:, init_col);
        end
        err = inf; iter = 0; c = zeros(n_out,1);
        while err > tol && iter < max_iter
            iter = iter + 1;
            w = X_res' * u; if norm(w)<eps_val, break; end; w = w/norm(w);
            t = X_res * w; if norm(t)<eps_val, break; end;
            c = Y_res' * t; if norm(c)<eps_val, break; end; c = c/norm(c);
            u_new = Y_res * c; if norm(u_new)<eps_val, break; end;
            if u_new'*u < 0, u_new = -u_new; c = -c; end
            err = norm(u_new-u)/max(norm(u),eps_val); u = u_new;
        end
        w = X_res' * u;
        if norm(w) < eps_val
            w = zeros(n_x,1); t = zeros(N,1); p = zeros(n_x,1); q = zeros(n_out,1);
        else
            w = w/norm(w); t = X_res*w;
            if norm(t) < eps_val
                p = zeros(n_x,1); q = zeros(n_out,1);
            else
                p = X_res'*t/(t'*t);
                q = Y_res'*u/max(u'*u,eps_val);
            end
        end
        T(:,a)=t; U(:,a)=u; W(:,a)=w; P(:,a)=p; C(:,a)=c; Q(:,a)=q;
        if norm(t)>eps_val, X_res = X_res - t*p'; end
        if norm(u)>eps_val, Y_res = Y_res - u*q'; end
    end

    %% 3. Q conditioning
    sv_Q = svd(Q); smin = min(sv_Q); smax = max(sv_Q);
    if smin < eps_Q*max(smax,eps_val)
        warning('build_dypls: Q ill-conditioned: cond(Q)=%.3e.', smax/smin);
    end

    %% 4. Projections
    R_y = Y0 - U * Q';
    PtW = P'*W;
    if rcond(PtW) < 1e-12, R_x = W*pinv(PtW); else, R_x = W/PtW; end
    QTC = Q'*C;
    R_u = C/QTC;

    %% 5. MIMO ARX identification + state-space realization
    % theta layout (length 2*A_lv*n_arx):
    %   [ alpha_11..alpha_1n, alpha_21..alpha_2n, ..., alpha_A1..alpha_An,
    %     beta_11..beta_1n,   beta_21..beta_2n,   ..., beta_A1..beta_An ]
    % theta((b-1)*n+j)        = alpha_ab,j (u_a(k) on u_b(k-j))
    % theta(A*n + (b-1)*n+j)  = beta_ab,j  (u_a(k) on t_b(k-j))
    n_reg = 2 * A_lv * n_arx;
    n_state = 2*n_arx - 1;
    nx = A_lv * n_state;

    if isfield(Params,'arx_lambda') && ~isempty(Params.arx_lambda)
        lambda_base = Params.arx_lambda;
    else
        lambda_base = 1e-6;
    end

    Theta_cell = cell(A_lv, 1);
    % We build the full A_bar dense (not blkdiag), B_bar (nx x A_lv),
    % C_bar (A_lv x nx), E_bar (nx x A_lv).

    A_bar = zeros(nx, nx);
    B_bar = zeros(nx, A_lv);
    C_bar = zeros(A_lv, nx);
    E_bar = zeros(nx, A_lv);

    % --- First pass: fit theta for each output channel a ---
    for a = 1:A_lv
        n_samp = N - n_arx;
        Phi = zeros(n_samp, n_reg);
        Y_arx = zeros(n_samp, 1);
        row = 0;
        for k = (n_arx+1):N
            row = row + 1;
            phi_row = zeros(1, n_reg);
            for b = 1:A_lv
                u_lag = U(k-1:-1:k-n_arx, b);  % u_b(k-1..k-n)
                t_lag = T(k-1:-1:k-n_arx, b);  % t_b(k-1..k-n)
                phi_row((b-1)*n_arx + (1:n_arx)) = u_lag';
                phi_row(A_lv*n_arx + (b-1)*n_arx + (1:n_arx)) = t_lag';
            end
            Phi(row,:) = phi_row;
            Y_arx(row) = U(k, a);
        end
        PhiTPhi = Phi'*Phi;
        reg_scale = trace(PhiTPhi) / n_reg;
        if reg_scale < eps_val, reg_scale = 1; end
        theta = (PhiTPhi + lambda_base * reg_scale * eye(n_reg)) \ (Phi' * Y_arx);
        if any(~isfinite(theta))
            error('build_dypls: theta not finite at LV %d.', a);
        end
        Theta_cell{a} = theta;

        % --- Fill A_bar, B_bar, C_bar, E_bar for this output channel ---
        idx_a = (a-1)*n_state + (1:n_state);

        for b = 1:A_lv
            idx_b = (b-1)*n_state + (1:n_state);
            alpha_ab = theta((b-1)*n_arx + (1:n_arx));
            beta_ab  = theta(A_lv*n_arx + (b-1)*n_arx + (1:n_arx));

            % u_a(k) on u_b(k-1..k-n):  A_bar(idx_a(1), idx_b(1:n)) = alpha_ab'
            A_bar(idx_a(1), idx_b(1:n_arx)) = alpha_ab';

            % u_a(k) on t_b(k-2..k-n): state positions n+1..2n-1 hold t_b(k-2..k-n)
            A_bar(idx_a(1), idx_b(n_arx+1:2*n_arx-1)) = beta_ab(2:end)';

            % u_a(k) on t_b(k-1): B_bar(idx_a(1), b) = beta_ab(1)
            B_bar(idx_a(1), b) = beta_ab(1);
        end

        % Output: u_a(k) is the first state position (after one-step shift)
        C_bar(a, idx_a(1)) = 1;
        E_bar(idx_a(1), a) = 1;
    end

    % --- Shift registers (per-channel, same as before) ---
    for a = 1:A_lv
        idx_a = (a-1)*n_state + (1:n_state);
        if n_arx > 1
            A_bar(idx_a(2:n_arx), idx_a(1:n_arx-1)) = eye(n_arx-1);
        end
        if n_arx > 2
            A_bar(idx_a(n_arx+2:2*n_arx-1), idx_a(n_arx+1:2*n_arx-2)) = eye(n_arx-2);
        end
        B_bar(idx_a(n_arx+1), a) = 1;  % new t_a(k-1) enters the t_a shift register
    end

    %% 6. Augmented incremental model (unchanged formulas)
    Model.barA = A_bar;
    Model.barB = B_bar;
    Model.barC = C_bar;
    Model.barE = E_bar;

    nu = A_lv;
    Model.A_aug  = [ A_bar,       B_bar;
                     zeros(nu, nx), eye(nu) ];
    Model.B_aug  = [ B_bar;
                     eye(nu)       ];
    Model.C_aug1 = [ C_bar,       zeros(nu, nu) ];
    Model.C_aug2 = [ zeros(nu,nx), eye(nu)       ];

    %% 7. Store PLS / ARX
    Model.T = T; Model.U = U; Model.W = W; Model.P = P;
    Model.C = C; Model.Q = Q;
    Model.R_x = R_x; Model.R_u = R_u; Model.R_y = R_y;
    Model.Theta = Theta_cell;

    Model.X_res_final = X_res;
    Model.Y_res_final = Y_res;

    Model.N = N; Model.n_arx = n_arx; Model.A_lv = A_lv; Model.n_out = n_out;
    Model.n_state = n_state;

    %% 8. Standardization / steady-state
    if isfield(Identification_data,'X_mean_dev')
        Model.X_mean_dev = Identification_data.X_mean_dev;
        Model.Y_mean_dev = Identification_data.Y_mean_dev;
        Model.X_std = Identification_data.X_std;
        Model.Y_std = Identification_data.Y_std;
    else
        Model.X_mean_dev = Identification_data.X_mean;
        Model.Y_mean_dev = Identification_data.Y_mean;
        Model.X_std = Identification_data.X_std;
        Model.Y_std = Identification_data.Y_std;
    end
    if isfield(Identification_data,'X_mean_abs')
        Model.bar_v = Identification_data.X_mean_abs;
        Model.bar_y = Identification_data.Y_mean_abs;
    else
        Model.bar_v = Identification_data.u_ss + Model.X_mean_dev;
        Model.bar_y = Identification_data.y_ss + Model.Y_mean_dev;
    end
    Model.u_ss = Identification_data.u_ss;
    Model.y_ss = Identification_data.y_ss;
    Model.S_v = Model.X_std;
    Model.S_y = Model.Y_std;

    %% 9. Model fidelity check
    Y0_simulation = latent_model_simulation(N, A_lv, n_arx, Theta_cell, T, U, Q);
    if plot_model
        figure(99); clf;
        plot(0:N-1, Y0_simulation, 'LineWidth', 1.2); hold on;
        plot(0:N-1, Y0, '--', 'LineWidth', 1.2);
        xlabel('Sample k'); ylabel('Standardized deviation output');
        title('DyPLS-MIMO-ARX model vs. true output (open-loop)');
        legend('Model output', 'True output'); grid on; drawnow;
    end

    %% 10. Summary
    fprintf('build_dypls: DyPLS-MIMO-ARX modeling completed (U-deflation, deviation form).\n');
    fprintf('  Latent channels A_lv   : %d\n', A_lv);
    fprintf('  ARX order n_arx        : %d\n', n_arx);
    fprintf('  Per-channel state dim  : %d\n', n_state);
    fprintf('  Total latent state dim : %d\n', nx);
    fprintf('  ARX regressor dim      : %d  (2*A*n = %d)\n', n_reg, 2*A_lv*n_arx);
    fprintf('  ||R_y||_F              : %.4e\n', norm(R_y,'fro'));
    fprintf('  cond(Q''C)              : %.4e\n', cond(QTC));
    fprintf('  cond(P''W)              : %.4e\n', cond(PtW));
end

% ========================================================================
function Y0_simulation = latent_model_simulation(N, A_lv, n_arx, Theta_cell, T, U, Q)
    % MIMO ARX free-run
    U_sim = zeros(N, A_lv);
    U_sim(1:n_arx, :) = U(1:n_arx, :);
    for k = (n_arx+1):N
        for a = 1:A_lv
            theta = Theta_cell{a};
            u_pred = 0;
            for b = 1:A_lv
                idx_u = (b-1)*n_arx + (1:n_arx);
                idx_t = A_lv*n_arx + idx_u;
                u_lag = U_sim(k-1:-1:k-n_arx, b);
                t_lag = T(k-1:-1:k-n_arx, b);
                u_pred = u_pred + theta(idx_u)'*u_lag + theta(idx_t)'*t_lag;
            end
            U_sim(k, a) = u_pred;
        end
    end
    Y0_simulation = U_sim * Q';
end
