function [Phi_x, x0_contrib, Phi_y, y0_const, Phi_v, v0_const] = ...
    build_lifting_matrices(N, nx, n_out, n_in, A_lv, barA, barB, barC, ...
        S_y, Q, S_v, P_pls, X_mean_dev, Y_mean_dev, x0)
% =========================================================================
% build_lifting_matrices -- N-step lifting matrices (deviation form).
%
% COORDINATE CONVENTION:
%   All signals below are DEVIATION quantities:
%     x = latent state (deviation from training zero state)
%     t = latent input score
%     u = latent output score
%     y_dev = y_abs - y_ss,  v_dev = v_abs - u_ss
%
% N-STEP LIFTING (aligned with Eq. (FHOCP), i in 0..N-1):
%   t_vec = col(t(0),...,t(N-1)),  size A_lv*N
%
%   x(k)     = Phi_x(k) t_vec + barA^k x0,   k = 1..N
%   y_dev(k) = Y_mean_dev + S_y Q Cbar x(k), k = 1..N
%   v_dev(k) = X_mean_dev + S_v P t(k),      k = 0..N-1
%
% X_mean_dev, Y_mean_dev are the means of the DEVIATION signals
% (≈ 0 for zero-mean PRBS excitation). They are NOT bar_v, bar_y.
% =========================================================================
    nx_ = nx;
    n_state_total = A_lv * N;   % t(0)...t(N-1)

    Phi_x   = zeros(nx_ * N, n_state_total);
    Gamma_x = zeros(nx_ * N, nx_);

    % A_pow{m} = barA^{m-1}
    A_pow = cell(N + 1, 1);
    A_pow{1} = eye(nx_);
    for i = 2:N+1
        A_pow{i} = barA * A_pow{i - 1};
    end

    % Phi_x 块 (k, j) 对应 x(k+1) ← t(j):  j = 0..k  (k = 0..N-1)
    for k = 0:N-1
        r = k*nx_+1 : (k+1)*nx_;          % x(k+1)
        for j = 0:k
            c = j*A_lv+1 : (j+1)*A_lv;    % t(j)
            Phi_x(r, c) = A_pow{k - j + 1} * barB;   % barA^{k-j}
        end
        Gamma_x(r, :) = A_pow{k + 2};      % barA^{k+1}
    end
    x0_contrib = Gamma_x * x0;

    C_y      = S_y * Q * barC;
    Phi_y    = zeros(n_out * N, n_state_total);
    y0_const = zeros(n_out * N, 1);
    for k = 1:N
        r  = (k-1)*n_out+1 : k*n_out;
        sr = (k-1)*nx_+1    : k*nx_;
        Phi_y(r, :)  = C_y * Phi_x(sr, :);
        y0_const(r)  = Y_mean_dev + C_y * x0_contrib(sr);
    end

    M_v      = S_v * P_pls;
    Phi_v    = zeros(n_in * N, n_state_total);
    v0_const = zeros(n_in * N, 1);
    for k = 0:N-1
        r = k*n_in+1 : (k+1)*n_in;
        c = k*A_lv+1 : (k+1)*A_lv;
        Phi_v(r, c) = M_v;
        v0_const(r) = X_mean_dev;
    end
end
