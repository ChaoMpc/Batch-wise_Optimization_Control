function [y_meas, x_next, y_true] = ...
    cryst_plant(x, v, Params, is_noisy, noise_sample)
% =========================================================================
% cryst_plant
%
% Nonlinear crystallization plant.
%
% State:
%   x = [T; Tj; C; mu0; mu1; mu2; mu3]
%
% Input:
%   v = [Tj_sp; Fc; Fanti; Fsol; Nstir; Fseed]
%
% If v is empty:
%   measure the current state without integration.
%
% Otherwise:
%   hold v constant over one sampling interval;
%   integrate the physical state;
%   measure the END-of-interval state.
%
% Measurement noise does not affect the physical state integration.
%
% The physical equations and numerical protections follow the uploaded code.
% =========================================================================

    if nargin < 4 || isempty(is_noisy)
        is_noisy = true;
    end

    if nargin < 5
        noise_sample = [];
    end

    x = x(:);

    if numel(x) ~= 7 || any(~isfinite(x))
        error('cryst_plant: x must be a finite 7-element vector.');
    end

    %% Measurement-only mode
    if isempty(v)
        x_next = x;

        [y_meas, y_true] = cryst_measure( ...
            x_next, Params, is_noisy, noise_sample);
        return;
    end

    %% Validate input
    v = v(:);

    if numel(v) ~= 6 || any(~isfinite(v))
        error('cryst_plant: v must be a finite 6-element vector.');
    end

    % Fractional power expressions require physically meaningful values.
    if v(2) < 0 || v(5) < 0
        error('cryst_plant: Fc and Nstir must be nonnegative.');
    end

    if v(3)+v(4)+10 <= 0
        error('cryst_plant: invalid antisolvent/feed denominator.');
    end

    %% Integration settings
    Ts_int = Params.Ts_int;
    n_steps = round(Params.Ts/Ts_int);

    if n_steps < 1 || ...
       abs(n_steps*Ts_int-Params.Ts) > 1e-10*max(1,Params.Ts)
        error('cryst_plant: invalid integration time settings.');
    end

    %% Physical constants
    rho_c = Params.rho_c;
    k_v = Params.k_v;

    b = Params.b;
    g = Params.g;

    kb = exp(Params.ln_kb);
    kg = exp(Params.ln_kg);

    m_w = Params.m_w;
    C_pw = Params.C_pw;
    DHc = Params.DHc;

    F_c0 = Params.F_c0;
    N0 = Params.N0;

    nu_F = Params.nu_F;
    nu_N = Params.nu_N;

    V_sol = Params.V_sol;
    C_feed = Params.C_feed;

    k_anti = Params.k_anti;
    k_seed = Params.k_seed;

    tau_inv = 1/Params.tau_j;

    %% Constant input over the sampling interval
    Tj_sp = v(1);
    Fc = v(2);
    F_anti = v(3);
    F_sol = v(4);
    N_stir = v(5);
    F_seed = v(6);

    UA = Params.UA0 * ...
         (Fc/F_c0)^nu_F * ...
         (N_stir/N0)^nu_N;

    anti_ratio = F_anti/(F_anti+F_sol+10);

    dC_feed_base = F_sol*1e-3/m_w;
    seed_rate = F_seed*k_seed;

    %% Scalar states
    T = x(1);
    Tj = x(2);
    C = x(3);

    mu0 = x(4);
    mu1 = x(5);
    mu2 = x(6);
    mu3 = x(7);

    %% Explicit Euler integration
    for j = 1:n_steps
        C_star = 1.5846e-5*T^2 - 9.0567e-3*T + 1.3066;
        C_star = C_star*(1-k_anti*anti_ratio);

        dC_sat = C-C_star;

        if dC_sat > 0
            G = kg*dC_sat^g;
            B_rate = kb*dC_sat^b;

            B_rate = min(B_rate, Params.nucleation_rate_max);
        else
            G = 0;
            B_rate = 0;
        end

        dmu0 = B_rate + seed_rate;
        dmu1 = G*mu0;
        dmu2 = 2*G*mu1;
        dmu3 = 3*G*mu2;

        dC_cryst = ...
            -rho_c*k_v*V_sol*(3*G*mu2)/m_w;

        dC_dt = ...
            dC_cryst + (C_feed-C)*dC_feed_base;

        dT_dt = ...
            (UA*(Tj-T)-DHc*m_w*dC_cryst)/(m_w*C_pw);

        dTj_dt = (Tj_sp-Tj)*tau_inv;

        T = T + dT_dt*Ts_int;
        Tj = Tj + dTj_dt*Ts_int;
        C = C + dC_dt*Ts_int;

        mu0 = mu0 + dmu0*Ts_int;
        mu1 = mu1 + dmu1*Ts_int;
        mu2 = mu2 + dmu2*Ts_int;
        mu3 = mu3 + dmu3*Ts_int;

        % Preserve original numerical protections.
        C = max(C, Params.concentration_floor);

        mu0 = max(mu0, 0);
        mu1 = max(mu1, 0);
        mu2 = max(mu2, 0);
        mu3 = max(mu3, 0);
    end

    x_next = [T; Tj; C; mu0; mu1; mu2; mu3];

    if any(~isfinite(x_next))
        error('cryst_plant: numerical instability detected.');
    end

    %% Measure the state after integration
    [y_meas, y_true] = cryst_measure( ...
        x_next, Params, is_noisy, noise_sample);
end
