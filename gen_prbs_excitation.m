function [Identification_data, Params] = gen_prbs_excitation(Params)
% =========================================================================
% gen_prbs_excitation -- Generate multi-channel PRBS excitation data.
%
% Coordinate conventions (IMPORTANT):
%   X_dev = X_abs - u_ss              (deviation input)
%   Y_dev = Y_abs - y_ss              (deviation output)
%   X0    = (X_dev - X_mean_dev)/X_std
%   Y0    = (Y_dev - Y_mean_dev)/Y_std
%
% Stored statistics:
%   X_mean_dev, Y_mean_dev : means of the DEVIATION signals
%   X_std,      Y_std      : std of the DEVIATION signals
%   X_mean_abs = u_ss + X_mean_dev : absolute mean (paper's bar_v)
%   Y_mean_abs = y_ss + Y_mean_dev : absolute mean (paper's bar_y)
%
% Online absolute form (paper notation):
%   t       = R_x' * S_v^{-1} * (v_abs - bar_v)
%   v_abs   = bar_v + S_v * P * t
%   y_abs   = bar_y + S_y * Q * u + S_y * r
% =========================================================================

    N     = Params.N;
    Ts    = Params.Ts;
    n_in  = Params.n_in;
    n_out = Params.n_out;

    u_ss = Params.u_ss(:);
    y_ss = Params.y_ss(:);

    prbs_amp   = Params.prbs_amp(:);
    prbs_bits  = Params.prbs_bits;
    prbs_clock = Params.prbs_clock;

    %% 1. Random seed
    if isfield(Params, 'rng_seed') && ~isempty(Params.rng_seed)
        if exist('rng', 'builtin')
            rng(Params.rng_seed);
        else
            rand('state',  Params.rng_seed);
            randn('state', Params.rng_seed);
        end
        base_seed = Params.rng_seed;
    else
        base_seed = 1;
    end

    %% 2. Validation
    [sys_ny, sys_nu] = size(Params.sys_d);
    if sys_nu ~= n_in || sys_ny ~= n_out
        error('gen_prbs_excitation: sys_d size inconsistent.');
    end

    %% 3. PRBS excitation
    X_dev = zeros(N, n_in);
    for i = 1:n_in
        amp_i = prbs_amp(i);
        X_dev(:, i) = prbs_signal(N, prbs_bits, -amp_i, amp_i, ...
                                  prbs_clock, base_seed + i);
    end

    %% 4. Disturbance input (truncated to [-3, 3] sigma)
    noise_std = sqrt(Params.noise_var);
    d_raw = randn(N, 1);
    d_raw = max(min(d_raw, 3.0), -3.0);   % was 2.0
    d = noise_std * d_raw;

    %% 5. Simulate the INCREMENTAL plant (deviation model)
    time = (0:N - 1)' * Ts;
    y_dev_nom  = lsim(Params.sys_d,   X_dev, time);
    y_dev_dist = lsim(Params.sys_d_d, d,     time);
    y_dev      = y_dev_nom + y_dev_dist;

    %% 6. Absolute physical trajectories
    X_abs = X_dev + ones(N, 1) * u_ss';
    Y_abs = y_dev + ones(N, 1) * y_ss';

    %% 7. Standardization statistics (computed on DEVIATION signals)
    X_mean_dev = mean(X_dev, 1)';
    X_std      = std(X_dev, 0, 1)';
    Y_mean_dev = mean(y_dev, 1)';
    Y_std      = std(y_dev, 0, 1)';

    X_std(X_std < 1e-12) = 1;
    Y_std(Y_std < 1e-12) = 1;

    X0 = (X_dev - ones(N, 1) * X_mean_dev') ./ (ones(N, 1) * X_std');
    Y0 = (y_dev - ones(N, 1) * Y_mean_dev') ./ (ones(N, 1) * Y_std');

    %% 8. Absolute means (paper's bar_v, bar_y)
    X_mean_abs = u_ss + X_mean_dev;
    Y_mean_abs = y_ss + Y_mean_dev;

    %% 9. Package
    Identification_data.X          = X_dev;
    Identification_data.Y          = y_dev;
    Identification_data.X_abs      = X_abs;
    Identification_data.Y_abs      = Y_abs;
    Identification_data.X0         = X0;
    Identification_data.Y0         = Y0;

    % Deviation-space statistics
    Identification_data.X_mean     = X_mean_dev;   % legacy alias
    Identification_data.Y_mean     = Y_mean_dev;   % legacy alias
    Identification_data.X_std      = X_std;
    Identification_data.Y_std      = Y_std;
    Identification_data.X_mean_dev = X_mean_dev;
    Identification_data.Y_mean_dev = Y_mean_dev;

    % Absolute-space statistics (paper's bar_v, bar_y)
    Identification_data.X_mean_abs = X_mean_abs;
    Identification_data.Y_mean_abs = Y_mean_abs;

    % Steady-state operating point
    Identification_data.u_ss = u_ss;
    Identification_data.y_ss = y_ss;

    Identification_data.N      = N;
    Identification_data.Ts     = Ts;
    Identification_data.n_in   = n_in;
    Identification_data.n_out  = n_out;
    Identification_data.time   = time;

    fprintf(['gen_prbs_excitation: PRBS data generated ' ...
             '(Noise var=%.1e, std=%.4f, clipped 3-sigma).\n'], ...
             Params.noise_var, noise_std);
end

function u = prbs_signal(N, n_bits, amp_low, amp_high, clock_period, seed)
    switch n_bits
        case  2,  taps = [2 1]; case  3,  taps = [3 2]; case  4,  taps = [4 3];
        case  5,  taps = [5 3]; case  6,  taps = [6 5]; case  7,  taps = [7 6];
        case  8,  taps = [8 6 5 4]; case  9,  taps = [9 5]; case 10,  taps = [10 7];
        case 11,  taps = [11 9]; case 12,  taps = [12 11 10 4]; case 13,  taps = [13 12 11 8];
        case 14,  taps = [14 13 12 2]; case 15,  taps = [15 14]; case 16,  taps = [16 15 13 4];
        case 17,  taps = [17 14]; case 18,  taps = [18 11]; case 19,  taps = [19 18 17 14];
        case 20,  taps = [20 17]; case 21,  taps = [21 19]; case 22,  taps = [22 21];
        case 23,  taps = [23 18]; case 24,  taps = [24 23 22 17]; case 25,  taps = [25 22];
        case 26,  taps = [26 25 24 20]; case 27,  taps = [27 26 25 22]; case 28,  taps = [28 25];
        case 29,  taps = [29 27]; case 30,  taps = [30 29 28 7];
        otherwise, error('prbs_signal: invalid n_bits.');
    end

    register = zeros(1, n_bits);
    s = max(1, mod(seed, 2^n_bits - 1));
    for k = 1:n_bits
        register(k) = bitand(bitshift(s, -(k-1)), 1);
    end
    if all(register == 0), register(1) = 1; end

    n_bits_needed = ceil(N / clock_period) + 1;
    bits = zeros(n_bits_needed, 1);

    for k = 1:n_bits_needed
        bits(k) = register(end);
        fb = 0;
        for t = taps, fb = bitxor(fb, register(t)); end
        register(2:end) = register(1:end-1);
        register(1)     = fb;
    end

    u_bits = kron(bits, ones(clock_period, 1));
    u_bits = u_bits(1:N);
    u = amp_low + (amp_high - amp_low) * u_bits;
end
