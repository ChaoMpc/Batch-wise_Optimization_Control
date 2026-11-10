function Params = init_params()
% =========================================================================
% init_params
%
% Crystallization Example 2: physical parameters and simulation settings.
%
% This function only defines parameters.
% Reference trajectories are generated separately by generate_reference().
%
% Time convention:
%   V(:,j)   acts over [(j-1)*Ts, j*Ts].
%   Y(:,k+1) is the output at k*Ts.
% =========================================================================

    Params = struct();

    %% 1. Project version
    Params.code_version = 'cryst_modular_v1';
    Params.parameter_source = 'uploaded crystallization code';

    %% 2. Random seeds
    Params.rng_seed = 2024;
    Params.rng_seed_val = 2025;

    % Separate calibration and online experiment seeds.
    Params.rng_seed_unc = 2026;
    Params.rng_seed_online = 2024;

    %% 3. Time and batch settings
    Params.Ts = 1.0;
    Params.Ts_int = 0.005;

    Params.Ble = 120;
    Params.ctrl_start = 60;

    Params.N_train_batches = 300;
    Params.N_val_batches = 5;
    Params.N_unc_batches = 5;
    Params.N_batches = 20;

    Params.N = Params.N_train_batches * Params.Ble;
    Params.N_val = Params.N_val_batches * Params.Ble;
    Params.N_ctrl = Params.Ble - Params.ctrl_start;

    %% 4. Physical dimensions and channel names
    Params.n_in = 6;
    Params.n_out = 6;

    Params.n_in_opt = 5;
    Params.opt_input_idx = (1:5)';

    Params.input_names = { ...
        'Tj_sp', 'Fc', 'Fanti', 'Fsol', 'Nstir', 'Fseed' ...
    };

    Params.output_names = { ...
        'T1', 'T2', 'T3', 'C_ATR', 'C_NIR', 'C_rho' ...
    };

    %% 5. Physical input and output bounds
    Params.u_phys_min = [ ...
        250.0;
         10.0;
          0.0;
          0.0;
         50.0;
          0.0 ...
    ];

    Params.u_phys_max = [ ...
        330.0;
        100.0;
         10.0;
         10.0;
        500.0;
          2.0 ...
    ];

    % Preserve the original code's numerical bounds.
    % Only C_ATR has an active engineering safety range.
    Params.y_phys_min = [ ...
        -1e6;
        -1e6;
        -1e6;
         0.01;
        -1e6;
        -1e6 ...
    ];

    Params.y_phys_max = [ ...
         1e6;
         1e6;
         1e6;
         0.04;
         1e6;
         1e6 ...
    ];

    Params.H_v = [ ...
        eye(Params.n_in_opt);
       -eye(Params.n_in_opt) ...
    ];

    Params.h_v = [ ...
        Params.u_phys_max(Params.opt_input_idx);
       -Params.u_phys_min(Params.opt_input_idx) ...
    ];

    Params.H_y = [eye(Params.n_out); -eye(Params.n_out)];

    Params.h_y = [ ...
        Params.y_phys_max;
       -Params.y_phys_min ...
    ];

    %% 6. Physical crystallization parameters
    Params.rho_c = 1296;
    Params.k_v = 0.24;

    Params.b = 6.2;
    Params.g = 1.5;

    Params.ln_kb = 35.0;
    Params.ln_kg = -2.5;

    Params.m_w = 0.45;
    Params.C_pw = 4.18;
    Params.DHc = 20e3;

    Params.UA0 = 0.12 * 60;
    Params.tau_j = 2.0;

    Params.F_c0 = 1.0;
    Params.N0 = 400;

    Params.nu_F = 0.5;
    Params.nu_N = 0.2;

    Params.V_sol = 0.5e-3;
    Params.C_feed = 0.03;

    Params.k_anti = 0.8;
    Params.k_seed = 5e-4;

    %% 7. Numerical protections inherited from the original plant
    Params.nucleation_rate_max = 1e6;
    Params.concentration_floor = 1e-5;

    %% 8. Physical initial state
    % x = [T; Tj; C; mu0; mu1; mu2; mu3]
    Params.x0_nom = [ ...
        305.0;
        305.0;
        0.0256;
        1e-4;
        1e-3;
        1e-2;
        1e-3 ...
    ];

    %% 9. Baseline input profile
    % These values reproduce the original init_params.m.
    Params.baseline.Tj_initial = 305.0;
    Params.baseline.Tj_drop = 10.0;

    Params.baseline.Fc = 60.0;

    Params.baseline.Fanti_initial = 5.0;%4
    Params.baseline.Fanti_increase = 1.0;

    Params.baseline.Fsol = 5.0; %5
    Params.baseline.Nstir = 300.0;

    Params.baseline.Fseed_initial = 1.0;
    Params.baseline.Fseed_steps = 5;
    Params.baseline.Fseed_after = 0.0;

    %% 10. Measurement noise
    % These are UNIFORM noise half-widths, not Gaussian standard deviations.
    Params.measurement_noise_halfwidth = [ ...
        0.015;    % 原 0.03
        0.020;    % 原 0.04
        0.025;    % 原 0.05
        1e-4;     % 原 2e-4
        1.5e-4;   % 原 3e-4
        1.25e-4   % 原 2.5e-4
    ];
    Params.measurement_noise_type = 'uniform';

    %% 11. Identification and model selection
    % Preserve the original code's default model dimensions.
    Params.A_lv = 2;
    Params.n_arx = 2;
    Params.arx_lambda = 1e-3;

    Params.A_lv_cv_range = [2];
    Params.n_arx_cv_range = [2, 3, 4];

    Params.nipals_tol = 1e-9;
    Params.nipals_max_iter = 500;

    Params.pls_rank_tol = 1e-10;
    Params.pls_sigma_min = 1e-8;

    % Preserve the current scaling configuration.
    Params.W_y = ones(Params.n_out, 1);

    % Excitation amplitude:
    % the original generator uses 50% of the available distance to bounds.
    Params.excitation_room_fraction = 0.5;
    Params.excitation_amplitude_fraction_min = 0.2;

    Params.prbs_periods = [5; 5; 8; 8; 5; 8];

    % Compatibility for auxiliary open-loop step experiments.
    % The main adaptive excitation generator does not use this field.
    Params.prbs_amp = ...
        Params.excitation_room_fraction * ...
        min( ...
            [305; 60; 4; 5; 300; 1] - Params.u_phys_min, ...
            Params.u_phys_max - [305; 60; 4; 5; 300; 1]);

    %% 12. Point-to-point tracking
    Params.track_idx = 4;
    Params.Kp = [90, 120];

    Params.y_ref_phys = [0.022; 0.015];

    Params.Q_weight = 100 * eye(Params.n_out);
    Params.R_weight = 1e-8;

    Params.delta_t_max = 0.2;
    Params.ddt_max = 0.5;

    %% 13. Uncertainty and Tube configuration
    Params.safety_factor = 1.0;
    Params.q_level = 1.0;

    Params.Q_lqr = 1;
    Params.R_lqr = 1;

    % Future Tube stage uses a full MIMO invariant ellipsoid.
    % Alpha must satisfy sqrt(1+alpha)*rho(A_K) < 1.
    Params.tube_beta = 0.6;
    Params.tube_alpha_grid = [0.01, 0.03, 0.1, 0.3, 1.0];
    Params.tube_beta_grid = [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0, 100.0]; % 【新增】
    Params.tube_shape_regularization = 1e-12;

    %% 14. Online initial-state deviations
    % The old code contains enable_init_perturb=true but does not actually
    % implement an initial-state perturbation in run_ilmpc.
    % Preserve the baseline experiment by specifying zero amplitude.
    Params.enable_init_perturb = true;
    Params.init_perturb_halfwidth = zeros(7, 1);

    %% 15. Injected concentration disturbance
    % Preserve the actual uploaded controller implementation.
    Params.disturbance.enabled = true;
    Params.disturbance.batches = [15, 16];
    Params.disturbance.input_columns = 80:85;
    Params.disturbance.state_idx = 3;
    Params.disturbance.increment = 0.0002;

    %% 16. Comparison method
    Params.xue_reference_lambda = -0.6;

% --- Oh (文献12) PTP-ILMPC 参数 ---
Params.oh_prediction_horizon = 60; % 预测时域 p
Params.oh_control_horizon = 60;    % 控制时域 m
Params.Q_weight = 100;             % 追踪权重
Params.S_weight = 0.0;             % 绝对输入惩罚 u'Su
Params.R_weight = 1e-4;            % 增量惩罚 delta_u'R delta_u
Params.P_weight = 1e-4;            % 批次增量惩罚 Delta_u'P Delta_u
Params.E_weight = 1e2;             % 松弛变量惩罚

    %% 17. QP numerical settings
    Params.qp_max_iter = 5000;
    Params.qp_tol_x = 1e-10;
    Params.qp_feas_tol = 1e-8;

    %% 18. Plotting
    Params.plot_results = true;
    Params.compare_batches = [15, 16, 17, 20];

    %% 19. Basic consistency checks
    if abs(Params.Ts/Params.Ts_int - ...
           round(Params.Ts/Params.Ts_int)) > 1e-10
        error('init_params: Ts must be an integer multiple of Ts_int.');
    end

    if Params.ctrl_start < Params.n_arx || ...
       Params.ctrl_start >= Params.Ble
        error('init_params: invalid control start.');
    end

    if any(Params.u_phys_min >= Params.u_phys_max) || ...
       any(Params.y_phys_min >= Params.y_phys_max)
        error('init_params: invalid physical bounds.');
    end

    if numel(Params.Kp) ~= numel(Params.y_ref_phys) || ...
       any(Params.Kp <= Params.ctrl_start) || ...
       any(Params.Kp > Params.Ble)
        error('init_params: invalid tracking specification.');
    end

    fprintf('init_params: physical parameters initialized.\n');
    fprintf('  Batch length: %d steps; control segment: %d steps.\n', ...
        Params.Ble, Params.N_ctrl);
    fprintf('  Physical inputs: %d; optimized inputs: %d; outputs: %d.\n', ...
        Params.n_in, Params.n_in_opt, Params.n_out);
end
