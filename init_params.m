function Params = init_params()
% =========================================================================
% init_params -- Wood-Berry DyPLS-IL-MPC parameters
% =========================================================================

if exist('OCTAVE_VERSION', 'builtin'), pkg load control; end

%% 1. Continuous-time plant
[numP_1, denP_1] = pade2(1);
[numP_3, denP_3] = pade2(3);
[numP_7, denP_7] = pade2(7);

G11 = tf( 12.8 * numP_1, conv([16.7 1], denP_1));
G12 = tf(-18.9 * numP_3, conv([21.0 1], denP_3));
G21 = tf(  6.6 * numP_7, conv([10.9 1], denP_7));
G22 = tf(-19.4 * numP_3, conv([14.4 1], denP_3));
G   = [G11, G12; G21, G22];
Gd  = [tf(0.5, [10 1]); tf(0.5, [10 1])];

%% 2. Discretization
Params.Ts      = 1;
Params.sys_d   = c2d(G,  Params.Ts);
Params.sys_d_d = c2d(Gd, Params.Ts);

%% 3. Data generation
Params.N             = 4000;
Params.N_val         = 2000;
Params.N_test        = 2000;
Params.rng_seed_val  = 99;
Params.rng_seed_test = 999;
Params.Ble           = 60;
Params.N_batches     = 15;
Params.verbose_ilmpc = 1;

%% 4. System dimensions
Params.n_in  = 2;
Params.n_out = 2;

%% 5. DyPLS parameters
Params.A_lv       = 2;
Params.n_arx      = 8;
Params.arx_lambda = 1e-2;
Params.eps_Q      = 1e-6;   % rank / conditioning threshold for Q

%% 6. Steady-state operating point
Params.u_ss = [1.95; 1.71];
Params.y_ss = [96.25; 0.50];

%% 7. Physical constraints
Params.u_phys_min = [1.8525; 1.6245];
Params.u_phys_max = [2.0475; 1.7955];
Params.y_phys_min = [95.0; 0.00];
Params.y_phys_max = [97.5; 1.00];

Params.u_dev_min = Params.u_phys_min - Params.u_ss;
Params.u_dev_max = Params.u_phys_max - Params.u_ss;
Params.y_dev_min = Params.y_phys_min - Params.y_ss;
Params.y_dev_max = Params.y_phys_max - Params.y_ss;

Params.H_v = [ eye(Params.n_in); -eye(Params.n_in) ];
Params.H_y = [ eye(Params.n_out); -eye(Params.n_out) ];
Params.h_v_phys = [Params.u_dev_max; -Params.u_dev_min];
Params.h_y_phys = [Params.y_dev_max; -Params.y_dev_min];

%% 8. Reference points
Params.Kp         = [20, 40, 60];
Params.y_ref_phys = [96.5, 96.7, 96.6];
Params.Ref_points = Params.y_ref_phys(:) - Params.y_ss(1);

%% 8. Reference points
%Params.Kp         = [10, 20, 30, 40, 50, 60];   % 从 [20,40,60] 增加到 6 个
%Params.y_ref_phys = [96.4, 96.5, 96.6, 96.7, 96.6, 96.5];  % 有升有降
%Params.Ref_points = Params.y_ref_phys(:) - Params.y_ss(1);

%% 9. Cost weights
Params.Q_weight = 2.0 * eye(Params.n_out);
Params.R_weight = 0.1 * eye(Params.A_lv);

%% 10. Horizons
Params.P_horizon = Params.Ble;
Params.M_horizon = Params.Ble;

%% 11. Incremental constraints
Params.delta_t_max = 0.2;
Params.ddt_max     = 0.2;

%% 12. PRBS excitation
Params.prbs_amp   = 0.05 * Params.u_ss;
Params.prbs_bits  = 10;
Params.prbs_clock = 3;
Params.noise_var  = 1e-4;

%% 13. Uncertainty-set estimation
Params.safety_factor = 1.0;
Params.q_level       = 0.95;
Params.unc_method    = 'quantile';

%% 14. Tube feedback (structured LQR, N-step iteration)
Params.tube_beta    = 1.0;
Params.tube_use_lqr = true;    % ← 保留反馈
Params.lqr_Q_scale  = 1.0;     % 结构化 Q 的统一缩放
Params.lqr_R_scale  = 1.0;
Params.tube_N       = 60;      % 可选；不写则用 Params.Ble
Params.init_R_weight_scale = 20.0;   % Initial_Feasible 的 R 放大倍数

%% 15. Tests
Params.run_init_state_test  = true;
Params.init_state_variation = 0.15;
Params.run_mc_test          = true;
Params.n_mc_runs            = 5;
Params.run_model_mismatch   = false;
Params.mismatch_gain        = 1.05;

%% 16. Plotting
Params.plot_excitation = false;
Params.plot_model      = false;
Params.plot_results    = true;
Params.rng_seed        = 26;
end

function [num, den] = pade2(tau)
    num = [ tau^2, -6 * tau, 12 ];
    den = [ tau^2,  6 * tau, 12 ];
end
