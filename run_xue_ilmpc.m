function Data = run_xue_ilmpc( ...
    Problem, Initial, Optimal, Scenario, Model, Tube, UncSets, Params)
% =========================================================================
% run_xue_ilmpc
%
% 严谨复现文献：
% Xue, S., Zhao, Z., & Liu, F. (2022). Latent variable point-to-point
% iterative learning model predictive control via reference trajectory updating.
%
% 特点：
%   1. Unconstrained optimization (Eq. 31 解析解)
%   2. Reference trajectory updating via \Phi matrix (Eq. 15 & 17)
%   3. Kalman observer feedback (在 simulate_control_batch 中实现)
%   4. 物理限幅防发散的 Anti-windup 机制
% =========================================================================

    B = Params.N_batches;
    N = Params.Ble;
    H = Problem.H;
    no = Params.n_out;
    nt = Problem.nt;

    Q = Params.Q_weight;
    if isscalar(Q), Q = Q*eye(no); end

    W = chol((Q+Q')/2);
    Wbar = kron(speye(H),W);

    Lfull = Wbar*Problem.Lift.Phi_y;

    % 无约束二次型 Hessian 矩阵 (对应文献 G^T Q G + R)
    Hessian = 2*(Lfull'*Lfull+Problem.Pprox);
    Hessian = full((Hessian+Hessian')/2);

    % Cholesky 分解，用于快速求解线性方程组
    [Rh,flag] = chol(Hessian);
    if flag ~= 0
        error('run_xue_ilmpc: unconstrained Hessian is not positive definite.');
    end

    % =========================================================================
    % Kalman 观测器设计 (用于批内状态反馈)
    % 对应文献中利用 Kalman filter 估计状态 \hat{x}
    % =========================================================================
    nx = Model.nx;
    n_out_lv = Model.A_lv; % 输出潜变量维度

    % 假设过程噪声与测量噪声协方差
    Q_kal = 1e-3 * eye(nx);       % 过程噪声
    R_kal = 1e-2 * eye(n_out_lv); % 测量噪声

    % 求解稳态 Kalman 增益 L
    % 使用对偶 dlqr 避免 dare 报错，保证观测器误差收敛
    try
        [K_dual, ~, ~] = dlqr(Model.barA', Model.barC', Q_kal, R_kal);
        L_kal = K_dual';
    catch
        warning('run_xue_ilmpc: dlqr 失败，退回简单增益。');
        L_kal = 0.1 * Model.barC';
    end

    % 为了不修改 simulate_control_batch 的函数接口，将 Kalman 相关矩阵塞入 Params
    Params.xue_Kalman_L = L_kal;
    Params.xue_A = Model.barA;
    Params.xue_B = Model.barB;
    Params.xue_C = Model.barC;
    Params.xue_K_feedback = Tube.K; % 文献 Eq.31 中通过反馈增益修正

    Data = struct();
    Data.method = 'Xue adapted comparison';

    Data.Y_hist = zeros(no,N+1,B);
    Data.Y_true_hist = zeros(no,N+1,B);

    Data.V_hist = zeros(Params.n_in,N,B);
    Data.V_request_hist = zeros(Params.n_in,N,B);

    Data.X_hist = zeros(7,N+1,B);

    Data.Y_nominal_hist = zeros(no,N+1,B);
    Data.V_nominal_hist = zeros(Params.n_in,N,B);

    Data.reference_hist = zeros(no,H,B);
    Data.T_hist = zeros(Model.A_lv,H,B);

    Data.QP_times = zeros(B,1);

    Data.Err_true_act = zeros(B,1);
    Data.Err_opt_act = zeros(B,1);
    Data.Err_opt_nom = zeros(B,1);

    Data.Batch = cell(B,1);

    t_previous = Initial.t;

    % 初始化参考轨迹（保证初始时刻关键点就是目标设定值）
    reference = Params.y_base(:,Problem.Lift.output_columns);
    ref_vec_init = reference(:);
    ref_vec_init(Problem.selected_rows) = Problem.selected_refs;
    reference = reshape(ref_vec_init, no, H);

    lambda = Params.xue_reference_lambda;

    for b = 1:B
        % =========================================================================
        % 1. 参考轨迹更新 (Eq. 15 & 17: \Phi 矩阵屏蔽关键点误差)
        % =========================================================================
        if b > 1
            previous_actual = Data.Y_hist(:,Problem.Lift.output_columns,b-1);
            err_latent = reference - previous_actual;

            err_vec = err_latent(:);
            err_vec(Problem.selected_rows) = 0;
            err_latent = reshape(err_vec, no, H);

            reference = reference + lambda * err_latent;
        end

        % =========================================================================
        % 2. 无约束优化解析解求解 (Eq. 31)
        % =========================================================================
        residual_previous = Wbar * ...
            (Problem.Lift.Phi_y*t_previous + ...
             Problem.Lift.y_offset-reference(:));

        f = 2*Lfull'*residual_previous;

        timer = tic();
        d = -(Rh\(Rh'\f)); % d = -H^{-1} f
        Data.QP_times(b) = toc(timer);

        if numel(d) ~= nt || any(~isfinite(d))
            error('run_xue_ilmpc: invalid unconstrained solution.');
        end

        t_new = t_previous+d;

        % =========================================================================
        % 3. Anti-windup (防止无约束优化导致 t_previous 无限发散)
        % =========================================================================
        T_mat = reshape(t_new, Model.A_lv, H);

        u_min_mat = zeros(Model.A_lv, H);
        u_max_mat = zeros(Model.A_lv, H);
        for k = 1:H
            col_idx = Problem.Lift.input_columns(k);
            vdev_min = Params.u_phys_min(Model.input_idx) - Params.u_base(Model.input_idx, col_idx);
            vdev_max = Params.u_phys_max(Model.input_idx) - Params.u_base(Model.input_idx, col_idx);

            u_min_mat(:,k) = Model.R_x' * ((vdev_min - Model.X_mean) ./ Model.X_std);
            u_max_mat(:,k) = Model.R_x' * ((vdev_max - Model.X_mean) ./ Model.X_std);
        end

        T_mat = max(min(T_mat, u_max_mat), u_min_mat);
        t_new = T_mat(:);

        Nominal = evaluate_nominal_trajectory(Problem,t_new);

        % =========================================================================
        % 4. 批内执行 (包含 Kalman 状态反馈)
        % =========================================================================
        Batch = simulate_control_batch( ...
            Nominal,'xue',Scenario,b,Model,Tube,UncSets,Params);

        Data.Y_hist(:,:,b) = Batch.Y;
        Data.Y_true_hist(:,:,b) = Batch.Y_true;
        Data.V_hist(:,:,b) = Batch.V;
        Data.V_request_hist(:,:,b) = Batch.V_request;
        Data.X_hist(:,:,b) = Batch.X;

        Ynom = Params.y_base;
        Vnom = Params.u_base;
        Ynom(:,Problem.Lift.output_columns) = Nominal.Y;
        Vnom(Model.input_idx,Problem.Lift.input_columns) = Nominal.V;

        Data.Y_nominal_hist(:,:,b) = Ynom;
        Data.V_nominal_hist(:,:,b) = Vnom;

        Data.reference_hist(:,:,b) = reference;
        Data.T_hist(:,:,b) = Nominal.T;

        Ycontrol = Batch.Y(:,Problem.Lift.output_columns);
        yvector = Ycontrol(:);
        actual_selected = yvector(Problem.selected_rows);

        Data.Err_true_act(b) = norm(actual_selected-Problem.selected_refs);
        Data.Err_opt_act(b) = norm(actual_selected-Optimal.selected_outputs);
        Data.Err_opt_nom(b) = norm(Nominal.selected_outputs-Optimal.selected_outputs);

        Data.Batch{b} = Batch;

        fprintf(['Xue batch %2d/%2d: act/ref %.4e; ', ...
                 'act/opt %.4e; nom/opt %.4e; ', ...
                 'clipped channel-samples %d.\n'], ...
            b,B,Data.Err_true_act(b), ...
            Data.Err_opt_act(b),Data.Err_opt_nom(b), ...
            sum(Batch.clip_flags(:)));

        t_previous = t_new;
    end
    Data.final_t = t_previous;
end
