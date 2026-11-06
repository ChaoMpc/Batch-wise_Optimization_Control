function Batch = simulate_control_batch( ...
    Nominal, method, Scenario, batch_index, Model, Tube, UncSets, Params)
% =========================================================================
% simulate_control_batch
%
% 共用非线性执行引擎。
%
% 说明：
%   - 'tube' 分支：实现了带物理硬约束保证的 Tube-MPC（无滤波，直接使用测量状态）
%   - 'xue'  分支：实现了无约束优化 + Kalman 观测器反馈 + 物理限幅 (文献32)
%   - 'oh'   分支：实现了 ILO 观测器反馈 + 物理限幅 (文献12)
% =========================================================================

    if ~any(strcmp(method,{'tube','xue','oh'}))
        error('simulate_control_batch: unknown method.');
    end

    N = Params.Ble;
    c0 = Params.ctrl_start;
    H = N-c0;

    ni = Params.n_in_opt;
    no = Params.n_out;
    A  = Model.A_lv;

    b = batch_index;

    Y        = zeros(no,N+1);
    Ytrue    = zeros(no,N+1);
    V        = zeros(Params.n_in,N);
    Vrequest = zeros(Params.n_in,N);
    X        = zeros(7,N+1);

    x = Scenario.X0(:,b);
    X(:,1) = x;

    [Y(:,1),Ytrue(:,1)] = cryst_measure( ...
        x,Params,true,Scenario.noise(:,1,b));

    %% 1. 固定 baseline 段
    for j = 1:c0
        x = x+Scenario.state_jump(:,j,b);
        v = Params.u_base(:,j);
        V(:,j)        = v;
        Vrequest(:,j) = v;
        [Y(:,j+1),x,Ytrue(:,j+1)] = cryst_plant( ...
            x,v,Params,true,Scenario.noise(:,j+1,b));
        X(:,j+1) = x;
    end

    History = initialize_score_history(Y,V,Model,Params);

    % === 隔离区 1：状态观测器初始状态 ===
    if strcmp(method, 'xue')
        x_hat_kal = assemble_latent_state(History, Model);
    elseif strcmp(method, 'oh')
        x_hat_oh = assemble_latent_state(History, Model);
    end

    %% 2. 控制段诊断缓存
    Trequest = zeros(A,H);
    Tapplied = zeros(A,H);

    Eta = zeros(A,H);
    R   = zeros(no,H);

    Epre  = zeros(Model.nx,H);
    Epost = zeros(Model.nx,H);

    eta_box       = false(1,H);
    r_box         = false(1,H);
    eta_ellipsoid = false(1,H);
    r_ellipsoid   = false(1,H);

    tube_pre       = zeros(1,H);
    tube_post      = zeros(1,H);
    tube_pre_obs   = zeros(1,H);
    tube_post_obs  = zeros(1,H);

    clip_flags = false(ni,H);
    violation_flags = false(ni,H);

    Lp     = chol(Tube.P,'lower');
    Lp_obs = chol(Tube.P_obs,'lower');
    Leta   = chol(UncSets.W_eta,'lower');
    Lr     = chol(UncSets.W_r,'lower');

    %% 3. 控制段
    for k = 1:H
        j = c0+k;

        % =====================================================================
        % 状态观测：完全隔离
        % =====================================================================
        if strcmp(method, 'tube')
            % 【Tube 专属】直接使用测量数据重组状态
            x_observed = assemble_latent_state(History,Model);

        elseif strcmp(method, 'xue')
            % 【Xue 专属】Kalman 滤波：测量更新 (Correction)
            y_meas_lv = History.U(:,1);
            y_pred_lv = Params.xue_C * x_hat_kal;

            x_observed = x_hat_kal + Params.xue_Kalman_L * (y_meas_lv - y_pred_lv);
            x_hat_kal = x_observed;

        elseif strcmp(method, 'oh')
            % 【Oh 专属】ILO (迭代学习观测器)：测量更新
            y_meas_lv = History.U(:,1);
            y_pred_lv = Params.oh_C * x_hat_oh;

            x_observed = x_hat_oh + Params.oh_L * (y_meas_lv - y_pred_lv);
            x_hat_oh = x_observed;
        end

        e_pre = x_observed - Nominal.X_pre(:,k);

        Epre(:,k)       = e_pre;
        tube_pre(k)     = sum((Lp\e_pre).^2);
        tube_pre_obs(k) = sum((Lp_obs\e_pre).^2);

        t_nominal = Nominal.T(:,k);

        % =====================================================================
        % 控制律：完全隔离
        % =====================================================================
        if strcmp(method,'tube')
            % 【Tube 专属】使用 Tube 增益
            t_request = t_nominal + Tube.K * e_pre;
        elseif strcmp(method,'xue')
            % 【Xue 专属】使用 Xue 的等效反馈增益
            t_request = t_nominal + Params.xue_K_feedback * e_pre;
        elseif strcmp(method,'oh')
            % 【Oh 专属】使用 Oh 的等效反馈增益
            t_request = t_nominal + Params.oh_K_feedback * e_pre;
        end

        v_request = Params.u_base(:,j);
        v_request(Model.input_idx) = ...
            Params.u_base(Model.input_idx,j) + ...
            Model.X_mean + ...
            diag(Model.X_std)*Model.P*t_request;

        lower = Params.u_phys_min;
        upper = Params.u_phys_max;

        raw_violation = max([0; v_request-upper; lower-v_request]);


        % =====================================================================
        % 饱和与越界处理：完全隔离
        % =====================================================================
        if strcmp(method,'tube')
            % 【Tube 专属】绝不妥协，越界直接报错
            if raw_violation > 1e-8
                error(['Tube input request violates physical bounds: ', ...
                       'batch %d, input column %d, violation %.3e.'], ...
                      b,j,raw_violation);
            end
            v_applied = v_request;

        elseif strcmp(method, 'xue')
            % 【Xue 专属】完全不作物理限幅，按照无约束请求执行
            v_applied = v_request;

            % 最小存活保护（防止结晶数值模型崩溃）
            if v_applied(2) < 0.1, v_applied(2) = 0.1; end % Fc 必须非负
            if v_applied(5) < 1.0, v_applied(5) = 1.0; end % Nstir 必须非负
            if v_applied(3) + v_applied(4) + 10 <= 0
                v_applied(3) = 1.0; v_applied(4) = 1.0;
            end

        else
            % 【Oh 专属】保留原本的软约束截断逻辑
            v_applied = max(min(v_request,upper),lower);
        end

        clip_flags(:,k) = ...
            abs(v_applied(Model.input_idx) - v_request(Model.input_idx)) > 1e-10;

        Vrequest(:,j) = v_request;
        V(:,j)        = v_applied;

        % === 新增：统计实际施加的输入是否越界 ===
        v_opt_applied = v_applied(Model.input_idx);
        violation_flags(:,k) = ...
            (v_opt_applied > Params.u_phys_max(Model.input_idx) + 1e-8) | ...
            (v_opt_applied < Params.u_phys_min(Model.input_idx) - 1e-8);

        vdev = v_applied(Model.input_idx) - Params.u_base(Model.input_idx,j);
        t_applied = Model.R_x' * ((vdev-Model.X_mean)./Model.X_std);

        u_prediction = Model.barC * ...
            (Model.barA*x_observed+Model.barB*t_applied);

        % 物理系统演化
        x = x+Scenario.state_jump(:,j,b);
        [y_next,x,y_true] = cryst_plant( ...
            x,v_applied,Params,true,Scenario.noise(:,j+1,b));

        Y(:,j+1)     = y_next;
        Ytrue(:,j+1) = y_true;
        X(:,j+1)     = x;

        [History,Observation] = update_score_history( ...
            History,y_next,v_applied,j,Model,Params);

        % =====================================================================
        % 时间更新（预测下一步先验状态）：完全隔离
        % =====================================================================
        if strcmp(method, 'xue')
            % Kalman 时间更新
            x_hat_kal = Params.xue_A * x_hat_kal + Params.xue_B * t_applied;
        elseif strcmp(method, 'oh')
            % ILO 时间更新
            x_hat_oh = Params.oh_A * x_hat_oh + Params.oh_B * t_applied;
        end

        eta = Observation.u_next-u_prediction;
        r   = Observation.r_next;

        x_after = assemble_latent_state(History,Model);
        e_post  = x_after-Nominal.X_post(:,k);

        Trequest(:,k) = t_request;
        Tapplied(:,k) = Observation.t_applied;

        Eta(:,k) = eta;
        R(:,k)   = r;

        Epost(:,k)       = e_post;
        tube_post(k)     = sum((Lp\e_post).^2);
        tube_post_obs(k) = sum((Lp_obs\e_post).^2);

        eta_box(k) = all(abs(eta) <= UncSets.eta_bound*(1+1e-10));
        r_box(k)   = all(abs(r) <= UncSets.r_bound*(1+1e-10));

        eta_ellipsoid(k) = sum((Leta\eta).^2) <= 1+1e-8;
        r_ellipsoid(k)   = sum((Lr\r).^2)   <= 1+1e-8;
    end

    Batch = struct();
    Batch.method = method;
    Batch.Y      = Y;
    Batch.Y_true = Ytrue;
    Batch.V      = V;
    Batch.V_request = Vrequest;
    Batch.X      = X;
    Batch.T_request = Trequest;
    Batch.T_applied = Tapplied;
    Batch.eta = Eta;
    Batch.r   = R;
    Batch.e_pre  = Epre;
    Batch.e_post = Epost;
    Batch.eta_box_in = eta_box;
    Batch.r_box_in   = r_box;
    Batch.eta_ellipsoid_in = eta_ellipsoid;
    Batch.r_ellipsoid_in   = r_ellipsoid;
    Batch.tube_pre_quadratic      = tube_pre;
    Batch.tube_post_quadratic     = tube_post;
    Batch.tube_pre_obs_quadratic  = tube_pre_obs;
    Batch.tube_post_obs_quadratic = tube_post_obs;
    Batch.clip_flags = clip_flags;
    Batch.violation_flags = violation_flags; % <--- 新增
end
