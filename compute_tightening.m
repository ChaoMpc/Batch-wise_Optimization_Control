function Tube = compute_tightening( ...
    Model, Feedback, Invariant, UncSets, StateError, Params)
% =========================================================================
% compute_tightening
%
% 输入方向：
%   v_actual − v_nominal = S_v*P_pls*K*(e + eps_x)
%   e in E(P_tube),  eps_x in E(W_epsilon)
%
% 输出方向：
%   y_meas − y_nominal = S_y*Q*barC*e + S_y*r
%   e in E(P_tube),  r in E(W_r)
%
% 因为 e 与 eps_x 分别有界，输入 tightening 的两项支撑函数相加。
% 输出方向不变（y 只依赖真实状态误差 e 与 PLS 残差 r）。
%
% 另外构造 P_obs = P_tube ⊕ W_epsilon 的椭球外逼近，
% 供在线诊断使用（x_obs − x_nom = e + eps_x）。
% =========================================================================

    Sv = diag(Model.X_std);
    Sy = diag(Model.Y_std);

    P_tube = Invariant.P;
    W_eps  = StateError.W_epsilon;

    Mv = Sv*Model.P*Feedback.K;
    My = Sy*Model.Q*Model.barC;

    Wv     = Mv*P_tube*Mv';
    Wv_eps = Mv*W_eps*Mv';

    Wys    = My*P_tube*My';
    Wyr    = Sy*UncSets.W_r*Sy';

    Wv     = symmetric_local(Wv);
    Wv_eps = symmetric_local(Wv_eps);
    Wys    = symmetric_local(Wys);
    Wyr    = symmetric_local(Wyr);

    Hv = Params.H_v;
    Hy = Params.H_y;

    delta_v = zeros(size(Hv,1),1);
    delta_y = zeros(size(Hy,1),1);

    for j = 1:size(Hv,1)
        h = Hv(j,:);
        state_margin = sqrt(max(0,h*Wv*h'));
        eps_margin   = sqrt(max(0,h*Wv_eps*h'));
        delta_v(j)   = state_margin + eps_margin;
    end

    for j = 1:size(Hy,1)
        h = Hy(j,:);
        state_margin          = sqrt(max(0,h*Wys*h'));
        reconstruction_margin = sqrt(max(0,h*Wyr*h'));
        delta_y(j)            = state_margin + reconstruction_margin;
    end

    % 观测状态偏差 e_obs = e + eps_x 的椭球外逼近（诊断参考集）。
    gamma = 1.0;
    P_obs = (1+gamma)*P_tube + (1+1/gamma)*W_eps;
    P_obs = symmetric_local(P_obs);

    [~,flag] = chol(P_obs);
    if flag ~= 0
        P_obs = P_obs + 1e-10*eye(size(P_obs));
    end

    Tube = struct();

    Tube.state_convention = 'observed_score_history';

    Tube.K = Feedback.K;
    Tube.A_K = Feedback.A_K;
    Tube.spectral_radius = Feedback.spectral_radius;

    Tube.P     = P_tube;
    Tube.P_obs = P_obs;
    Tube.alpha = Invariant.alpha;
    Tube.tube_beta = Invariant.tube_beta;

    Tube.W_eta = UncSets.W_eta;
    Tube.W_r   = UncSets.W_r;

    Tube.W_d      = Invariant.W_d;
    Tube.W_e_eps  = Invariant.W_e_eps;
    Tube.W_initial = StateError.W_initial;

    % 关键：W_epsilon 非零，是 tube 反馈每个时刻的额外不确定通道。
    Tube.W_epsilon = W_eps;

    Tube.W_v     = Wv;
    Tube.W_v_eps = Wv_eps;
    Tube.W_y_state = Wys;
    Tube.W_y_reconstruction = Wyr;

    Tube.H_v = Hv;
    Tube.H_y = Hy;

    Tube.h_v = Params.h_v;
    Tube.h_y = Params.h_y;

    Tube.delta_v = delta_v;
    Tube.delta_y = delta_y;

    Tube.h_v_tight = Params.h_v-delta_v;
    Tube.h_y_tight = Params.h_y-delta_y;

    Tube.invariant_ratio = Invariant.invariant_ratio;
    Tube.initial_containment_ratio = ...
        Invariant.initial_containment_ratio;

    fprintf('\nPhysical constraint tightening\n');
    fprintf('  delta_v (state + eps): [%s]\n', ...
        num2str(delta_v','%.5e '));
    fprintf('  delta_y (state + r):   [%s]\n', ...
        num2str(delta_y','%.5e '));

    ni = Params.n_in_opt;
    no = Params.n_out;

    Hv_box = [eye(ni);-eye(ni)];
    Hy_box = [eye(no);-eye(no)];

    if isequal(size(Hv),size(Hv_box)) && ...
       norm(Hv-Hv_box,'fro') < 1e-12

        vmax = Params.u_phys_max(Model.input_idx);
        vmin = Params.u_phys_min(Model.input_idx);

        Tube.v_nominal_min = vmin+delta_v(ni+1:end);
        Tube.v_nominal_max = vmax-delta_v(1:ni);

        if any(Tube.v_nominal_min > Tube.v_nominal_max)
            error('compute_tightening: tightened input box is empty.');
        end
    end

    if isequal(size(Hy),size(Hy_box)) && ...
       norm(Hy-Hy_box,'fro') < 1e-12

        Tube.y_nominal_min = ...
            Params.y_phys_min+delta_y(no+1:end);

        Tube.y_nominal_max = ...
            Params.y_phys_max-delta_y(1:no);

        if any(Tube.y_nominal_min > Tube.y_nominal_max)
            error('compute_tightening: tightened output box is empty.');
        end

        tr = Params.track_idx;

        fprintf('  Tightened C_ATR interval: [%.8f, %.8f].\n', ...
            Tube.y_nominal_min(tr),Tube.y_nominal_max(tr));
    end
end

function S = symmetric_local(M)
    S = real((M+M')/2);
end
