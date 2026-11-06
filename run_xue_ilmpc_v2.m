function Data = run_xue_ilmpc_v2( ...
    Problem, Initial, Optimal, Scenario, Model, Tube, UncSets, Params)
% =========================================================================
% run_xue_ilmpc_v2
% 严格遵循 Xue (2022) 的纯无约束逻辑，不作任何限幅
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

    Hessian = 2*(Lfull'*Lfull + Problem.Pprox);
    Hessian = full((Hessian+Hessian')/2);
    [Rh,flag] = chol(Hessian);
    if flag ~= 0
        error('run_xue_ilmpc: unconstrained Hessian is not positive definite.');
    end

    L_kal = design_kalman_gain(Model, Params);
    Params.xue_Kalman_L = L_kal;
    Params.xue_A = Model.barA;
    Params.xue_B = Model.barB;
    Params.xue_C = Model.barC;
    Params.xue_K_feedback = 0.1 * Tube.K;

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

    reference = Params.y_base(:,Problem.Lift.output_columns);
    ref_vec_init = reference(:);
    ref_vec_init(Problem.selected_rows) = Problem.selected_refs;
    reference = reshape(ref_vec_init, no, H);

    lambda = Params.xue_reference_lambda;

    for b = 1:B
        if b > 1
            previous_actual = Data.Y_hist(:,Problem.Lift.output_columns,b-1);
            err_latent = reference - previous_actual;
            err_vec = err_latent(:);
            err_vec(Problem.selected_rows) = 0;
            err_latent = reshape(err_vec, no, H);
            reference = reference + lambda * err_latent;
        end

        residual_previous = Wbar * ...
            (Problem.Lift.Phi_y*t_previous + ...
             Problem.Lift.y_offset - reference(:));

        f = 2*Lfull'*residual_previous;
        timer = tic();
        d = -(Rh\(Rh'\f));
        Data.QP_times(b) = toc(timer);

        % 直接叠加，不作任何限制
        t_new = t_previous + d;
        Nominal = evaluate_nominal_trajectory(Problem,t_new);

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
                 'VIOLATED channel-samples %d.\n'], ... % <--- 改为大写 VIOLATED
            b,B,Data.Err_true_act(b), ...
            Data.Err_opt_act(b),Data.Err_opt_nom(b), ...
            sum(Batch.violation_flags(:)));         % <--- 读取 violation_flags

        t_previous = t_new;
    end
    Data.final_t = t_previous;
end
