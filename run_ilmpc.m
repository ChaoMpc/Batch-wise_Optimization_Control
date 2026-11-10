function Data = run_ilmpc( ...
    Problem, Initial, Optimal, Scenario, Model, Tube, UncSets, Params)
% =========================================================================
% run_ilmpc
%
% Fixed-model nominal proximal batch updates + within-batch Tube feedback.
% =========================================================================

    B = Params.N_batches;
    N = Params.Ble;
    H = Problem.H;

    Data = struct();

    Data.method = 'Proposed Tube';

    Data.Y_hist = zeros(Params.n_out,N+1,B);
    Data.Y_true_hist = zeros(Params.n_out,N+1,B);

    Data.V_hist = zeros(Params.n_in,N,B);
    Data.V_request_hist = zeros(Params.n_in,N,B);

    Data.X_hist = zeros(7,N+1,B);

    Data.Y_nominal_hist = zeros(Params.n_out,N+1,B);
    Data.V_nominal_hist = zeros(Params.n_in,N,B);

    Data.T_hist = zeros(Model.A_lv,H,B);
    Data.QP_times = zeros(B,1);

    Data.nominal_cost = zeros(B,1);
    Data.step_norm = zeros(B,1);

    Data.Err_true_act = zeros(B,1);
    Data.Err_opt_act = zeros(B,1);
    Data.Err_opt_nom = zeros(B,1);

    Data.Batch = cell(B,1);

    t_previous = Initial.t;

    for b = 1:B
        Update = solve_batch_update(Problem,t_previous,Params);
        Nominal = Update.Trajectory;

        Batch = simulate_control_batch( ...
            Nominal,'tube',Scenario,b,Model,Tube,UncSets,Params);

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

        Data.T_hist(:,:,b) = Nominal.T;
        Data.QP_times(b) = Update.QP.elapsed;

        Data.nominal_cost(b) = Update.new_cost;
        Data.step_norm(b) = Update.step_norm;

        Ycontrol = Batch.Y(:,Problem.Lift.output_columns);
        yvector = Ycontrol(:);

        actual_selected = yvector(Problem.selected_rows);
        nominal_selected = Nominal.selected_outputs;

        Data.Err_true_act(b) = ...
            norm(actual_selected-Problem.selected_refs);

        Data.Err_opt_act(b) = ...
            norm(actual_selected-Optimal.selected_outputs);

        Data.Err_opt_nom(b) = ...
            norm(nominal_selected-Optimal.selected_outputs);

        Data.Batch{b} = Batch;

        fprintf(['Tube batch %2d/%2d: act/ref %.4e; ', ...
                 'act/opt %.4e; nom/opt %.4e; ', ...
                 'eta coverage %.1f%%; Tube coverage %.1f%%.\n'], ...
            b,B,Data.Err_true_act(b), ...
            Data.Err_opt_act(b),Data.Err_opt_nom(b), ...
            100*mean(Batch.eta_ellipsoid_in), ...
            100*mean(Batch.tube_post_quadratic <= 1+1e-7));

        t_previous = Update.t;
    end

    Data.final_t = t_previous;
end
