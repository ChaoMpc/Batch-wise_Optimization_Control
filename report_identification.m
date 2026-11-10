function report_identification(Model, Prep, Validation, Params)
% =========================================================================
% report_identification
%
% Report model dimensions, PLS reconstruction evidence, and validation.
% =========================================================================

    A = Model.A_lv;
    H = Params.Ble-Params.ctrl_start;

    ni = numel(Params.opt_input_idx);
    no = Params.n_out;

    Xcov = Prep.X0'*Prep.X0/size(Prep.X0,1);

    T = Model.PLS.T;
    Tcov = T'*T/size(T,1);

    fprintf('\n========================================================\n');
    fprintf('Identification report\n');
    fprintf('========================================================\n');

    fprintf('Physical optimized input dimension: %d\n',ni);
    fprintf('Physical measured output dimension: %d\n',no);
    fprintf('Retained latent dimension: %d\n',A);

    fprintf('ARX order: %d\n',Model.n_arx);
    fprintf('State dimension: %d\n',Model.nx);

    fprintf('MIMO ARX coefficient count: %d\n', ...
        Model.ARX.n_parameters);

    fprintf('Open-loop spectral radius: %.8f\n', ...
        Model.open_loop_spectral_radius);

    fprintf('\nPLS cumulative reconstruction evidence:\n');

    for a = 1:A
        fprintf('  LV %d: EV_X=%8.3f%%, EV_Y=%8.3f%%\n', ...
            a,100*Model.PLS.EV_X(a),100*Model.PLS.EV_Y(a));
    end

    fprintf('  sigma_min(Q): %.6e\n',Model.PLS.sigma_min_Q);
    fprintf('  ||R_x''*P-I||_F: %.6e\n', ...
        Model.PLS.input_inverse_error);
    fprintf('  ||R_u''*Q-I||_F: %.6e\n', ...
        Model.PLS.output_inverse_error);

    fprintf('\nCovariance conditioning:\n');
    fprintf('  cond(X0''*X0/M): %.6e\n',cond(Xcov));
    fprintf('  cond(T''*T/M):   %.6e\n',cond(Tcov));

    fprintf('\nActual optimization horizon: %d\n',H);
    fprintf('  Direct input decision variables: %d\n',ni*H);
    fprintf('  Latent input decision variables: %d\n',A*H);

    fprintf('\nIndependent selection-set validation:\n');
    fprintf(' Channel       NRMSE      VAF(%%)     One-step VAF(%%)\n');

    for j = 1:no
        fprintf(' %-8s %11.5f %11.3f %19.3f\n', ...
            Params.output_names{j}, ...
            Validation.NRMSE(j), ...
            Validation.VAF(j), ...
            Validation.VAF_one_step(j));
    end

    fprintf('========================================================\n');
end
