function StateError = compute_state_error_set(Model, Params)
    if ~isfield(Params,'init_perturb_halfwidth')
        error('compute_state_error_set: initial perturbation bounds missing.');
    end
    if any(Params.init_perturb_halfwidth(:) ~= 0)
        error('compute_state_error_set: nonzero physical initial perturbations.');
    end

    n = Model.n_arx;
    A = Model.A_lv;
    no = Params.n_out;
    block = Model.n_state_per_lv;
    nx = Model.nx;
    halfwidth = Params.measurement_noise_halfwidth(:);

    Gscore = Model.R_u' * diag(1./Model.Y_std);
    dimension = n*no;
    Ginitial = zeros(nx,dimension);

    for a = 1:A
        idx = (a-1)*block+(1:block);
        for ell = 1:n
            cols = (ell-1)*no+(1:no);
            Ginitial(idx(ell),cols) = Gscore(a,:);
        end
    end

    noise_bounds = repmat(halfwidth,n,1);

    % ====== 退回：乘以维度 dimension ======
    primitive_shape = dimension * diag(noise_bounds.^2);

    W_initial = Ginitial*primitive_shape*Ginitial';
    W_initial = (W_initial+W_initial')/2;

    StateError = struct();
    StateError.state_convention = 'observed_score_history';
    StateError.W_initial = W_initial;
    StateError.G_initial_noise = Ginitial;
    StateError.primitive_noise_bound = noise_bounds;

    % ====== 退回：在线误差等于初始误差 ======
    StateError.W_epsilon = W_initial;

    fprintf('\nInitial state-error set (Original Box-to-Ellipsoid)\n');
    fprintf('  State dimension: %d.\n',nx);
    fprintf('  Primitive noise dimension: %d.\n',dimension);
    fprintf('  Initial shape trace: %.6e.\n',trace(W_initial));
    fprintf('  Online reconstruction error shape trace: %.6e.\n',trace(StateError.W_epsilon));
end
