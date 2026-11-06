function Scenario = generate_simulation_scenario(Params)
% =========================================================================
% generate_simulation_scenario
%
% Shared physical initial conditions, measurement noise, and disturbances.
%
% noise(:,k+1,b): additive physical measurement noise at physical time k.
% state_jump(:,j,b): added BEFORE applying input column j.
% =========================================================================

    saved_rng = rng();
    cleanup_rng = onCleanup(@() rng(saved_rng));

    rng(Params.rng_seed_online,'twister');

    N = Params.Ble;
    B = Params.N_batches;

    Scenario = struct();

    Scenario.seed = Params.rng_seed_online;
    Scenario.X0 = zeros(7,B);

    Scenario.noise = zeros(Params.n_out,N+1,B);
    Scenario.state_jump = zeros(7,N,B);

    perturbation = Params.init_perturb_halfwidth(:);

    if ~Params.enable_init_perturb
        perturbation = zeros(7,1);
    end

    for b = 1:B
        Scenario.X0(:,b) = ...
            Params.x0_nom + perturbation.*(2*rand(7,1)-1);

        Scenario.noise(:,:,b) = bsxfun(@times, ...
            Params.measurement_noise_halfwidth, ...
            2*rand(Params.n_out,N+1)-1);
    end

    if Params.disturbance.enabled
        batches = Params.disturbance.batches;
        columns = Params.disturbance.input_columns;
        state_idx = Params.disturbance.state_idx;

        if any(batches < 1) || any(batches ~= fix(batches)) || ...
           any(columns < 1) || any(columns > N) || ...
           any(columns ~= fix(columns)) || ...
           state_idx < 1 || state_idx > 7 || ...
           state_idx ~= fix(state_idx)
            error('generate_simulation_scenario: invalid disturbance indices.');
        end

        batches = batches(batches <= B);

        for b = batches(:)'
            Scenario.state_jump(state_idx,columns,b) = ...
                Params.disturbance.increment;
        end
    end

    Scenario.description = ...
        'Identical initial conditions, noise, and state jumps for both methods';
end
