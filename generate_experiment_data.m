function Data = generate_experiment_data(Params, n_batches, seed, label)
% =========================================================================
% generate_experiment_data
%
% Generate independent nonlinear plant experiments.
%
% Data.V : n_in x N x B
% Data.Y : n_out x (N+1) x B, measured
% Data.Y_true : n_out x (N+1) x B, noise-free sensor output
% Data.X : 7 x (N+1) x B
%
% X and Y_true are saved for simulation diagnostics.
% Identification uses measured Y only.
% =========================================================================

    if nargin < 4
        label = 'experiment';
    end

    if n_batches < 1 || n_batches ~= fix(n_batches)
        error('generate_experiment_data: invalid batch count.');
    end

    saved_rng = rng();
    cleanup_rng = onCleanup(@() rng(saved_rng));

    rng(seed, 'twister');

    N = Params.Ble;
    B = n_batches;

    Data = struct();

    Data.label = label;
    Data.seed = seed;
    Data.Ble = N;
    Data.N_batches = B;
    Data.ctrl_start = Params.ctrl_start;

    Data.V = zeros(Params.n_in,N,B);
    Data.Y = zeros(Params.n_out,N+1,B);
    Data.Y_true = zeros(Params.n_out,N+1,B);
    Data.X = zeros(7,N+1,B);

    fprintf('Generating %s data: %d batches.\n',label,B);

    for b = 1:B
        excitation_seed = seed + 1000*b;

        V = generate_excitation(Params,excitation_seed);

        noise = bsxfun(@times, ...
            Params.measurement_noise_halfwidth, ...
            2*rand(Params.n_out,N+1)-1);

        x = Params.x0_nom;

        Data.V(:,:,b) = V;
        Data.X(:,1,b) = x;

        [y0,y0_true] = cryst_measure(x,Params,true,noise(:,1));

        Data.Y(:,1,b) = y0;
        Data.Y_true(:,1,b) = y0_true;

        for j = 1:N
            [y,x,y_true] = cryst_plant( ...
                x,V(:,j),Params,true,noise(:,j+1));

            Data.X(:,j+1,b) = x;
            Data.Y(:,j+1,b) = y;
            Data.Y_true(:,j+1,b) = y_true;
        end

        if mod(b,25) == 0 || b == B
            fprintf('  %s: %d/%d batches completed.\n',label,b,B);
        end
    end

    Data.time_v = (0:N-1)*Params.Ts;
    Data.time_y = (0:N)*Params.Ts;

    if any(~isfinite(Data.V(:))) || ...
       any(~isfinite(Data.Y(:))) || ...
       any(~isfinite(Data.X(:)))
        error('generate_experiment_data: nonfinite generated data.');
    end
end
