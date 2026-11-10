function [Data_out, Params] = gen_validation_data(Params, data_type, TrainStats)
% =========================================================================
% gen_validation_data -- Generate validation or test data.
%
% If TrainStats is provided (as Identification_data from the training run),
% X0/Y0 are recomputed using the TRAINING statistics, in accordance with
% the paper statement:
%   "all subsequent training, validation, and online operations use
%    the training-set parameters bar_v, bar_y, S_v, S_y."
%
% The global RNG state is saved and restored so that subsequent closed-loop
% simulations are reproducible.
% =========================================================================

    if nargin < 2, data_type = 'val';  end
    if nargin < 3, TrainStats = [];    end

    % --- Save RNG state ---
    if exist('rng', 'builtin')
        rng_state_before = rng;
        has_rng = true;
    else
        rand_state_before  = rand('state');
        randn_state_before = randn('state');
        has_rng = false;
    end

    saved_seed = Params.rng_seed;
    saved_N    = Params.N;

    if strcmp(data_type, 'val')
        Params.rng_seed = Params.rng_seed_val;
        Params.N        = Params.N_val;
    elseif strcmp(data_type, 'test')
        Params.rng_seed = Params.rng_seed_test;
        Params.N        = Params.N_test;
    else
        error('gen_validation_data: Unknown data_type "%s".', data_type);
    end

    [Data_out, Params] = gen_prbs_excitation(Params);

    Params.rng_seed = saved_seed;
    Params.N        = saved_N;

    % --- Re-standardize using TRAINING statistics if available ---
    if ~isempty(TrainStats)
        Nv = size(Data_out.X, 1);

        if isfield(TrainStats, 'X_mean_dev')
            xm = TrainStats.X_mean_dev;
        else
            xm = TrainStats.X_mean;
        end
        if isfield(TrainStats, 'Y_mean_dev')
            ym = TrainStats.Y_mean_dev;
        else
            ym = TrainStats.Y_mean;
        end

        Data_out.X0 = (Data_out.X - ones(Nv,1) * xm') ./ (ones(Nv,1) * TrainStats.X_std');
        Data_out.Y0 = (Data_out.Y - ones(Nv,1) * ym') ./ (ones(Nv,1) * TrainStats.Y_std');

        % Store the training statistics inside the data struct for traceability
        Data_out.TrainStats = TrainStats;
    end

    % --- Restore RNG state ---
    if has_rng
        rng(rng_state_before);
    else
        rand('state',  rand_state_before);
        randn('state', randn_state_before);
    end
end
