function V = generate_excitation(Params, seed)
% =========================================================================
% generate_excitation
%
% Generate one batch of bounded adaptive random-step excitation.
%
% V(:,j) acts over [(j-1)*Ts, j*Ts].
% Before control start, execute the baseline exactly.
%
% This function restores the incoming RNG state.
% =========================================================================

    saved_rng = rng();
    cleanup_rng = onCleanup(@() rng(saved_rng));

    rng(seed, 'twister');

    N = Params.Ble;
    V = Params.u_base;

    fraction = Params.excitation_room_fraction;
    amp_min = Params.excitation_amplitude_fraction_min;

    if fraction <= 0 || fraction > 1
        error('generate_excitation: room fraction must lie in (0,1].');
    end

    if amp_min < 0 || amp_min > 1
        error('generate_excitation: invalid amplitude fraction.');
    end

    for ii = 1:numel(Params.opt_input_idx)
        channel = Params.opt_input_idx(ii);
        period = Params.prbs_periods(channel);

        if period < 1 || period ~= fix(period)
            error('generate_excitation: periods must be positive integers.');
        end

        room_up = ...
            Params.u_phys_max(channel) - Params.u_base(channel,:);

        room_down = ...
            Params.u_base(channel,:) - Params.u_phys_min(channel);

        first = Params.ctrl_start + 1;

        while first <= N
            last = min(first + period - 1, N);

            positive = rand() >= 0.5;
            amplitude_fraction = amp_min + (1-amp_min)*rand();

            if positive
                deviation = fraction * amplitude_fraction * ...
                            room_up(first:last);
            else
                deviation = -fraction * amplitude_fraction * ...
                            room_down(first:last);
            end

            V(channel,first:last) = ...
                Params.u_base(channel,first:last) + deviation;

            first = last + 1;
        end
    end

    up = bsxfun(@minus, V, Params.u_phys_max);
    lo = bsxfun(@minus, Params.u_phys_min, V);

    violation = max([0; up(:); lo(:)]);

    if violation > 1e-10
        error('generate_excitation: input bound violation %.3e.', violation);
    end

    if norm(V(:,1:Params.ctrl_start) - ...
            Params.u_base(:,1:Params.ctrl_start), 'fro') > 1e-12
        error('generate_excitation: baseline phase was modified.');
    end
end
