function [y_meas, y_true] = cryst_measure(x, Params, is_noisy, noise_sample)
% =========================================================================
% cryst_measure
%
% Measurement model for the crystallization plant.
%
% Inputs:
%   x            : physical state, 7 x 1
%   Params       : parameter structure
%   is_noisy     : true/false; default true
%   noise_sample : optional additive noise in PHYSICAL units, 6 x 1
%
% If noise_sample is supplied, it is used directly.
% Otherwise bounded uniform noise is generated when is_noisy=true.
%
% This function never advances the plant state.
% =========================================================================

    if nargin < 3 || isempty(is_noisy)
        is_noisy = true;
    end

    if nargin < 4
        noise_sample = [];
    end

    x = x(:);

    if numel(x) ~= 7 || any(~isfinite(x))
        error('cryst_measure: x must be a finite 7-element state.');
    end

    T = x(1);
    C = x(3);

    y_true = zeros(6, 1);

    y_true(1) = T;
    y_true(2) = T + 0.08;
    y_true(3) = T - 0.06;

    y_true(4) = max(0, C);
    y_true(5) = max(0, C + 2.0e-4 + 1.5e-4*(T-300));
    y_true(6) = max(0, C - 1.5e-4 - 1.0e-4*(T-300));

    if is_noisy
        if isempty(noise_sample)
            halfwidth = Params.measurement_noise_halfwidth(:);

            if numel(halfwidth) ~= 6 || ...
               any(~isfinite(halfwidth)) || any(halfwidth < 0)
                error('cryst_measure: invalid measurement noise bounds.');
            end

            noise_sample = halfwidth .* (2*rand(6,1)-1);
        else
            noise_sample = noise_sample(:);

            if numel(noise_sample) ~= 6 || ...
               any(~isfinite(noise_sample))
                error('cryst_measure: invalid explicit noise sample.');
            end
        end

        y_meas = max(0, y_true + noise_sample);
    else
        y_meas = y_true;
    end
end
