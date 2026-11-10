function [Reference, Params] = generate_reference(Params)
% =========================================================================
% generate_reference
%
% Generate the baseline trajectory used by the uploaded code.
%
% Outputs:
%   Reference.V : 6 x N, baseline physical inputs
%   Reference.Y : 6 x (N+1), noise-free baseline measurements
%   Reference.X : 7 x (N+1), baseline physical states
%
% Also attaches compatibility fields to Params:
%   u_base, y_base, x_base, Ref_points,
%   u_base_rep, y_base_rep.
%
% No random numbers are used.
% =========================================================================

    N = Params.Ble;
    ni = Params.n_in;
    no = Params.n_out;

    s = (0:N-1)/max(N-1,1);

    V = zeros(ni,N);

    V(1,:) = Params.baseline.Tj_initial ...
             - Params.baseline.Tj_drop*s.^2;

    V(2,:) = Params.baseline.Fc;

    V(3,:) = Params.baseline.Fanti_initial ...
             + Params.baseline.Fanti_increase*s.^1.5;

    V(4,:) = Params.baseline.Fsol;
    V(5,:) = Params.baseline.Nstir;
    V(6,:) = Params.baseline.Fseed_after;

    n_seed = min(N, Params.baseline.Fseed_steps);
    V(6,1:n_seed) = Params.baseline.Fseed_initial;

    % A baseline outside the physical input domain is a configuration error.
    upper_violation = bsxfun(@minus, V, Params.u_phys_max);
    lower_violation = bsxfun(@minus, Params.u_phys_min, V);

    max_input_violation = max([ ...
        0;
        upper_violation(:);
        lower_violation(:) ...
    ]);

    if max_input_violation > 1e-10
        error(['generate_reference: baseline input violates physical ', ...
               'bounds by %.3e.'], max_input_violation);
    end

    X = zeros(7,N+1);
    Y = zeros(no,N+1);

    X(:,1) = Params.x0_nom;

    % Initial measurement: no state advance.
    [Y(:,1), x_initial] = cryst_plant( ...
        X(:,1), [], Params, false);

    if norm(x_initial-X(:,1),inf) ~= 0
        error('generate_reference: measurement-only mode changed state.');
    end

    for j = 1:N
        [Y(:,j+1), X(:,j+1)] = cryst_plant( ...
            X(:,j), V(:,j), Params, false);
    end

    if any(~isfinite(X(:))) || any(~isfinite(Y(:)))
        error('generate_reference: nonfinite baseline trajectory.');
    end

    %% Reference structure
    Reference = struct();

    Reference.V = V;
    Reference.Y = Y;
    Reference.X = X;

    Reference.time_y = (0:N)*Params.Ts;

    % Input application times: start of each sampling interval.
    Reference.time_v = (0:N-1)*Params.Ts;

    Reference.control_input_columns = ...
        Params.ctrl_start+1:N;

    Reference.control_output_columns = ...
        Params.ctrl_start+2:N+1;

    Reference.tracking_times = Params.Kp(:)*Params.Ts;
    Reference.tracking_targets = Params.y_ref_phys(:);

    %% Attach baseline fields used throughout subsequent stages
    Params.u_base = V;
    Params.y_base = Y;
    Params.x_base = X;

    Params.Ref_points = ...
        Params.y_ref_phys(:) ...
        - Y(Params.track_idx,Params.Kp+1)';

    Params.u_base_rep = V(Params.opt_input_idx,:)';

    % Compatibility with the uploaded code:
    % this field stores post-interval baseline outputs.
    Params.y_base_rep = Y(:,2:end)';

    %% Diagnostic results
    Reference.input_min = min(V,[],2);
    Reference.input_max = max(V,[],2);

    Reference.output_min = min(Y,[],2);
    Reference.output_max = max(Y,[],2);

    y_upper = bsxfun(@minus,Y,Params.y_phys_max);
    y_lower = bsxfun(@minus,Params.y_phys_min,Y);

    Reference.max_output_violation = max([ ...
        0;
        y_upper(:);
        y_lower(:) ...
    ]);

    fprintf('generate_reference: baseline trajectory generated.\n');
    fprintf('  Input size: %d x %d.\n',size(V,1),size(V,2));
    fprintf('  Output size: %d x %d.\n',size(Y,1),size(Y,2));

    for j = 1:numel(Params.Kp)
        fprintf(['  k=%3d: baseline C_ATR=%.8f; ', ...
                 'target=%.8f; target deviation=%+.8f.\n'], ...
            Params.Kp(j), ...
            Y(Params.track_idx,Params.Kp(j)+1), ...
            Params.y_ref_phys(j), ...
            Params.Ref_points(j));
    end
end
