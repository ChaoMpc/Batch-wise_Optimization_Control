function Prep = fit_preprocessing(Data, Params)
% =========================================================================
% fit_preprocessing
%
% Fit training-only centering/scaling on the control segment.
%
% Same-time pairing:
%   input V(:,j), measured output Y(:,j)
% correspond to input application time j-1.
%
% Dynamic ARX subsequently predicts Y(:,j+1).
% =========================================================================

    columns = Params.ctrl_start+1:Params.Ble;
    input_idx = Params.opt_input_idx(:);

    X = zeros(numel(columns)*Data.N_batches,numel(input_idx));
    Y = zeros(numel(columns)*Data.N_batches,Params.n_out);

    row = 0;

    for b = 1:Data.N_batches
        Vdev = Data.V(input_idx,columns,b) - ...
               Params.u_base(input_idx,columns);

        Ydev = Data.Y(:,columns,b) - ...
               Params.y_base(:,columns);

        rows = row + (1:numel(columns));

        X(rows,:) = Vdev';
        Y(rows,:) = Ydev';

        row = row + numel(columns);
    end

    Prep = struct();

    Prep.X_mean = mean(X,1)';
    Prep.Y_mean = mean(Y,1)';

    Prep.X_std = std(X,0,1)';
    Prep.Y_std = std(Y,0,1)';

    % Preserve the current all-one configuration.
    % Nontrivial weighting is not silently mixed into physical scaling.
    if isfield(Params,'W_y') && ...
       any(abs(Params.W_y(:)-1) > 1e-12)
        error(['fit_preprocessing: this version requires W_y=ones. ', ...
               'Use explicit PLS weighting if needed.']);
    end

    if any(Prep.X_std < 1e-12)
        error('fit_preprocessing: a retained physical input is unexcited.');
    end

    if any(Prep.Y_std < 1e-12)
        error('fit_preprocessing: an output has negligible training variance.');
    end

    Prep.input_idx = input_idx;
    Prep.training_columns = columns;
    Prep.n_samples = size(X,1);

    Prep.X0 = bsxfun(@rdivide, ...
        bsxfun(@minus,X,Prep.X_mean'),Prep.X_std');

    Prep.Y0 = bsxfun(@rdivide, ...
        bsxfun(@minus,Y,Prep.Y_mean'),Prep.Y_std');

    fprintf('fit_preprocessing: %d control-segment samples.\n', ...
        Prep.n_samples);
end
