function Standard = apply_preprocessing(Data, Prep, Params)
% =========================================================================
% apply_preprocessing
%
% Apply fixed TRAINING means/scales to complete batches.
%
% X0: N x n_in_opt x B
% Y0: (N+1) x n_out x B
%
% No new mean or standard deviation is estimated.
% =========================================================================

    N = Params.Ble;
    B = Data.N_batches;
    ni = numel(Prep.input_idx);

    Standard = struct();

    Standard.X0 = zeros(N,ni,B);
    Standard.Y0 = zeros(N+1,Params.n_out,B);
    Standard.N_batches = B;

    for b = 1:B
        Xdev = Data.V(Prep.input_idx,:,b) - ...
               Params.u_base(Prep.input_idx,:);

        Ydev = Data.Y(:,:,b) - Params.y_base;

        Standard.X0(:,:,b) = bsxfun(@rdivide, ...
            bsxfun(@minus,Xdev',Prep.X_mean'),Prep.X_std');

        Standard.Y0(:,:,b) = bsxfun(@rdivide, ...
            bsxfun(@minus,Ydev',Prep.Y_mean'),Prep.Y_std');
    end

    if any(~isfinite(Standard.X0(:))) || ...
       any(~isfinite(Standard.Y0(:)))
        error('apply_preprocessing: nonfinite standardized data.');
    end
end
