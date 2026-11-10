function Validation = validate_model(Model, Standard, Data, Params)
% =========================================================================
% validate_model
%
% Free-run validation：
%   从控制段开始前测量历史初始化，递归预测整个控制段。
%
% One-step validation：
%   使用测量得到的先前输出历史进行一步预测。
%
% 补充计算：输入重构的 VAF 和 NRMSE（评价 PLS 降维对输入的保留度）
% =========================================================================

    N = Params.Ble;
    c0 = Params.ctrl_start;
    H = N - c0;

    A = Model.A_lv;
    n = Model.n_arx;
    B = Data.N_batches;

    if c0 < n
        error('validate_model: insufficient pre-control history.');
    end

    Scores = project_scores(Standard, Model);

    % --- 1. 输入重构评估 (Input Reconstruction) ---
    ni = size(Standard.X0, 2);
    SSE_in = zeros(B, ni);
    SST_in = zeros(B, ni);

    for b = 1:B
        X_b = Standard.X0(:,:,b); % N x ni
        T_b = X_b * Model.R_x;    % N x A
        X_hat = T_b * Model.P';   % N x ni

        err_in = X_b - X_hat;
        SSE_in(b,:) = sum(err_in.^2, 1);

        X_centered = bsxfun(@minus, X_b, mean(X_b,1));
        SST_in(b,:) = sum(X_centered.^2, 1);
    end

    % 防止除零
    SST_in(SST_in < 1e-12) = 1e-12;

    Input_VAF = 100 * (1 - sum(SSE_in,1) ./ sum(SST_in,1));
    Input_NRMSE = sqrt(sum(SSE_in,1) ./ sum(SST_in,1));

    % --- 2. 输出自由运行 (Free-run) 评估 ---
    Yfree = zeros(Params.n_out, H, B);
    Yone = zeros(Params.n_out, H, B);

    Eta = zeros(H*B, A);
    R = zeros(H*B, Params.n_out);

    SSE_free = zeros(B, Params.n_out);
    SSE_one = zeros(B, Params.n_out);

    NRMSE = zeros(B, Params.n_out);
    VAF = zeros(B, Params.n_out);
    VAF_one = zeros(B, Params.n_out);

    n_reg = 2*A*n;

    for b = 1:B
        Tb = Scores.T(:, :, b);
        Uobs = Scores.U(:, :, b);

        Ufree = zeros(N+1, A);

        % 控制开始前使用测量输出得分。
        Ufree(1:c0+1, :) = Uobs(1:c0+1, :);

        Uone_control = zeros(H, A);

        for j = c0+1:N
            phi_free = zeros(n_reg, 1);
            phi_one = zeros(n_reg, 1);

            for a = 1:A
                iu = (a-1)*n + (1:n);
                it = A*n + iu;

                phi_free(iu) = Ufree(j:-1:j-n+1, a);
                phi_one(iu) = Uobs(j:-1:j-n+1, a);

                phi_free(it) = Tb(j:-1:j-n+1, a);
                phi_one(it) = Tb(j:-1:j-n+1, a);
            end

            predicted_free = Model.Theta_matrix.' * phi_free;
            predicted_one = Model.Theta_matrix.' * phi_one;

            Ufree(j+1, :) = predicted_free.';

            k = j - c0;
            Uone_control(k, :) = predicted_one.';

            row = (b-1)*H + k;
            Eta(row, :) = Uobs(j+1, :) - predicted_one.';
        end

        output_columns = c0+2:N+1;

        Yfree0 = Ufree(output_columns, :) * Model.Q.';
        Yone0 = Uone_control * Model.Q.';

        Yfree_dev = bsxfun(@plus, ...
            bsxfun(@times, Yfree0, Model.Y_std.'), Model.Y_mean.');

        Yone_dev = bsxfun(@plus, ...
            bsxfun(@times, Yone0, Model.Y_std.'), Model.Y_mean.');

        Ybase = Params.y_base(:, output_columns).';

        Yfree_phys = Yfree_dev + Ybase;
        Yone_phys = Yone_dev + Ybase;

        Ytrue = Data.Y(:, output_columns, b).';

        Yfree(:, :, b) = Yfree_phys.';
        Yone(:, :, b) = Yone_phys.';

        rows = (b-1)*H + (1:H);
        R(rows, :) = Scores.R_y(output_columns, :, b);

        for out = 1:Params.n_out
            target = Ytrue(:, out);
            err_free = target - Yfree_phys(:, out);
            err_one = target - Yone_phys(:, out);

            SSE_free(b, out) = sum(err_free.^2);
            SSE_one(b, out) = sum(err_one.^2);

            centered = target - mean(target);
            denom = max(norm(centered), sqrt(H)*1e-12);

            NRMSE(b, out) = norm(err_free) / denom;

            target_variance = max(var(target, 1), 1e-24);

            VAF(b, out) = 100 * (1 - var(err_free, 1) / target_variance);
            VAF_one(b, out) = 100 * (1 - var(err_one, 1) / target_variance);
        end
    end

    Validation = struct();

    Validation.Y_free = Yfree;
    Validation.Y_one_step = Yone;

    Validation.output_columns = c0+2:N+1;
    Validation.time = (c0+1:N) * Params.Ts;

    Validation.NRMSE_per_batch = NRMSE;
    Validation.VAF_per_batch = VAF;
    Validation.VAF_one_step_per_batch = VAF_one;

    Validation.NRMSE = mean(NRMSE, 1);
    Validation.VAF = mean(VAF, 1);
    Validation.VAF_one_step = mean(VAF_one, 1);

    Validation.RMSE = sqrt(sum(SSE_free, 1) / (H*B));
    Validation.RMSE_one_step = sqrt(sum(SSE_one, 1) / (H*B));

    Validation.eta_samples = Eta;
    Validation.r_samples = R;

    % 把输入的评价指标存起来
    Validation.input_VAF = Input_VAF;
    Validation.input_NRMSE = Input_NRMSE;

    Validation.finite_prediction = ...
        all(isfinite(Yfree(:))) && all(isfinite(Yone(:)));

    if ~Validation.finite_prediction
        error('validate_model: free-run prediction became nonfinite.');
    end
end
