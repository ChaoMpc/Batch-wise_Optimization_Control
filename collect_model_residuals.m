function Residuals = collect_model_residuals(Data, Model, Prep, Params)
% =========================================================================
% collect_model_residuals
%
% 状态约定：
%   observed output scores + known applied input score histories.
%
% eta：
%   next observed output score - MIMO ARX one-step prediction.
%
% r：
%   standardized measured output - Q*observed output score.
%
% 回归向量与 fit_latent_arx.m 完全一致：包含全部 A 个潜变量通道。
% =========================================================================

    Standard = apply_preprocessing(Data, Prep, Params);
    Scores = project_scores(Standard, Model);

    N = Params.Ble;
    c0 = Params.ctrl_start;
    H = N - c0;

    A = Model.A_lv;
    n = Model.n_arx;
    block = Model.n_state_per_lv;
    B = Data.N_batches;

    if c0 < n
        error('collect_model_residuals: insufficient initial history.');
    end

    Eta = zeros(H * B, A);
    R = zeros(H * B, Params.n_out);
    Xinitial = zeros(Model.nx, B);

    realization_error = 0;

    for b = 1:B
        T = Scores.T(:, :, b);   % N   x A
        U = Scores.U(:, :, b);   % N+1 x A

        for j = c0+1:N
            x = zeros(Model.nx, 1);
            phi = zeros(2*A*n, 1);

            for a = 1:A
                idx = (a-1)*block + (1:block);
                iu = (a-1)*n + (1:n);
                it = A*n + iu;

                % 状态 x_a(j)
                x(idx(1:n)) = U(j:-1:j-n+1, a);

                if n > 1
                    x(idx(n+1:end)) = T(j-1:-1:j-n+1, a);
                end

                % 统一 MIMO ARX 回归向量
                phi(iu) = U(j:-1:j-n+1, a);
                phi(it) = T(j:-1:j-n+1, a);
            end

            if j == c0+1
                Xinitial(:, b) = x;
            end

            prediction_arx = Model.Theta_matrix.' * phi;

            prediction_ss = Model.barC * ...
                (Model.barA * x + Model.barB * T(j, :).');

            realization_error = max(realization_error, ...
                norm(prediction_arx - prediction_ss, inf));

            row = (b-1)*H + (j - c0);

            Eta(row, :) = U(j+1, :) - prediction_arx.';
            R(row, :) = Scores.R_y(j+1, :, b);
        end
    end

    if any(~isfinite(Eta(:))) || any(~isfinite(R(:)))
        error('collect_model_residuals: nonfinite residuals.');
    end

    if realization_error > 1e-8
        error('collect_model_residuals: ARX/state-space realization mismatch %.3e.', ...
            realization_error);
    end

    Residuals = struct();

    Residuals.label = Data.label;
    Residuals.eta = Eta;
    Residuals.r = R;
    Residuals.x_initial = Xinitial;
    Residuals.e_initial = bsxfun(@minus, Xinitial, Model.x0_ctrl);

    Residuals.n_samples = size(Eta, 1);
    Residuals.n_batches = B;
    Residuals.realization_error = realization_error;
end
