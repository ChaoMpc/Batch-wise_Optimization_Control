function ARX = fit_latent_arx(Scores, n_arx, Params)
% =========================================================================
% fit_latent_arx
%
% 统一 MIMO ARX 建模：在 PLS 提取完成并得到全部潜变量得分后进行。
%
% 对第 a 个输出潜变量：
%   u_a(j+1) 由全部 A 个潜变量的历史
%       u_b(j),...,u_b(j-n+1)
%       t_b(j),...,t_b(j-n+1)
%   联合解释。
%
% Theta 排列：
%   [U 通道 1 历史;
%    ...
%    U 通道 A 历史;
%    T 通道 1 历史;
%    ...
%    T 通道 A 历史]
%
% Theta_matrix 的第 a 列对应输出潜变量 u_a 的 MIMO ARX 参数。
% =========================================================================

    n = n_arx;
    A = size(Scores.T, 2);
    B = Scores.N_batches;
    N = Params.Ble;

    if size(Scores.U, 2) ~= A
        error('fit_latent_arx: Scores.U and Scores.T have inconsistent latent dimensions.');
    end

    if n < 1 || n ~= fix(n) || n > Params.ctrl_start
        error('fit_latent_arx: invalid ARX order.');
    end

    % 需要前 n 个历史样本，因此从 ctrl_start+n 开始。
    first_input_column = Params.ctrl_start + n;

    if first_input_column > N
        error('fit_latent_arx: insufficient control-segment samples.');
    end

    count_per_batch = N - first_input_column + 1;
    ns = B * count_per_batch;
    nr = 2 * A * n;

    Phi = zeros(ns, nr);
    Target = zeros(ns, A);

    row = 0;

    for b = 1:B
        Tb = Scores.T(:, :, b);   % N   x A
        Ub = Scores.U(:, :, b);   % N+1 x A

        for j = first_input_column:N
            row = row + 1;

            phi = zeros(nr, 1);

            % 统一构造全部 A 个潜变量的回归向量。
            for a = 1:A
                iu = (a - 1) * n + (1:n);
                it = A * n + iu;

                phi(iu) = Ub(j:-1:j-n+1, a);
                phi(it) = Tb(j:-1:j-n+1, a);
            end

            Phi(row, :) = phi.';
            Target(row, :) = Ub(j+1, :);
        end
    end

    lambda = Params.arx_lambda;

    if ~isscalar(lambda) || ~isfinite(lambda) || lambda < 0
        error('fit_latent_arx: invalid ridge parameter.');
    end

    if lambda > 0
        Theta_matrix = ...
            [Phi; sqrt(lambda) * eye(nr)] \ ...
            [Target; zeros(nr, A)];
    else
        if rank(Phi) < nr
            error('fit_latent_arx: unregularized regression is rank deficient.');
        end

        Theta_matrix = Phi \ Target;
    end

    residual = Target - Phi * Theta_matrix;

    ARX = struct();

    ARX.A_lv = A;
    ARX.n_arx = n;
    ARX.lambda = lambda;

    ARX.Theta_matrix = Theta_matrix;
    ARX.Theta = cell(A, 1);

    for a = 1:A
        ARX.Theta{a} = Theta_matrix(:, a);
    end

    ARX.training_innovation = residual;
    ARX.n_samples = ns;
    ARX.n_parameters = nr * A;

    ARX.regression_condition = cond(Phi);
    ARX.training_rmse = sqrt(mean(residual.^2, 1));

    % 删除了包含 @(a) 函数句柄的 regressor_layout 字段，避免 save 报错。
    % 回归向量的排列逻辑已在上方注释中完整说明，并由后续 realizer 和
    % residual 采集函数严格遵守。
end
