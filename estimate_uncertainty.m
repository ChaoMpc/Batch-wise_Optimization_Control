function UncSets = estimate_uncertainty(Model, Identification_data, Params)
% =========================================================================
% estimate_uncertainty
%
% 统一适配 Example 1 和 Example 2 的有限样本不确定性估计。
% 直接利用 Model 中的 U, T, Theta 序列计算 MIMO ARX 的预测误差。
% =========================================================================

    N = Params.N;
    A_lv = Params.A_lv;
    n_arx = Params.n_arx;
    n_out = Model.n_out;

    U_lv = Model.U;
    T_lv = Model.T;
    Theta = Model.Theta;

    if isempty(U_lv) || isempty(T_lv), error('estimate_uncertainty: U/T empty.'); end
    if size(U_lv,1) ~= size(T_lv,1), error('estimate_uncertainty: U/T sample mismatch.'); end

    N_data = size(U_lv,1);
    N_use  = min([N, Model.N, N_data]);

    if N_use <= n_arx, error('estimate_uncertainty: N_use <= n_arx.'); end

    % 1. 获取分位数与安全系数
    if isfield(Params,'safety_factor') && ~isempty(Params.safety_factor)
        safety_factor = Params.safety_factor;
    else
        safety_factor = 1.2;
    end
    if isfield(Params,'q_level') && ~isempty(Params.q_level)
        q_level = Params.q_level;
    else
        q_level = 0.95;
    end

    % 2. MIMO ARX innovations eta 计算
    n_reg_expect = 2 * A_lv * n_arx;
    eta_all = NaN(N_use, A_lv);
    for a = 1:A_lv
        theta = Theta{a}(:);
        if length(theta) ~= n_reg_expect
            error('estimate_uncertainty: Theta length mismatch at LV %d.', a);
        end
        for k = (n_arx+1):N_use
            u_pred = 0;
            for b = 1:A_lv
                idx_u = (b-1)*n_arx + (1:n_arx);
                idx_t = A_lv*n_arx + idx_u;
                u_lag = U_lv(k-1:-1:k-n_arx, b);
                t_lag = T_lv(k-1:-1:k-n_arx, b);
                u_pred = u_pred + theta(idx_u)'*u_lag + theta(idx_t)'*t_lag;
            end
            eta_all(k, a) = U_lv(k, a) - u_pred;
        end
    end
    eta_valid = eta_all(n_arx+1:N_use,:);
    eta_valid = eta_valid(all(isfinite(eta_valid),2),:);

    % 3. PLS 静态残差 r 计算
    R_y_model = Model.R_y(1:N_use, :);
    r_valid = R_y_model(all(isfinite(R_y_model),2),:);

    % 4. 边界与形状矩阵
    eta_q     = column_quantile(abs(eta_valid), q_level)';
    eta_bound = safety_factor * eta_q;
    eta_bound(eta_bound < 1e-8) = 1e-8;

    if A_lv == n_out
        r_bound = 1e-8 * ones(n_out,1);
        W_r = zeros(n_out);
    else
        r_q = column_quantile(abs(r_valid), q_level)';
        r_bound = safety_factor * r_q;
        r_bound(r_bound < 1e-8) = 1e-8;
        % 注意这里乘了维度包络系数
        W_r = numel(r_bound) * diag(r_bound.^2);
    end

    % 5. 状态估计误差集的形状矩阵 W_x
    R_u = Model.R_u;
    Sigma_u = R_u'*W_r*R_u;
    Sigma_u = 0.5*(Sigma_u+Sigma_u');
    n_state = 2*n_arx - 1;
    nx = A_lv*n_state;
    W_x = zeros(nx,nx);
    for a = 1:A_lv
        for b = 1:A_lv
            idx_a = (a-1)*n_state + (1:n_arx);
            idx_b = (b-1)*n_state + (1:n_arx);
            W_x(idx_a,idx_b) = Sigma_u(a,b) * eye(n_arx);
        end
    end
    W_x = 0.5*(W_x+W_x');

    % 6. 打包
    UncSets.eta_bound = eta_bound;
    UncSets.r_bound = r_bound;
    UncSets.eta = eta_all;
    UncSets.r = R_y_model;
    UncSets.W_r = W_r;
    UncSets.W_x = W_x;
    UncSets.Sigma_u = Sigma_u;
    UncSets.R_u = R_u;
    UncSets.safety_factor = safety_factor;
    UncSets.q_level = q_level;
    UncSets.n_eta_samples = size(eta_valid,1);
    UncSets.n_r_samples = size(r_valid,1);

    fprintf('estimate_uncertainty: bounds estimated.\n');
    fprintf('  eta_bound : '); fprintf('%.4e ', eta_bound); fprintf('\n');
    fprintf('  r_bound   : '); fprintf('%.4e ', r_bound);   fprintf('\n');

    if norm(W_r,'fro')>0
        fprintf('  W_r (diag): '); fprintf('%.4e ', diag(W_r)); fprintf('\n');
    else
        fprintf('  W_r       : 0 (A_lv == n_out, no truncation)\n');
    end
end

function q = column_quantile(X, p)
    Xs = sort(X,1); n = size(Xs,1);
    if n == 1, q = Xs; return; end
    idx = p*(n-1)+1; lo = floor(idx); hi = ceil(idx); w = idx-lo;
    if lo == hi, q = Xs(lo,:); else, q = (1-w)*Xs(lo,:) + w*Xs(hi,:); end
end
