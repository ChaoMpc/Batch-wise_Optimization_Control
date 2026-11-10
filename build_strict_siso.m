function Model = build_strict_siso(Identification_data, Params)
% =========================================================================
% 严格拟合 SISO-ARX 模型，不包含任何跨通道耦合。
% =========================================================================
    X0 = Identification_data.X0; Y0 = Identification_data.Y0;
    N = size(X0, 1); A_lv = Params.A_lv; n_arx = Params.n_arx; n_out = Params.n_out;

    % --- 1. PLS 提取 (与 MIMO 完全一致，保证在同一潜空间) ---
    T = zeros(N, A_lv); U = zeros(N, A_lv);
    W = zeros(size(X0,2), A_lv); P = zeros(size(X0,2), A_lv);
    C = zeros(n_out, A_lv); Q = zeros(n_out, A_lv);
    X_res = X0; Y_res = Y0;
    for a = 1:A_lv
        u = Y_res(:, 1);
        for iter = 1:1000
            w = X_res' * u; w = w / norm(w);
            t = X_res * w;
            c = Y_res' * t; c = c / norm(c);
            u_new = Y_res * c;
            if u_new'*u < 0, u_new = -u_new; c = -c; end
            if norm(u_new-u)/norm(u) < 1e-6, u = u_new; break; end
            u = u_new;
        end
        p = X_res'*t/(t'*t); q = Y_res'*u/(u'*u);
        T(:,a)=t; U(:,a)=u; W(:,a)=w; P(:,a)=p; C(:,a)=c; Q(:,a)=q;
        X_res = X_res - t*p'; Y_res = Y_res - u*q';
    end
    R_y = Y0 - U * Q';
    R_x = W / (P'*W); R_u = C / (Q'*C);

    % --- 2. 严格 SISO-ARX 拟合 ---
    n_reg_mimo = 2 * A_lv * n_arx;  % 长度对齐要求，用于兼容验证函数
    n_state = 2*n_arx - 1;
    nx = A_lv * n_state;

    A_bar = zeros(nx, nx); B_bar = zeros(nx, A_lv); C_bar = zeros(A_lv, nx); E_bar = zeros(nx, A_lv);
    Theta_cell = cell(A_lv, 1);

    for a = 1:A_lv
        % 仅取本通道 (a) 的数据，不涉及其他通道
        Phi = zeros(N-n_arx, 2*n_arx);
        Y_arx = zeros(N-n_arx, 1);
        row = 0;
        for k = (n_arx+1):N
            row = row + 1;
            u_lag = U(k-1:-1:k-n_arx, a);
            t_lag = T(k-1:-1:k-n_arx, a);
            Phi(row,:) = [u_lag', t_lag'];
            Y_arx(row) = U(k, a);
        end

        lambda = Params.arx_lambda;
        reg_scale = trace(Phi'*Phi) / (2*n_arx);
        if reg_scale < 1e-12, reg_scale = 1; end
        theta_siso = (Phi'*Phi + lambda * reg_scale * eye(2*n_arx)) \ (Phi' * Y_arx);

        % 构造稀疏 Theta 兼容原验证函数（跨通道系数严格置 0）
        theta_mimo = zeros(n_reg_mimo, 1);
        theta_mimo((a-1)*n_arx + (1:n_arx)) = theta_siso(1:n_arx);
        theta_mimo(A_lv*n_arx + (a-1)*n_arx + (1:n_arx)) = theta_siso(n_arx+1:end);
        Theta_cell{a} = theta_mimo;

        % 状态空间实现 (仅对角块，无跨通道)
        idx_a = (a-1)*n_state + (1:n_state);
        alpha = theta_siso(1:n_arx);
        beta = theta_siso(n_arx+1:end);

        A_bar(idx_a(1), idx_a(1:n_arx)) = alpha';
        if n_arx > 1
            A_bar(idx_a(1), idx_a(n_arx+1:2*n_arx-1)) = beta(2:end)';
        end
        B_bar(idx_a(1), a) = beta(1);
        C_bar(a, idx_a(1)) = 1;
        E_bar(idx_a(1), a) = 1;

        for ell = 2:n_arx
            A_bar(idx_a(ell), idx_a(ell-1)) = 1;
        end
        if n_arx > 1
            B_bar(idx_a(n_arx+1), a) = 1;
            for ell = n_arx+2:2*n_arx-1
                A_bar(idx_a(ell), idx_a(ell-1)) = 1;
            end
        end
    end

    % 打包
    Model = struct();
    Model.barA = A_bar; Model.barB = B_bar; Model.barC = C_bar; Model.barE = E_bar;
    Model.T = T; Model.U = U; Model.Q = Q; Model.R_u = R_u; Model.R_y = R_y; Model.R_x = R_x;
    Model.P = P; Model.W = W; Model.C = C;
    Model.A_lv = A_lv; Model.n_arx = n_arx; Model.n_out = n_out; Model.n_state = n_state;
    Model.Theta = Theta_cell;
    Model.X_std = Identification_data.X_std; Model.Y_std = Identification_data.Y_std;
    Model.X_mean_dev = Identification_data.X_mean_dev; Model.Y_mean_dev = Identification_data.Y_mean_dev;
    Model.y_ss = Identification_data.y_ss; Model.u_ss = Identification_data.u_ss;
end
