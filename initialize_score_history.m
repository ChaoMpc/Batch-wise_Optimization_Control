function History = initialize_score_history(Y, V, Model, Params)
% =========================================================================
% initialize_score_history
%
% 在物理时间 ctrl_start 初始化得分历史。
%
% Y(:,ctrl_start+1) 是当前可用测量。
% V(:,ctrl_start) 是最后一个已施加的物理输入。
% =========================================================================

    c0 = Params.ctrl_start;
    n = Model.n_arx;
    A = Model.A_lv;

    if c0 < n
        error('initialize_score_history: insufficient pre-control history.');
    end

    History.U = zeros(A, n);
    History.T = zeros(A, n);

    for ell = 1:n
        % 当前输出时间 c0, c0-1, ...
        yc = c0 + 2 - ell;

        ydev = Y(:, yc) - Params.y_base(:, yc);
        ystd = (ydev - Model.Y_mean) ./ Model.Y_std;

        History.U(:, ell) = Model.R_u.' * ystd;

        % 先前已施加输入 c0, c0-1, ...
        vc = c0 + 1 - ell;

        vdev = V(Model.input_idx, vc) - ...
               Params.u_base(Model.input_idx, vc);

        vstd = (vdev - Model.X_mean) ./ Model.X_std;

        History.T(:, ell) = Model.R_x.' * vstd;
    end

    if any(~isfinite(History.U(:))) || any(~isfinite(History.T(:)))
        error('initialize_score_history: nonfinite score histories.');
    end
end
