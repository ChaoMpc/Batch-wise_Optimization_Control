function L_kal = design_kalman_gain(Model, Params)
    % 强行加载 control 包
    try
        pkg load control;
    catch
    end
% design_kalman_gain
% 稳态 Kalman 增益。用对偶 dlqr 比 dare 更稳健；
% 若 (barA, barC) 有不可检测模态，自动放大 Q_proc 重试。

    nx       = Model.nx;
    n_out_lv = Model.A_lv;

    Q_proc = 1e-2 * eye(nx);          % 比原来 1e-3 大，改善可检测性
    R_meas = 1e-2 * eye(n_out_lv);

    A = Model.barA;
    C = Model.barC;

    % 逐次放大 Q_proc，直到 (barA, barC) 可检测
    [detectable, ~] = check_detectability(A, C);
    if ~detectable
        for amp = [1e-1, 1, 10]
            Q_proc = amp * eye(nx);
            [detectable, ~] = check_detectability(A, C);
            if detectable, break; end
        end
    end

    if ~detectable
        warning('design_kalman_gain: (barA, barC) 不可检测，返回零增益。');
        L_kal = zeros(nx, n_out_lv);
        return;
    end


    try
        [K_dual, ~, ~] = dlqr(A', C', Q_proc, R_meas);
        L_kal = K_dual';
    catch ME
        warning('design_kalman_gain: dlqr 失败(%s)，退回零增益。', ME.message);
        L_kal = zeros(nx, n_out_lv);
    end

    if any(~isfinite(L_kal(:)))
        warning('design_kalman_gain: 非有限增益，退回零增益。');
        L_kal = zeros(nx, n_out_lv);
    end

    fprintf('design_kalman_gain: ||L||_F=%.4e，Q_proc 缩放=%.1e\n', ...
        norm(L_kal,'fro'), trace(Q_proc)/nx);
end

function [ok, n_unobs] = check_detectability(A, C)
    lam = eig(A);
    unstable = find(abs(lam) >= 1 - 1e-10);
    n_unobs = 0;
    for j = 1:numel(unstable)
        M = [lam(unstable(j))*eye(size(A,1)) - A; C];
        if rank(M) < size(A,1)
            n_unobs = n_unobs + 1;
        end
    end
    ok = (n_unobs == 0);
end
