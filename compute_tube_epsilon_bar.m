function eps_bar = compute_tube_epsilon_bar(track_idx, Tube, Model, UncSets, Params)
% 依据 Theorem 3 / Remark 9 计算单个跟踪点在 steady-state 的支撑上界。
% 假设 Q_trk 为标量对角（Q_trk = q * I），否则按完整加权矩阵处理。

    Sy   = diag(Model.Y_std);
    Qout = Model.Q;
    Pe   = Tube.P;
    Wr   = UncSets.W_r;
    Cy   = Sy * Qout * Model.barC;         % 物理输出误差对 e 的映射

    for j = 1:numel(track_idx)
        trk = track_idx(j);
        row = Cy(trk, :);

        % 两椭球 Minkowski 和的支撑
        s1 = sqrt(max(0, row * Pe * row'));                 % e in E(Pe)
        s2 = sqrt(max(0, Sy(trk,trk)^2 * Wr(trk,trk)));     % r in E(Wr)

        % Q_trk 加权
        if isscalar(Params.Q_weight)
            q = Params.Q_weight;
        else
            q = Params.Q_weight(trk, trk);
        end

        eps_bar(j) = sqrt(q) * (s1 + s2);
    end

    fprintf('[epsilon_bar] '); fprintf('%.4e  ', eps_bar); fprintf('\n');
end
