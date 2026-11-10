function Invariant = compute_invariant_ellipsoid( ...
    Model, Feedback, UncSets, StateError, Params)
% =========================================================================
% compute_invariant_ellipsoid
% 包含 RPI 椭球计算，以及 alpha / beta 的 2D 敏感度网格扫描
% =========================================================================

    AK = Feedback.A_K;
    nx = size(AK,1);

    % 1. 预计算基础扰动矩阵
    We_eta = Model.barE * UncSets.W_eta * Model.barE';
    BK = Model.barB * Feedback.K;
    We_eps = BK * StateError.W_epsilon * BK';

    We_eta = symmetric_local(We_eta);
    We_eps = symmetric_local(We_eps);

    W0 = symmetric_local(StateError.W_initial);

    rho = max(abs(eig(AK)));
    rho2 = rho^2;
    alpha_limit = 1/rho2 - 1;

    % 选定用于实际控制的 beta
    beta = Params.tube_beta;
    Wd = (1+beta)*We_eps + (1+1/beta)*We_eta;
    Wd = symmetric_local(Wd);

    % =====================================================================
    % 敏感度分析扫描 (Sensitivity Scan)
    % 扫描不同的 alpha 和 beta，记录 W_e 和 P_e 的 trace
    % =====================================================================
    % 1. 基于稳定性条件过滤 alpha 候选值
    alpha_candidates = Params.tube_alpha_grid(:)';
    alpha_candidates = alpha_candidates(alpha_candidates > 0 & (1+alpha_candidates)*rho2 < 1-1e-10);

    % 2. 从 init_params 读取 beta 的扫描网格
    beta_candidates = Params.tube_beta_grid(:)';

    num_a = numel(alpha_candidates);
    num_b = numel(beta_candidates);

    sens_mat_traceWe = zeros(num_a, num_b);
    sens_mat_tracePe = zeros(num_a, num_b);

    for ia = 1:num_a
        alpha = alpha_candidates(ia);
        As = sqrt(1+alpha)*AK;
        for ib = 1:num_b
            b_test = beta_candidates(ib);
            We_test = (1+b_test)*We_eps + (1+1/b_test)*We_eta;
            sens_mat_traceWe(ia, ib) = trace(We_test);

            % 求解此时的 P 矩阵
            try
                Qshape = (1+1/alpha)*We_test + 1e-12*eye(nx);
                P_test = dlyap(As, Qshape);
                if any(~isfinite(P_test(:)))
                    sens_mat_tracePe(ia, ib) = NaN;
                else
                    sens_mat_tracePe(ia, ib) = trace(P_test);
                end
            catch
                sens_mat_tracePe(ia, ib) = NaN;
            end
        end
    end

    Sens = struct();
    Sens.alpha_candidates = alpha_candidates;
    Sens.beta_candidates = beta_candidates;
    Sens.trace_We = sens_mat_traceWe;
    Sens.trace_Pe = sens_mat_tracePe;
    % =====================================================================

    % --- 正常的主循环（寻找最优 alpha） ---
    candidates = unique([Params.tube_alpha_grid(:)', alpha_limit*[0.02,0.1,0.3,0.6]]);
    candidates = candidates(isfinite(candidates) & candidates > 0 & (1+candidates)*rho2 < 1-1e-10);

    if isempty(candidates), error('No valid alpha candidate.'); end

    relative_regularization = Params.tube_shape_regularization;
    numerical_scale = max([trace(Wd)/nx, trace(W0)/nx, realmin]);
    regularization = max(relative_regularization*numerical_scale, 1e-30);

    best = [];
    best_score = Inf;

    for j = 1:numel(candidates)
        alpha = candidates(j);
        As = sqrt(1+alpha)*AK;
        Qshape = (1+1/alpha)*Wd + regularization*eye(nx);

        try
            P = dlyap(As,Qshape);
            P = symmetric_local(P);
            [L,flag] = chol(P,'lower');
            if flag ~= 0, continue; end

            image_shape = (1+alpha)*AK*P*AK' + (1+1/alpha)*Wd;
            invariant_ratio = max(eig(symmetric_local(L\image_shape/L')));
            initial_ratio = max(eig(symmetric_local(L\W0/L')));
            lyapunov_residual = P-As*P*As'-Qshape;
            relative_residual = norm(lyapunov_residual,'fro') / max(norm(P,'fro'),realmin);

            if any(~isfinite([invariant_ratio,initial_ratio,relative_residual])) || ...
               invariant_ratio > 1+1e-7 || initial_ratio > 1+1e-7 || relative_residual > 1e-8
                continue;
            end

            score = trace(P);
            if score < best_score
                best_score = score;
                best = struct();
                best.P = P;
                best.alpha = alpha;
                best.initial_scale = 1;
                best.invariant_ratio = invariant_ratio;
                best.initial_containment_ratio = initial_ratio;
                best.lyapunov_relative_residual = relative_residual;
                best.rcond_P = rcond(P);
            end
        catch
        end
    end

    if isempty(best), error('compute_invariant_ellipsoid: no candidate passed.'); end

    Invariant = best;
    Invariant.W_d = Wd;
    Invariant.W_e_eta = We_eta;
    Invariant.W_e_eps = We_eps;
    Invariant.tube_beta = beta;
    Invariant.shape_regularization = regularization;

    % 将敏感度分析结果存入返回值
    Invariant.Sens = Sens;

    fprintf('\nRPI ellipsoid construction (Decoupled W0)\n');
    fprintf('  Selected alpha: %.4f | Selected beta: %.4f\n', Invariant.alpha, Invariant.tube_beta);
    fprintf('  trace(W_e): %.4e | trace(P_e): %.4e\n', trace(Wd), trace(Invariant.P));
end

function S = symmetric_local(M)
    S = real((M+M')/2);
end
