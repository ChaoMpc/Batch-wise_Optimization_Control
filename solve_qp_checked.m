function [x,Report] = solve_qp_checked( ...
    H,f,A,b,x0,Params,label)
% =========================================================================
% solve_qp_checked
%
% Solve:
%   min 0.5*x'*H*x + f'*x
%   s.t. A*x <= b
%
% Uses GNU Octave qp.
% Requires a feasible initial point.
% =========================================================================

    if nargin < 7
        label = 'QP';
    end

    f = f(:);
    b = b(:);
    x0 = x0(:);

    n = numel(f);

    if ~isequal(size(H),[n,n]) || ...
       size(A,2) ~= n || size(A,1) ~= numel(b) || ...
       numel(x0) ~= n
        error('solve_qp_checked: invalid dimensions for %s.',label);
    end

    inputs = {H,f,A,b,x0};

    for j = 1:numel(inputs)
        if any(~isfinite(inputs{j}(:)))
            error('solve_qp_checked: nonfinite data in %s.',label);
        end
    end

    if norm(H-H','fro') > 1e-9*max(1,norm(H,'fro'))
        error('solve_qp_checked: nonsymmetric Hessian in %s.',label);
    end

    H = full((H+H')/2);

    eigenvalues = eig(H);
    hscale = max(norm(H,2),realmin);

    if min(eigenvalues) < -1e-10*hscale
        error('solve_qp_checked: nonconvex Hessian in %s.',label);
    end

    % 为了防止多约束下的求解器卡死，适当加入一点点正则化
    min_eig = min(eigenvalues);
    if min_eig < 1e-6 * hscale
        H = H + (1e-6 * hscale - min_eig) * eye(n);
    end

    %% Constraint scaling
    row_scale = max(abs(A),[],2);
    zero_rows = row_scale == 0;

    if any(b(zero_rows) < 0)
        error('solve_qp_checked: impossible constant constraint in %s.',label);
    end

    row_scale(zero_rows) = 1;

    As = bsxfun(@rdivide,A,row_scale);
    bs = b./row_scale;

    As = As(~zero_rows,:);
    bs = bs(~zero_rows);

    initial_violation = max([0;As*x0-bs]);

    if initial_violation > Params.qp_feas_tol
        error(['solve_qp_checked: infeasible warm start in %s; ', ...
               'scaled violation %.3e.'],label,initial_violation);
    end

    %% Scale the complete objective by one positive scalar
    objective_scale = max([ ...
        norm(H,inf),norm(f,inf),realmin ...
    ]);

    Hs = H/objective_scale;
    fs = f/objective_scale;

    options.MaxIter = Params.qp_max_iter;
    options.TolX = Params.qp_tol_x;

    timer = tic();

    try
        [x,~,info,lambda] = qp( ...
            x0,Hs,fs,[],[], ...
            -Inf(n,1),Inf(n,1), ...
            -Inf(size(bs)),full(As),bs,options);
    catch ME
        error('solve_qp_checked: %s raised an error: %s',label,ME.message);
    end

    elapsed = toc(timer);

    % ====== 新增：针对 status=3 (最大迭代次数) 的宽容处理 ======
    if info.info ~= 0
        if info.info == 3 && ~isempty(x) && all(isfinite(x))
            % 检查它返回的这个解是否满足约束
            check_violation = max([0;As*x-bs]);
            if check_violation <= Params.qp_feas_tol
                fprintf('  [Warning] %s reached max iterations, but a feasible solution was accepted.\n', label);
            else
                error(['solve_qp_checked: %s reached max iterations AND is infeasible (violation %.3e).'], ...
                      label, check_violation);
            end
        else
            error(['solve_qp_checked: %s did not return a convex optimum; ', ...
                   'qp status=%d, iterations=%d.'], ...
                  label,info.info,info.solveiter);
        end
    end

    if isempty(x) || any(~isfinite(x))
        error('solve_qp_checked: invalid solution returned for %s.',label);
    end

    scaled_violation = max([0;As*x-bs]);
    raw_violation = max([0;A*x-b]);

    if scaled_violation > Params.qp_feas_tol
        error(['solve_qp_checked: %s violates constraints; ', ...
               'scaled violation %.3e.'],label,scaled_violation);
    end

    objective = 0.5*x'*H*x+f'*x;
    initial_objective = 0.5*x0'*H*x0+f'*x0;

    objective_tol = 1e-8*max([ ...
        abs(objective),abs(initial_objective),objective_scale ...
    ]);

    if objective > initial_objective+objective_tol
        error('solve_qp_checked: %s worsened its feasible warm start.',label);
    end

    Report = struct();

    Report.label = label;
    Report.status = info.info;
    Report.iterations = info.solveiter;
    Report.elapsed = elapsed;

    Report.objective = objective;
    Report.initial_objective = initial_objective;

    Report.max_scaled_violation = scaled_violation;
    Report.max_raw_violation = raw_violation;

    Report.objective_scale = objective_scale;
    Report.lambda_scaled = lambda;
end
