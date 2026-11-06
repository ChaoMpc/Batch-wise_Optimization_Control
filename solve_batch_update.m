function Update = solve_batch_update(Problem,t_previous,Params)
% =========================================================================
% solve_batch_update
%
% Proximal nominal batch update:
%
% min_d ||L*(t_previous+d)+g||^2 + d'*Pprox*d
%
% subject to:
%   A*(t_previous+d) <= b
%   |d| <= delta_t_max
%   |D*d| <= ddt_max
%
% Zero update is feasible whenever t_previous is feasible.
% =========================================================================

    t_previous = t_previous(:);

    if numel(t_previous) ~= Problem.nt
        error('solve_batch_update: invalid previous decision dimensions.');
    end

    L = Problem.L;
    g = Problem.g;
    D = Problem.D;

    nt = Problem.nt;

    previous_residual = L*t_previous+g;

    H = 2*(L'*L+Problem.Pprox);
    f = 2*L'*previous_residual;

    Aupdate = [ ...
        Problem.A;
        eye(nt);
       -eye(nt);
        D;
       -D ...
    ];

    bupdate = [ ...
        Problem.b-Problem.A*t_previous;
        Problem.delta_t_max*ones(nt,1);
        Problem.delta_t_max*ones(nt,1);
        Problem.ddt_max*ones(size(D,1),1);
        Problem.ddt_max*ones(size(D,1),1) ...
    ];

    [d,Report] = solve_qp_checked( ...
        H,f,Aupdate,bupdate,zeros(nt,1),Params,'nominal batch update');

    t_new = t_previous+d;

    old_cost = previous_residual'*previous_residual;
    new_residual = L*t_new+g;
    new_cost = new_residual'*new_residual;

    proximal_cost = d'*Problem.Pprox*d;

    % Exact optimization implies:
    % J(t_new)+proximal_cost <= J(t_previous).
    numerical_tol = 1e-7*max([old_cost,new_cost,1e-12]);

    descent_residual = new_cost+proximal_cost-old_cost;

    if descent_residual > numerical_tol
        error(['solve_batch_update: proximal descent check failed; ', ...
               'residual %.3e.'],descent_residual);
    end

    Update = struct();

    Update.t = t_new;
    Update.d = d;
    Update.q = D*d;

    Update.old_cost = old_cost;
    Update.new_cost = new_cost;
    Update.proximal_cost = proximal_cost;

    Update.descent_residual = descent_residual;
    Update.step_norm = norm(d);
    Update.max_abs_step = max(abs(d));

    Update.QP = Report;
    Update.Trajectory = evaluate_nominal_trajectory(Problem,t_new);
end
