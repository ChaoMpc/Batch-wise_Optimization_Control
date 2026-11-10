function Optimal = compute_optimal_reachable( ...
    Problem,initial_t,Params)
% =========================================================================
% compute_optimal_reachable
%
% Solve the fixed limiting nominal problem:
%   min ||L*t+g||^2
%   s.t. A*t <= b.
%
% No proximal penalty or artificial Hessian regularization is added.
% =========================================================================

    L = Problem.L;
    g = Problem.g;

    H = 2*(L'*L);
    f = 2*L'*g;

    [t,Report] = solve_qp_checked( ...
        H,f,Problem.A,Problem.b,initial_t,Params, ...
        'optimal reachable reference');

    Trajectory = evaluate_nominal_trajectory(Problem,t);

    Optimal = struct();

    Optimal.t = t;
    Optimal.Trajectory = Trajectory;
    Optimal.QP = Report;

    Optimal.cost = Trajectory.tracking_cost;

    Optimal.selected_outputs = Trajectory.selected_outputs;
    Optimal.selected_references = Trajectory.selected_references;

    Optimal.reference_error = ...
        Optimal.selected_outputs-Optimal.selected_references;

    fprintf('\nOptimal reachable selected outputs\n');

    for j = 1:numel(Optimal.selected_outputs)
        fprintf('  Component %d: target %.9g; optimum %.9g; error %+.3e.\n', ...
            j,Optimal.selected_references(j), ...
            Optimal.selected_outputs(j),Optimal.reference_error(j));
    end

    fprintf('  Optimal weighted tracking cost: %.8e.\n',Optimal.cost);
end
