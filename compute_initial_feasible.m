function Initial = compute_initial_feasible(Problem,Tube,Params)
% =========================================================================
% compute_initial_feasible
%
% Select a reproducible low-norm initial feasible trajectory.
% Uses the Stage 3 LP certificate as a feasible warm start.
% =========================================================================

    if ~isfield(Tube,'nominal_feasible') || ~Tube.nominal_feasible || ...
       ~isfield(Tube,'feasible_seed') || isempty(Tube.feasible_seed)
        error('compute_initial_feasible: Stage 3 feasible seed is missing.');
    end

    seed = Tube.feasible_seed(:);

    if numel(seed) ~= Problem.nt
        error('compute_initial_feasible: incompatible feasible-seed size.');
    end

    H = 2*eye(Problem.nt);
    f = zeros(Problem.nt,1);

    [t,Report] = solve_qp_checked( ...
        H,f,Problem.A,Problem.b,seed,Params,'initial feasible trajectory');

    Initial = struct();

    Initial.t = t;
    Initial.Trajectory = evaluate_nominal_trajectory(Problem,t);
    Initial.QP = Report;

    fprintf('Initial feasible trajectory: tracking cost %.8e.\n', ...
        Initial.Trajectory.tracking_cost);
end
