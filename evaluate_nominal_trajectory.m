function Trajectory = evaluate_nominal_trajectory(Problem,t)
% =========================================================================
% evaluate_nominal_trajectory
%
% Evaluate one fixed nominal control-segment trajectory.
% All physical inputs/outputs are absolute physical quantities.
% =========================================================================

    t = t(:);

    if numel(t) ~= Problem.nt || any(~isfinite(t))
        error('evaluate_nominal_trajectory: invalid decision vector.');
    end

    Lift = Problem.Lift;
    H = Problem.H;

    Trajectory = struct();

    Trajectory.t_vector = t;
    Trajectory.T = reshape(t,Problem.nu,H);

    Trajectory.X_pre = reshape( ...
        Lift.Phi_x_pre*t+Lift.x_pre_offset,Problem.nx,H);

    Trajectory.X_post = reshape( ...
        Lift.Phi_x_post*t+Lift.x_post_offset,Problem.nx,H);

    Trajectory.V = reshape( ...
        Lift.Phi_v*t+Lift.v_offset,Problem.n_in,H);

    Trajectory.Y = reshape( ...
        Lift.Phi_y*t+Lift.y_offset,Problem.n_out,H);

    y_vector = Trajectory.Y(:);

    Trajectory.selected_outputs = y_vector(Problem.selected_rows);
    Trajectory.selected_references = Problem.selected_refs;

    Trajectory.tracking_residual = Problem.L*t+Problem.g;
    Trajectory.tracking_cost = ...
        Trajectory.tracking_residual'*Trajectory.tracking_residual;

    Trajectory.max_constraint_violation = ...
        max([0;Problem.A*t-Problem.b]);

    Trajectory.input_columns = Lift.input_columns;
    Trajectory.output_columns = Lift.output_columns;
end
