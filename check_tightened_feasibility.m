function Feasibility = check_tightened_feasibility(Model, Tube, Params)
% =========================================================================
% check_tightened_feasibility
%
% Linear feasibility check for the complete fixed nominal trajectory.
%
% No tracking objective and no batch-update restrictions are imposed.
% This function does not depend on Stage 4 lifting/optimization files.
% =========================================================================

    A = Model.barA;
    B = Model.barB;
    C = Model.barC;

    nx = Model.nx;
    nu = Model.A_lv;

    H = Params.Ble-Params.ctrl_start;
    nvar = nu*H;

    ni = Params.n_in_opt;
    no = Params.n_out;

    Sv = diag(Model.X_std);
    Sy = diag(Model.Y_std);

    input_map = Sv*Model.P;
    output_map = Sy*Model.Q*C;

    nv_constraints = size(Tube.H_v,1);
    ny_constraints = size(Tube.H_y,1);

    nrows = H*(nv_constraints+ny_constraints);

    Aineq = zeros(nrows,nvar);
    bineq = zeros(nrows,1);

    state_map = zeros(nx,nvar);
    state_offset = Model.x0_ctrl;

    row = 0;

    for k = 1:H
        input_column = Params.ctrl_start+k;
        output_column = input_column+1;

        decision_columns = (k-1)*nu+(1:nu);

        % Input at current interval.
        vmap = zeros(ni,nvar);
        vmap(:,decision_columns) = input_map;

        voffset = ...
            Params.u_base(Model.input_idx,input_column) + Model.X_mean;

        % Predict state/output at the END of the interval.
        state_map = A*state_map;
        state_map(:,decision_columns) = ...
            state_map(:,decision_columns)+B;

        state_offset = A*state_offset;

        ymap = output_map*state_map;

        yoffset = ...
            Params.y_base(:,output_column) + ...
            Model.Y_mean + output_map*state_offset;

        rv = row+(1:nv_constraints);

        Aineq(rv,:) = Tube.H_v*vmap;
        bineq(rv) = Tube.h_v_tight-Tube.H_v*voffset;

        row = row+nv_constraints;

        ry = row+(1:ny_constraints);

        Aineq(ry,:) = Tube.H_y*ymap;
        bineq(ry) = Tube.h_y_tight-Tube.H_y*yoffset;

        row = row+ny_constraints;
    end

    row_scale = max(abs(Aineq),[],2);
    zero_rows = row_scale == 0;

    Feasibility = struct();

    Feasibility.feasible = false;
    Feasibility.n_variables = nvar;
    Feasibility.n_constraints = nrows;
    Feasibility.errnum = NaN;
    Feasibility.status = NaN;
    Feasibility.max_scaled_violation = Inf;
    Feasibility.t_feasible = [];

    if any(bineq(zero_rows) < 0)
        Feasibility.message = 'An input-independent constraint is infeasible.';
        return;
    end

    row_scale(zero_rows) = 1;

    As = bsxfun(@rdivide,Aineq,row_scale);
    bs = bineq./row_scale;

    As = As(~zero_rows,:);
    bs = bs(~zero_rows);

    options.msglev = 0;
    options.presol = 1;

    [t,~,errnum,extra] = glpk( ...
        zeros(nvar,1),As,bs, ...
        -Inf(nvar,1),Inf(nvar,1), ...
        repmat('U',1,size(As,1)), ...
        repmat('C',1,nvar),1,options);

    Feasibility.errnum = errnum;

    if isfield(extra,'status')
        Feasibility.status = extra.status;
    end

    if errnum == 0 && ...
       any(Feasibility.status == [2,5]) && ...
       ~isempty(t) && all(isfinite(t))

        violation = max([0;As*t-bs]);

        Feasibility.max_scaled_violation = violation;
        Feasibility.feasible = violation <= Params.qp_feas_tol;

        if Feasibility.feasible
            Feasibility.t_feasible = t;
            Feasibility.message = 'Complete tightened nominal trajectory is feasible.';
        else
            Feasibility.message = 'LP solution violates the scaled constraints.';
        end
    else
        Feasibility.message = sprintf( ...
            'GLPK did not certify feasibility: errnum=%d, status=%g.', ...
            errnum,Feasibility.status);
    end

    fprintf('\nTightened trajectory feasibility\n');
    fprintf('  %s\n',Feasibility.message);
    fprintf('  Variables: %d; constraints: %d.\n',nvar,nrows);
end
