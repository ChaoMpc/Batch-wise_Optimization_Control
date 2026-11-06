function Data = run_oh_ilmpc( ...
    Problem, Initial, Optimal, Scenario, Model, Tube, UncSets, Params)
% =========================================================================
% run_oh_ilmpc
%
% Latent-variable point-to-point ILMPC comparison controller.
% =========================================================================

    %#ok<INUSD>

    %% =====================================================================
    % 1. Required fields
    % ======================================================================

    required_params = { ...
        'N_batches', 'Ble', 'ctrl_start', 'n_out', 'n_in', ...
        'u_base', 'y_base', 'u_phys_min', 'u_phys_max', ...
        'y_phys_min', 'y_phys_max', 'Kp', 'track_idx', 'y_ref_phys'};

    for ii = 1:numel(required_params)
        if ~isfield(Params,required_params{ii})
            error('run_oh_ilmpc: Params.%s is required.', required_params{ii});
        end
    end

    required_model = { ...
        'A_lv', 'n_arx', 'nx', 'barA', 'barB', 'barC', ...
        'X_mean', 'X_std', 'Y_mean', 'Y_std', ...
        'P', 'Q', 'R_x', 'R_u', 'input_idx', 'n_state_per_lv'};

    for ii = 1:numel(required_model)
        if ~isfield(Model,required_model{ii})
            error('run_oh_ilmpc: Model.%s is required.', required_model{ii});
        end
    end

    if ~isfield(Scenario,'X0') || ~isfield(Scenario,'noise') || ~isfield(Scenario,'state_jump')
        error(['run_oh_ilmpc: Scenario must contain X0, ', 'noise, and state_jump.']);
    end

    %% =====================================================================
    % 2. Basic dimensions
    % ======================================================================

    N       = Params.Ble;
    Batches = Params.N_batches;
    c0      = Params.ctrl_start;

    no  = Params.n_out;
    ni  = Params.n_in;
    nlv = Model.A_lv;
    nx  = Model.nx;

    H = N-c0;

    input_idx = Model.input_idx(:);
    nopt = numel(input_idx);

    if c0 < 1 || c0 >= N
        error('run_oh_ilmpc: ctrl_start must satisfy 1 <= ctrl_start < Ble.');
    end

    if H < 1
        error('run_oh_ilmpc: no controlled interval is available.');
    end

    if size(Params.u_base,1) ~= ni || size(Params.u_base,2) < N
        error('run_oh_ilmpc: invalid u_base dimensions.');
    end

    if size(Params.y_base,1) ~= no || size(Params.y_base,2) < N+1
        error('run_oh_ilmpc: invalid y_base dimensions.');
    end

    if numel(Params.u_phys_min) ~= ni || numel(Params.u_phys_max) ~= ni
        error('run_oh_ilmpc: invalid physical input bounds.');
    end

    if numel(Params.y_phys_min) ~= no || numel(Params.y_phys_max) ~= no
        error('run_oh_ilmpc: invalid physical output bounds.');
    end

    if any(input_idx < 1) || any(input_idx > ni)
        error('run_oh_ilmpc: Model.input_idx is out of range.');
    end

    if size(Model.barA,1) ~= nx || size(Model.barA,2) ~= nx
        error('run_oh_ilmpc: barA dimension mismatch.');
    end

    if size(Model.barB,1) ~= nx || size(Model.barB,2) ~= nlv
        error('run_oh_ilmpc: barB dimension mismatch.');
    end

    if size(Model.barC,1) ~= nlv || size(Model.barC,2) ~= nx
        error('run_oh_ilmpc: barC dimension mismatch.');
    end

    if size(Model.P,1) ~= nopt || size(Model.P,2) ~= nlv
        error('run_oh_ilmpc: P dimension mismatch.');
    end

    if size(Model.Q,1) ~= no || size(Model.Q,2) ~= nlv
        error('run_oh_ilmpc: Q dimension mismatch.');
    end

    Xmean = Model.X_mean(:);
    Xstd  = Model.X_std(:);
    Ymean = Model.Y_mean(:);
    Ystd  = Model.Y_std(:);

    if numel(Xmean) ~= nopt || numel(Xstd) ~= nopt
        error('run_oh_ilmpc: X_mean/X_std dimension mismatch.');
    end

    if numel(Ymean) ~= no || numel(Ystd) ~= no
        error('run_oh_ilmpc: Y_mean/Y_std dimension mismatch.');
    end

    if any(abs(Xstd) < eps) || any(abs(Ystd) < eps)
        error('run_oh_ilmpc: X_std or Y_std contains zero entries.');
    end

    expected_nx = nlv*(2*Model.n_arx-1);

    if nx ~= expected_nx
        error(['run_oh_ilmpc: Model.nx=%d is inconsistent with ', ...
               'A_lv=%d and n_arx=%d. Expected %d.'], ...
            nx,nlv,Model.n_arx,expected_nx);
    end

    if Model.n_state_per_lv ~= 2*Model.n_arx-1
        error('run_oh_ilmpc: n_state_per_lv is inconsistent with n_arx.');
    end

    %% =====================================================================
    % 3. QP parameters
    % ======================================================================

    prediction_horizon = 60;

    if isfield(Params,'oh_prediction_horizon') && ~isempty(Params.oh_prediction_horizon)
        prediction_horizon = round(Params.oh_prediction_horizon);
    end
    prediction_horizon = max(1,min(H,prediction_horizon));

    move_horizon = prediction_horizon;
    if isfield(Params,'oh_control_horizon') && ~isempty(Params.oh_control_horizon)
        move_horizon = round(Params.oh_control_horizon);
    end
    move_horizon = max(1,min(prediction_horizon,move_horizon));

    q_track = get_scalar_parameter(Params,'Q_weight',1);
    r_delta = get_scalar_parameter(Params,'R_weight',1e-6);
    p_delta = get_scalar_parameter(Params,'P_weight',1e-6);
    eps_weight = get_scalar_parameter(Params,'E_weight',1e2);

    if q_track < 0 || r_delta < 0 || p_delta < 0 || eps_weight <= 0
        error('run_oh_ilmpc: invalid QP weights.');
    end

    qp_feas_tol = 1e-8;
    if isfield(Params,'qp_feas_tol') && ~isempty(Params.qp_feas_tol)
        qp_feas_tol = Params.qp_feas_tol;
    end

    %% =====================================================================
    % 4. Latent/physical mappings
    % ======================================================================

    Mv = diag(Xstd)*Model.P;
    My = diag(Ystd)*Model.Q;

    if size(Mv,1) ~= nopt || size(Mv,2) ~= nlv
        error('run_oh_ilmpc: Mv dimension mismatch.');
    end
    if size(My,1) ~= no || size(My,2) ~= nlv
        error('run_oh_ilmpc: My dimension mismatch.');
    end

    %% =====================================================================
    % 5. Augmented incremental model
    % ======================================================================

    Ainc = [Model.barA, Model.barB; zeros(nlv,nx), eye(nlv)];
    Binc = [Model.barB; eye(nlv)];
    Cinc = [Model.barC, zeros(nlv,nlv)];
    nz = nx+nlv;

    if size(Ainc,1) ~= nz || size(Ainc,2) ~= nz
        error('run_oh_ilmpc: Ainc dimension mismatch.');
    end

    %% =====================================================================
    % 6. Physical box constraints
    % ======================================================================

    umin = Params.u_phys_min(input_idx);
    umax = Params.u_phys_max(input_idx);

    Hv = [eye(nopt); -eye(nopt)];
    hv = [umax; -umin];

    Hy = [eye(no); -eye(no)];
    hy = [Params.y_phys_max(:); -Params.y_phys_min(:)];

    %% =====================================================================
    % 7. Data allocation
    % ======================================================================

    Data = struct();
    Data.method = 'Latent PTP-ILMPC with consistent observer and constraints';

    Data.Y_hist = zeros(no,N+1,Batches);
    Data.Y_true_hist = zeros(no,N+1,Batches);
    Data.V_hist = zeros(ni,N,Batches);
    Data.V_request_hist = zeros(ni,N,Batches);
    Data.X_hist = zeros(size(Scenario.X0,1),N+1,Batches);
    Data.Y_nominal_hist = zeros(no,N+1,Batches);
    Data.V_nominal_hist = zeros(ni,N,Batches);
    Data.T_hist = zeros(nlv,H,Batches);
    Data.T_request_hist = zeros(nlv,H,Batches);
    Data.T_applied_hist = zeros(nlv,H,Batches);
    Data.delta_t_hist = zeros(nlv,H,Batches);
    Data.delta_t_norm = zeros(Batches,H);
    Data.input_change_norm = zeros(Batches,H);
    Data.QP_times = zeros(Batches,H);
    Data.tube_pre_quadratic = zeros(Batches,H);
    Data.tube_post_quadratic = zeros(Batches,H);
    Data.Err_true_act = NaN(Batches,1);
    Data.Err_opt_act = NaN(Batches,1);
    Data.Err_opt_nom = NaN(Batches,1);
    Data.Batch = cell(Batches,1);

    %% =====================================================================
    % 8. Initial previous trajectory
    % ======================================================================

    t_base = Model.R_x.'*((-Xmean)./Xstd);
    t_base = t_base(:);

    T_previous = repmat(t_base,1,N);
    Y_previous = Params.y_base;

    if isfield(Initial,'Trajectory') && isfield(Initial.Trajectory,'Y') && ~isempty(Initial.Trajectory.Y)
        ycols = get_output_columns(Problem,c0,H);
        Yinitial = Initial.Trajectory.Y;
        if size(Yinitial,1) ~= no || size(Yinitial,2) ~= numel(ycols)
            error('run_oh_ilmpc: Initial.Trajectory.Y dimensions are inconsistent.');
        end
        Y_previous(:,ycols) = Yinitial;
    end

    if isfield(Initial,'Trajectory') && isfield(Initial.Trajectory,'V') && ~isempty(Initial.Trajectory.V)
        vcols = get_input_columns(Problem,c0,H);
        Vinitial = Initial.Trajectory.V;
        if size(Vinitial,1) ~= nopt || size(Vinitial,2) ~= numel(vcols)
            error('run_oh_ilmpc: Initial.Trajectory.V dimensions are inconsistent.');
        end
        for jj = 1:numel(vcols)
            time_index = vcols(jj);
            vdev = Vinitial(:,jj) - Params.u_base(input_idx,time_index);
            T_previous(:,time_index) = physical_input_to_score(vdev,Model);
        end
    end

    %% =====================================================================
    % 9. Batch loop
    % ======================================================================

    total_timer = tic();

    fprintf('\n============================================================\n');
    fprintf('Starting adapted Oh ILMPC simulation\n');
    fprintf('Batches            : %d\n',Batches);
    fprintf('Control steps/batch: %d\n',H);
    fprintf('Prediction horizon : %d\n',prediction_horizon);
    fprintf('Control horizon    : %d\n',move_horizon);
    fprintf('============================================================\n');

    for batch = 1:Batches
        batch_timer = tic();

        % (基线仿真与数据初始化省略，保持原样)
        Y = zeros(no,N+1); Ytrue = zeros(no,N+1);
        V = zeros(ni,N); Vrequest = zeros(ni,N);
        X = zeros(size(Scenario.X0,1),N+1);
        Tcurrent = zeros(nlv,N); Trequest = zeros(nlv,H);
        Tapplied = zeros(nlv,H); deltaTcurrent = zeros(nlv,H);

        for j = 1:c0
            Tcurrent(:,j) = physical_input_to_score(zeros(nopt,1),Model);
        end

        x_phys = Scenario.X0(:,batch);
        X(:,1) = x_phys;
        [Y(:,1),Ytrue(:,1)] = cryst_measure(x_phys,Params,true,Scenario.noise(:,1,batch));

        for j = 1:c0
            x_phys = x_phys+Scenario.state_jump(:,j,batch);
            v = Params.u_base(:,j);
            V(:,j) = v; Vrequest(:,j) = v;
            [Y(:,j+1),x_phys,Ytrue(:,j+1)] = cryst_plant(x_phys,v,Params,true,Scenario.noise(:,j+1,batch));
            X(:,j+1) = x_phys;
            vdev = v(input_idx)-Params.u_base(input_idx,j);
            Tcurrent(:,j) = physical_input_to_score(vdev,Model);
        end

        z_current = make_augmented_state_from_trajectory(Y,Tcurrent,c0+1,Model,Params);
        last_qp_report = struct();

        %% Controlled intervals
        for kk = 1:H
            j = c0+kk;
            p = min(prediction_horizon,N-j+1);
            m = min(move_horizon,p);
            if p < 1, break; end
            nd = nlv*m;

            z_current = make_augmented_state_from_trajectory(Y,Tcurrent,j,Model,Params);
            z_previous = make_augmented_state_from_trajectory(Y_previous,T_previous,j,Model,Params);
            dz0 = z_current-z_previous;

            [Fz,Gz] = build_prediction_matrices(Ainc,Binc,p,m);
            CuF = kron(eye(p),Cinc)*Fz;
            CuG = kron(eye(p),Cinc)*Gz;

            Uprevious = zeros(nlv,p);
            for q = 1:p
                output_time = j+q;
                Uprevious(:,q) = physical_output_to_score(Y_previous(:,output_time), Params.y_base(:,output_time), Model);
            end

            u_previous_vec = Uprevious(:);
            u_offset = u_previous_vec+CuF*dz0;

            y_offset = zeros(no*p,1);
            for q = 1:p
                rows_y = (q-1)*no+(1:no);
                rows_u = (q-1)*nlv+(1:nlv);
                output_time = j+q;
                y_offset(rows_y) = Params.y_base(:,output_time) + Ymean + My*u_offset(rows_u);
            end

            Gy_phys = kron(eye(p),My)*CuG;
            [Zref,rref] = make_ptp_selector(Params,j,p,no);

            if isempty(Zref)
                Htrack = zeros(nd,nd); ftrack = zeros(nd,1);
            else
                Zref = full(Zref); rref = rref(:);
                Gref = Zref*Gy_phys;
                e0 = Zref*y_offset-rref;
                Htrack = q_track*(Gref.'*Gref);
                ftrack = q_track*(Gref.'*e0);
            end

            Hregularization = kron(eye(m), (r_delta+p_delta)*eye(nlv));
            Hdecision = Htrack+Hregularization;
            Hdecision = 0.5*(Hdecision+Hdecision.');
            fdecision = ftrack(:);

            TdiffF = zeros(p*nlv,nz);
            TdiffG = zeros(p*nlv,nd);
            for q = 1:p
                rows_z = (q-1)*nz+(1:nz);
                rows_t = (q-1)*nlv+(1:nlv);
                TdiffF(rows_t,:) = Fz(rows_z(nx+1:nz),:);
                TdiffG(rows_t,:) = Gz(rows_z(nx+1:nz),:);
            end

            Tprevious_future = T_previous(:,j:j+p-1);
            t_previous_vec = Tprevious_future(:);
            t_offset = t_previous_vec+TdiffF*dz0;

            Hlatent_v = Hv*Mv;
            nvrow = size(Hv,1);
            Ainput = zeros(p*nvrow,nd);
            binput = zeros(p*nvrow,1);
            for q = 1:p
                rows_v = (q-1)*nvrow+(1:nvrow);
                rows_t = (q-1)*nlv+(1:nlv);
                input_time = j+q-1;
                vbase_selected = Params.u_base(input_idx,input_time);
                hlatent_v = hv - Hv*vbase_selected - Hv*Xmean;
                Ainput(rows_v,:) = Hlatent_v*TdiffG(rows_t,:);
                binput(rows_v) = hlatent_v - Hlatent_v*t_offset(rows_t);
            end

            Hlatent_y = Hy*My;
            nyrow = size(Hy,1);
            neps = p*nyrow;
            Aoutput = zeros(p*nyrow,nd);
            boutput = zeros(p*nyrow,1);
            for q = 1:p
                rows_h = (q-1)*nyrow+(1:nyrow);
                rows_u = (q-1)*nlv+(1:nlv);
                output_time = j+q;
                ybase_q = Params.y_base(:,output_time);
                hlatent_y = hy - Hy*ybase_q - Hy*Ymean;
                Aoutput(rows_h,:) = Hlatent_y*CuG(rows_u,:);
                boutput(rows_h) = hlatent_y - Hlatent_y*u_offset(rows_u);
            end

            Aineq = [Ainput, zeros(size(Ainput,1),neps);
                     Aoutput, -eye(neps);
                     zeros(neps,nd), -eye(neps)];
            bineq = [binput; boutput; zeros(neps,1)];

            d0 = zeros(nd,1);
            dz_input0 = dz0(nx+1:nz);
            d0(1:nlv) = -dz_input0;

            eps0 = zeros(neps,1);
            y_warm = y_offset+Gy_phys*d0;
            for q = 1:p
                rows_y = (q-1)*no+(1:no);
                rows_h = (q-1)*nyrow+(1:nyrow);
                y0 = y_warm(rows_y);
                eps0(rows_h) = max(Hy*y0-hy,0)+1e-9;
            end
            x0 = [d0;eps0];

            input_warm_violation = max([0;Ainput*d0-binput]);
            if input_warm_violation > qp_feas_tol
                error('run_oh_ilmpc: no feasible input warm start at batch %d, time %d. Violation %.6e.', batch,j,input_warm_violation);
            end

            %% Solve QP
            Hqp = blkdiag(Hdecision+1e-10*eye(nd), eps_weight*eye(neps));
            Hqp = 0.5*(Hqp+Hqp.');
            fqp = [fdecision; zeros(neps,1)];

            timer_qp = tic();
            [solution,qp_report] = solve_qp_checked(Hqp, fqp, Aineq, bineq, x0, Params, 'adapted Oh comparison');
            Data.QP_times(batch,kk) = toc(timer_qp);
            last_qp_report = qp_report;

            % =========================================================
            % 新增：同时打印当前步 QP 时间和当前批次累积总时间
            % =========================================================
            fprintf('  Oh: Batch %2d/%2d | Step %2d/%2d | QP: %6.2fs | Accum: %6.1fs\n', ...
                batch, Batches, kk, H, Data.QP_times(batch,kk), toc(batch_timer));
            fflush(stdout);
            % =========================================================

            solution = solution(:);
            delta_sequence = solution(1:nd);
            delta_now = delta_sequence(1:nlv);
            Data.delta_t_norm(batch,kk) = norm(delta_now,inf);

            t_current_previous = Tcurrent(:,j-1);
            dt_previous = T_previous(:,j)-T_previous(:,j-1);
            t_request = t_current_previous + dt_previous + delta_now;
            t_request = t_request(:);

            v_request = Params.u_base(:,j);
            v_request(input_idx) = Params.u_base(input_idx,j) + Xmean + Mv*t_request;

            input_change = norm(v_request(input_idx) - Params.u_base(input_idx,j),inf);
            Data.input_change_norm(batch,kk) = input_change;

            raw_input_violation = max([0; v_request(input_idx)-umax; umin-v_request(input_idx)]);
            if raw_input_violation > 1e-7
                error('run_oh_ilmpc: requested physical input violates bounds by %.6e at batch %d, time %d.', raw_input_violation,batch,j);
            end

            v_apply = v_request;
            Vrequest(:,j) = v_request;
            V(:,j) = v_apply;
            Trequest(:,kk) = t_request;

            x_phys = x_phys+Scenario.state_jump(:,j,batch);
            [Y(:,j+1),x_phys,Ytrue(:,j+1)] = cryst_plant(x_phys, v_apply, Params, true, Scenario.noise(:,j+1,batch));
            X(:,j+1) = x_phys;

            vdev_applied = v_apply(input_idx)-Params.u_base(input_idx,j);
            t_applied = physical_input_to_score(vdev_applied,Model);
            Tapplied(:,kk) = t_applied;
            Tcurrent(:,j) = t_applied;

            deltaTcurrent(:,kk) = (t_applied-Tcurrent(:,j-1)) - (T_previous(:,j)-T_previous(:,j-1));
        end

        %% Store batch data
        Data.Y_hist(:,:,batch) = Y;
        Data.Y_true_hist(:,:,batch) = Ytrue;
        Data.V_hist(:,:,batch) = V;
        Data.V_request_hist(:,:,batch) = Vrequest;
        Data.X_hist(:,:,batch) = X;
        Data.T_hist(:,:,batch) = Tcurrent(:,c0+1:N);
        Data.T_request_hist(:,:,batch) = Trequest;
        Data.T_applied_hist(:,:,batch) = Tapplied;
        Data.delta_t_hist(:,:,batch) = deltaTcurrent;

        Ynom = Params.y_base; Vnom = Params.u_base;
        if isfield(Initial,'Trajectory') && isfield(Initial.Trajectory,'Y')
            ycols = get_output_columns(Problem,c0,H);
            if size(Initial.Trajectory.Y,2) == numel(ycols)
                Ynom(:,ycols) = Initial.Trajectory.Y;
            end
        end
        if isfield(Initial,'Trajectory') && isfield(Initial.Trajectory,'V')
            vcols = get_input_columns(Problem,c0,H);
            if size(Initial.Trajectory.V,2) == numel(vcols)
                Vnom(input_idx,vcols) = Initial.Trajectory.V;
            end
        end
        Data.Y_nominal_hist(:,:,batch) = Ynom;
        Data.V_nominal_hist(:,:,batch) = Vnom;

        if isfield(Problem,'Lift') && isfield(Problem.Lift,'output_columns')
            output_columns = Problem.Lift.output_columns;
            yvector = Y(:,output_columns)(:);
            selected_rows = Problem.selected_rows(:);
            selected_refs = Problem.selected_refs(:);

            if ~isempty(selected_rows) && numel(selected_refs) == numel(selected_rows)
                Data.Err_true_act(batch) = norm(yvector(selected_rows)-selected_refs);
            end
            if isfield(Optimal,'selected_outputs') && ~isempty(Optimal.selected_outputs)
                Data.Err_opt_act(batch) = norm(yvector(selected_rows)-Optimal.selected_outputs(:));
            end
        end
        Data.Err_opt_nom(batch) = NaN;

        BatchRecord = struct();
        BatchRecord.method = Data.method;
        BatchRecord.qp_report = last_qp_report;
        BatchRecord.T_request = Trequest;
        BatchRecord.T_applied = Tapplied;
        BatchRecord.delta_t = deltaTcurrent;
        BatchRecord.clip_flags = false(nopt,H);
        BatchRecord.eta_box_in = true(1,H);
        BatchRecord.r_box_in = true(1,H);
        BatchRecord.eta_ellipsoid_in = true(1,H);
        BatchRecord.r_ellipsoid_in = true(1,H);
        BatchRecord.tube_pre_quadratic = zeros(1,H);
        BatchRecord.tube_post_quadratic = zeros(1,H);
        Data.tube_pre_quadratic(batch,:) = BatchRecord.tube_pre_quadratic;
        Data.tube_post_quadratic(batch,:) = BatchRecord.tube_post_quadratic;
        Data.Batch{batch} = BatchRecord;

        fprintf('\n>>> Oh batch %2d/%2d SUMMARY: tracking error %.4e | max delta-t %.4e | batch total time %.2f s\n\n', ...
            batch, Batches, Data.Err_true_act(batch), ...
            max(Data.delta_t_norm(batch,:)), toc(batch_timer));

        Y_previous = Y;
        T_previous = Tcurrent;
    end

    fprintf('\n============================================================\n');
    fprintf('Adapted Oh ILMPC simulation completed.\n');
    fprintf('Total simulation time: %.3f s\n',toc(total_timer));
    fprintf('============================================================\n');
end

% =========================================================================
% 以下为辅助函数
% =========================================================================
function value = get_scalar_parameter(Params,name,default_value)
    if ~isfield(Params,name) || isempty(Params.(name))
        value = default_value; return;
    end
    raw = Params.(name);
    if ~isscalar(raw), value = raw(1,1); else, value = raw; end
    if ~isfinite(value)
        error('run_oh_ilmpc: Params.%s is not finite.',name);
    end
end

function t = physical_input_to_score(vdev,Model)
    vdev = vdev(:); Xmean = Model.X_mean(:); Xstd = Model.X_std(:);
    if numel(vdev) ~= numel(Xmean)
        error('physical_input_to_score: dimension mismatch.');
    end
    vstd = (vdev-Xmean)./Xstd; t = Model.R_x.'*vstd; t = t(:);
    if any(~isfinite(t)), error('physical_input_to_score: nonfinite result.'); end
end

function u = physical_output_to_score(y,ybase,Model)
    y = y(:); ybase = ybase(:); Ymean = Model.Y_mean(:); Ystd = Model.Y_std(:);
    if numel(y) ~= numel(Ymean) || numel(ybase) ~= numel(Ymean)
        error('physical_output_to_score: dimension mismatch.');
    end
    ystd = (y-ybase-Ymean)./Ystd; u = Model.R_u.'*ystd; u = u(:);
    if any(~isfinite(u)), error('physical_output_to_score: nonfinite result.'); end
end

function z = make_augmented_state_from_trajectory(Y,T,j,Model,Params)
    n = Model.n_arx; nlv = Model.A_lv; block = Model.n_state_per_lv; nx = Model.nx;
    if size(Y,1) ~= Params.n_out, error('make_augmented_state_from_trajectory: Y row dimension mismatch.'); end
    if size(T,1) ~= nlv, error('make_augmented_state_from_trajectory: T row dimension mismatch.'); end
    if j < n, error('make_augmented_state_from_trajectory: time index is too early for ARX order.'); end
    if j > size(Y,2), error('make_augmented_state_from_trajectory: Y has insufficient columns.'); end
    if j-1 > size(T,2), error('make_augmented_state_from_trajectory: T has insufficient input history.'); end
    x = zeros(nx,1);
    for a = 1:nlv
        idx = (a-1)*block+(1:block);
        for ell = 1:n
            output_time = j-ell+1;
            u_score = physical_output_to_score(Y(:,output_time), Params.y_base(:,output_time), Model);
            x(idx(ell)) = u_score(a);
        end
        if n > 1
            for ell = 1:n-1
                input_time = j-ell;
                x(idx(n+ell)) = T(a,input_time);
            end
        end
    end
    t_previous = T(:,j-1);
    z = [x;t_previous];
end

function [F,G] = build_prediction_matrices(A,B,p,m)
    nz = size(A,1); nu = size(B,2);
    F = zeros(nz*p,nz); G = zeros(nz*p,nu*m);
    for q = 1:p
        rows = (q-1)*nz+(1:nz); F(rows,:) = A^q;
        for ell = 1:min(q,m)
            cols = (ell-1)*nu+(1:nu); G(rows,cols) = A^(q-ell)*B;
        end
    end
end

function [Z,r] = make_ptp_selector(Params,j,p,no)
    Z = zeros(0,no*p); r = zeros(0,1);
    Kp = Params.Kp(:); track_idx = Params.track_idx(:); yref = Params.y_ref_phys;
    if isempty(Kp) || isempty(track_idx), return; end
    for q = 1:p
        output_time = j+q; hit = find(Kp == output_time,1);
        if isempty(hit), continue; end
        for ii = 1:numel(track_idx)
            out = track_idx(ii);
            if out < 1 || out > no, error('make_ptp_selector: track_idx out of range.'); end
            refval = get_reference_value(yref,hit,ii,out,numel(Kp),numel(track_idx),no);
            row = zeros(1,no*p); row((q-1)*no+out) = 1;
            Z = [Z;row]; r = [r;refval];
        end
    end
    r = r(:);
end

function refval = get_reference_value(yref,hit,ii,out,nkp,ntrack,no)
    if isvector(yref)
        yref = yref(:);
        if numel(yref) == nkp, refval = yref(hit);
        elseif numel(yref) == nkp*ntrack, refval = yref((hit-1)*ntrack+ii);
        else, error('run_oh_ilmpc: y_ref_phys vector dimension is inconsistent.'); end
    elseif size(yref,1) == no && size(yref,2) == nkp, refval = yref(out,hit);
    elseif size(yref,1) == nkp && size(yref,2) == no, refval = yref(hit,out);
    else, error('run_oh_ilmpc: y_ref_phys dimensions are inconsistent.'); end
end

function cols = get_output_columns(Problem,c0,H)
    if isfield(Problem,'Lift') && isfield(Problem.Lift,'output_columns') && ~isempty(Problem.Lift.output_columns)
        cols = Problem.Lift.output_columns(:).';
    else
        cols = c0+2:c0+H+1;
    end
end

function cols = get_input_columns(Problem,c0,H)
    if isfield(Problem,'Lift') && isfield(Problem.Lift,'input_columns') && ~isempty(Problem.Lift.input_columns)
        cols = Problem.Lift.input_columns(:).';
    else
        cols = c0+1:c0+H;
    end
end
