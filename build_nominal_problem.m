function Problem = build_nominal_problem(Model, Tube, Params)
% =========================================================================
% build_nominal_problem
%
% Fixed nominal optimization:
%   minimize ||L*t+g||^2
%   subject to A*t <= b.
%
% All offsets include:
%   physical baseline;
%   training mean;
%   fixed nominal initial-state contribution.
% =========================================================================

    if ~isfield(Tube,'state_convention') || ...
       ~strcmp(Tube.state_convention,'observed_score_history')
        error('build_nominal_problem: incompatible Tube state convention.');
    end

    if Params.A_lv ~= Model.A_lv || Params.n_arx ~= Model.n_arx
        error('build_nominal_problem: inconsistent model dimensions.');
    end

    Lift = build_lifting_matrices(Model,Params);

    H = Lift.H;
    nt = Lift.nt;
    nu = Model.A_lv;
    no = Params.n_out;

    %% 1. Fixed physical constraints
    Hv = kron(speye(H),Tube.H_v);
    Hy = kron(speye(H),Tube.H_y);

    Aineq = [ ...
        Hv*Lift.Phi_v;
        Hy*Lift.Phi_y ...
    ];

    bineq = [ ...
        repmat(Tube.h_v_tight,H,1)-Hv*Lift.v_offset;
        repmat(Tube.h_y_tight,H,1)-Hy*Lift.y_offset ...
    ];

    %% 2. Reference matrix
    refs = NaN(no,H);

    if isfield(Params,'y_ref_matrix') && ...
       ~isempty(Params.y_ref_matrix)

        refs = Params.y_ref_matrix;

        if ~isequal(size(refs),[no,H])
            error('build_nominal_problem: invalid y_ref_matrix dimensions.');
        end

    else
        kp = Params.Kp(:);
        ref = Params.y_ref_phys(:);
        tr = Params.track_idx;

        if numel(kp) ~= numel(ref) || ...
           numel(unique(kp)) ~= numel(kp)
            error('build_nominal_problem: invalid tracking specification.');
        end

        if ~isscalar(tr) || tr < 1 || tr > no || tr ~= fix(tr)
            error('build_nominal_problem: invalid track_idx.');
        end

        for j = 1:numel(kp)
            rel = kp(j)-Params.ctrl_start;

            if ~isfinite(rel) || rel ~= fix(rel) || ...
               rel < 1 || rel > H
                error('build_nominal_problem: invalid tracking time.');
            end

            refs(tr,rel) = ref(j);
        end
    end

    if any(isinf(refs(:))) || ~any(isfinite(refs(:)))
        error('build_nominal_problem: invalid or empty reference matrix.');
    end

    %% 3. Weighted masked objective
    Q = Params.Q_weight;

    if isscalar(Q)
        Q = Q*eye(no);
    end

    if ~isequal(size(Q),[no,no]) || any(~isfinite(Q(:)))
        error('build_nominal_problem: invalid tracking weight.');
    end

    if norm(Q-Q','fro') > 1e-10*max(1,norm(Q,'fro'))
        error('build_nominal_problem: tracking weight must be symmetric.');
    end

    Q = (Q+Q')/2;

    [~,flag] = chol(Q);

    if flag ~= 0
        error('build_nominal_problem: tracking weight must be positive definite.');
    end

    L = zeros(0,nt);
    g = zeros(0,1);

    selected_rows = zeros(0,1);
    selected_refs = zeros(0,1);

    for k = 1:H
        active = find(isfinite(refs(:,k)));

        if isempty(active)
            continue;
        end

        rows = (k-1)*no+active;
        W = chol(Q(active,active));

        L = [L;W*Lift.Phi_y(rows,:)];
        g = [g;W*(Lift.y_offset(rows)-refs(active,k))];

        selected_rows = [selected_rows;rows(:)];
        selected_refs = [selected_refs;refs(active,k)];
    end

    %% 4. Time differences of the batch update
    Dt = eye(H);

    if H > 1
        Dt = Dt-diag(ones(H-1,1),-1);
    end

    D = kron(sparse(Dt),speye(nu));

    R = Params.R_weight;

    if isscalar(R)
        R = R*eye(nu);
    end

    if ~isequal(size(R),[nu,nu]) || any(~isfinite(R(:))) || ...
       norm(R-R','fro') > 1e-10*max(1,norm(R,'fro'))
        error('build_nominal_problem: invalid proximal weight.');
    end

    R = (R+R')/2;

    [~,flag] = chol(R);

    if flag ~= 0
        error('build_nominal_problem: proximal weight must be positive definite.');
    end

    Rbar = kron(speye(H),R);
    Pprox = D'*Rbar*D;

    if ~isscalar(Params.delta_t_max) || ...
       ~isfinite(Params.delta_t_max) || Params.delta_t_max <= 0 || ...
       ~isscalar(Params.ddt_max) || ...
       ~isfinite(Params.ddt_max) || Params.ddt_max <= 0
        error('build_nominal_problem: invalid update bounds.');
    end

    Problem = struct();

    Problem.Lift = Lift;

    Problem.H = H;
    Problem.nt = nt;
    Problem.nu = nu;
    Problem.nx = Model.nx;
    Problem.n_in = Params.n_in_opt;
    Problem.n_out = no;

    Problem.A = Aineq;
    Problem.b = bineq;

    Problem.L = L;
    Problem.g = g;

    Problem.D = D;
    Problem.Rbar = Rbar;
    Problem.Pprox = Pprox;

    Problem.delta_t_max = Params.delta_t_max;
    Problem.ddt_max = Params.ddt_max;

    Problem.reference_matrix = refs;
    Problem.selected_rows = selected_rows;
    Problem.selected_refs = selected_refs;
end
