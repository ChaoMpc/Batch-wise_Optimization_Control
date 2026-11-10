function Feedback = design_tube_feedback(Model, Params)
% =========================================================================
% design_tube_feedback
%
% Full MIMO discrete-time LQR.
%
% dlqr returns the gain for input = -K_lqr*x.
% The Tube control convention is input = nominal + K*error,
% therefore K = -K_lqr.
% =========================================================================

    A = Model.barA;
    B = Model.barB;

    nx = size(A,1);
    nu = size(B,2);

    if isfield(Params,'Q_lqr_matrix') && ...
       ~isempty(Params.Q_lqr_matrix)
        Q = Params.Q_lqr_matrix;
    else
        Q = Params.Q_lqr*eye(nx);
    end

    if isfield(Params,'R_lqr_matrix') && ...
       ~isempty(Params.R_lqr_matrix)
        R = Params.R_lqr_matrix;
    else
        R = Params.R_lqr*eye(nu);
    end

    check_pd_local(Q,nx,'Q_lqr');
    check_pd_local(R,nu,'R_lqr');

    Q = (Q+Q')/2;
    R = (R+R')/2;

    try
        [K_lqr,S,~] = dlqr(A,B,Q,R);
    catch ME
        error('design_tube_feedback: full MIMO dlqr failed: %s',ME.message);
    end

    K = -K_lqr;
    AK = A+B*K;

    rho = max(abs(eig(AK)));

    if ~isfinite(rho) || rho >= 1-1e-10
        error(['design_tube_feedback: closed-loop matrix is not ', ...
               'numerically Schur stable; rho=%.12g.'],rho);
    end

    % Independent Riccati residual diagnostic.
    Riccati_residual = ...
        A'*S*A-S ...
        - A'*S*B*((R+B'*S*B)\(B'*S*A)) + Q;

    relative_residual = norm(Riccati_residual,'fro') / ...
                        max(1,norm(S,'fro'));

    if relative_residual > 1e-6
        error('design_tube_feedback: large Riccati residual %.3e.', ...
            relative_residual);
    end

    Feedback = struct();

    Feedback.K = K;
    Feedback.A_K = AK;

    Feedback.Q = Q;
    Feedback.R = R;
    Feedback.S = S;

    Feedback.spectral_radius = rho;
    Feedback.riccati_residual = relative_residual;

    fprintf('\nFull MIMO feedback\n');
    fprintf('  Gain size: %d x %d.\n',size(K,1),size(K,2));
    fprintf('  Closed-loop spectral radius: %.8f.\n',rho);
    fprintf('  Riccati relative residual: %.3e.\n',relative_residual);
end

function check_pd_local(M,n,label)
    if ~isequal(size(M),[n,n]) || any(~isfinite(M(:)))
        error('design_tube_feedback: invalid %s dimensions/data.',label);
    end

    if norm(M-M','fro') > 1e-10*max(1,norm(M,'fro'))
        error('design_tube_feedback: %s must be symmetric.',label);
    end

    [~,flag] = chol((M+M')/2);

    if flag ~= 0
        error('design_tube_feedback: %s must be positive definite.',label);
    end
end
