function PLS = fit_pls(X0, Y0, A_requested, Params)
% =========================================================================
% fit_pls
%
% Rank-aware static PLS with:
%   X deflation: X <- X - t*p'
%   Y deflation: Y <- Y - u*q'
%
% ARX fitting is performed only AFTER this function completes.
% =========================================================================

    if size(X0,1) ~= size(Y0,1)
        error('fit_pls: input/output sample counts differ.');
    end

    max_A = min([A_requested,size(X0,2),size(Y0,2)]);

    Xres = X0;
    Yres = Y0;

    T = [];
    U = [];
    W = [];
    C = [];
    P = [];
    Q = [];

    EVX = [];
    EVY = [];

    total_X = sum(X0(:).^2);
    total_Y = sum(Y0(:).^2);

    for a = 1:max_A
        [~,start_col] = max(sum(Yres.^2,1));
        u = Yres(:,start_col);

        accepted_iteration = false;

        for iter = 1:Params.nipals_max_iter
            wraw = Xres'*u;

            if norm(wraw) <= eps*max(1,norm(Xres,'fro')*norm(u))
                break;
            end

            w = wraw/norm(wraw);
            t = Xres*w;

            craw = Yres'*t;

            if norm(craw) <= eps*max(1,norm(Yres,'fro')*norm(t))
                break;
            end

            c = craw/norm(craw);
            unew = Yres*c;

            % Align all paired directions together.
            if unew'*u < 0
                w = -w;
                t = -t;
                c = -c;
                unew = -unew;
            end

            change = norm(unew-u)/max(norm(u),eps);
            u = unew;

            if change <= Params.nipals_tol
                accepted_iteration = true;
                break;
            end
        end

        if ~accepted_iteration
            break;
        end

        if t'*t <= eps || u'*u <= eps
            break;
        end

        p = Xres'*t/(t'*t);
        q = Yres'*u/(u'*u);

        Qcandidate = [Q,q];
        singular_values = svd(Qcandidate,'econ');

        if singular_values(end) < Params.pls_sigma_min || ...
           singular_values(end) <= ...
           Params.pls_rank_tol*singular_values(1)
            break;
        end

        T = [T,t];
        U = [U,u];

        W = [W,w];
        C = [C,c];

        P = [P,p];
        Q = Qcandidate;

        Xres = Xres-t*p';
        Yres = Yres-u*q';

        EVX(a) = 1-sum(Xres(:).^2)/total_X;
        EVY(a) = 1-sum(Yres(:).^2)/total_Y;
    end

    A_eff = size(P,2);

    if A_eff == 0
        error('fit_pls: no numerically valid PLS component was retained.');
    end

    % Right division solves the projection equations without inv().
    Rx = W/(P'*W);
    Ru = C/(Q'*C);

    input_inverse_error = norm(Rx'*P-eye(A_eff),'fro');
    output_inverse_error = norm(Ru'*Q-eye(A_eff),'fro');

    input_score_error = ...
        norm(X0*Rx-T,'fro')/max(norm(T,'fro'),eps);

    output_score_error = ...
        norm(Y0*Ru-U,'fro')/max(norm(U,'fro'),eps);

    if max([input_inverse_error,output_inverse_error, ...
            input_score_error,output_score_error]) > 1e-6
        error('fit_pls: score projection consistency check failed.');
    end

    PLS = struct();

    PLS.A_requested = A_requested;
    PLS.A_lv = A_eff;

    PLS.T = T;
    PLS.U = U;

    PLS.W = W;
    PLS.C = C;
    PLS.P = P;
    PLS.Q = Q;

    PLS.R_x = Rx;
    PLS.R_u = Ru;

    PLS.E_X = Xres;
    PLS.E_Y = Yres;

    PLS.EV_X = EVX;
    PLS.EV_Y = EVY;

    PLS.input_inverse_error = input_inverse_error;
    PLS.output_inverse_error = output_inverse_error;

    PLS.input_score_error = input_score_error;
    PLS.output_score_error = output_score_error;

    PLS.sigma_min_Q = min(svd(Q));
end
