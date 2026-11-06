function Scores = project_scores(Standard, PLS)
% =========================================================================
% project_scores
%
% Project complete batches with fixed PLS projection matrices.
%
% T: N x A x B
% U: (N+1) x A x B
% R_y: standardized output reconstruction residual
% =========================================================================

    N = size(Standard.X0,1);
    Ny = size(Standard.Y0,1);
    B = Standard.N_batches;
    A = PLS.A_lv;

    Scores = struct();

    Scores.T = zeros(N,A,B);
    Scores.U = zeros(Ny,A,B);
    Scores.R_y = zeros(size(Standard.Y0));
    Scores.N_batches = B;

    for b = 1:B
        Scores.T(:,:,b) = Standard.X0(:,:,b)*PLS.R_x;
        Scores.U(:,:,b) = Standard.Y0(:,:,b)*PLS.R_u;

        Scores.R_y(:,:,b) = ...
            Standard.Y0(:,:,b)-Scores.U(:,:,b)*PLS.Q';
    end
end
