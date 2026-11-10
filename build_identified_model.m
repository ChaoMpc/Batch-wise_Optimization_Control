function Model = build_identified_model( ...
    Prep, Standard_train, A_requested, n_arx, Params)
% =========================================================================
% build_identified_model
%
% 组装：
%   训练预处理；
%   rank-aware 静态 PLS；
%   PLS 完成后统一 MIMO ARX；
%   状态空间实现；
%   固定 nominal 初始得分历史。
%
% 注意：ARX 不在 PLS deflation 过程中拟合，而是在全部潜变量得分
% 提取完成后统一拟合。
% =========================================================================

    PLS = fit_pls(Prep.X0, Prep.Y0, A_requested, Params);

    Scores = project_scores(Standard_train, PLS);

    % PLS 已完成，使用全部潜变量得分统一拟合 MIMO ARX。
    ARX = fit_latent_arx(Scores, n_arx, Params);

    SS = realize_latent_model(ARX);

    Model = struct();

    Model.code_version = Params.code_version;

    Model.A_requested = A_requested;
    Model.A_lv = PLS.A_lv;
    Model.n_arx = n_arx;
    Model.nx = SS.nx;

    Model.X_mean = Prep.X_mean;
    Model.X_std = Prep.X_std;
    Model.Y_mean = Prep.Y_mean;
    Model.Y_std = Prep.Y_std;

    Model.input_idx = Prep.input_idx;

    Model.W = PLS.W;
    Model.C = PLS.C;
    Model.P = PLS.P;
    Model.Q = PLS.Q;

    Model.R_x = PLS.R_x;
    Model.R_u = PLS.R_u;

    Model.Theta = ARX.Theta;
    Model.Theta_matrix = ARX.Theta_matrix;

    Model.barA = SS.barA;
    Model.barB = SS.barB;
    Model.barC = SS.barC;
    Model.barE = SS.barE;

    Model.n_state_per_lv = SS.n_state_per_lv;
    Model.open_loop_spectral_radius = SS.open_loop_spectral_radius;

    Model.PLS = PLS;
    Model.ARX = ARX;

    % 固定 nominal 初始状态（物理时间 ctrl_start）。
    t0 = Model.R_x.' * ((-Model.X_mean) ./ Model.X_std);
    u0 = Model.R_u.' * ((-Model.Y_mean) ./ Model.Y_std);

    block = Model.n_state_per_lv;
    x0 = zeros(Model.nx, 1);

    for a = 1:Model.A_lv
        idx = (a-1)*block + (1:block);

        x0(idx(1:n_arx)) = u0(a) * ones(n_arx, 1);

        if n_arx > 1
            x0(idx(n_arx+1:end)) = t0(a) * ones(n_arx-1, 1);
        end
    end

    Model.x0_ctrl = x0;
    Model.baseline_input_score = t0;
    Model.baseline_output_score = u0;
end
