function T = report_pls_evidence_example2(Model, Prep, Validation, Val_data, Params)
% 输出三张表：逐成分方差、Q 奇异值、输入方差贡献，写入 CSV/TXT。

    A = Model.A_lv;
    PLS = Model.PLS;

    %% 表 A：逐成分方差（PLS 累进）与 Q 奇异值
    Qmat = PLS.Q;
    sv_Q = svd(Qmat, 'econ');
    sQmin = sv_Q(end);
    condQ = sv_Q(1)/sv_Q(end);

    fprintf('\n========= Table A. PLS evidence per component =========\n');
    fprintf('%3s %10s %10s %14s %10s\n', 'a', 'EV_X(%)', 'EV_Y(%)', 'sigma_min(Q_a)', 'rank(Q_a)');
    for a = 1:A
        Q_a = Qmat(:,1:a);
        fprintf('%3d %10.3f %10.3f %14.4e %10d\n', ...
            a, 100*PLS.EV_X(a), 100*PLS.EV_Y(a), svd(Q_a,'econ')(end), rank(Q_a));
    end
    fprintf('\nTotal sigma_Q = %.4e, cond(Q) = %.4e\n', sQmin, condQ);

    %% 表 B：无噪声输出空间有效秩
    Y_true_std_all = [];
    for b = 1:Val_data.N_batches
        Y_t = Val_data.Y_true(:, 1:Params.Ble, b);
        Y_dev = Y_t - Params.y_base(:, 1:Params.Ble);
        Y_std_b = Y_dev ./ Model.Y_std;

        Y_true_std_all = [Y_true_std_all, Y_std_b]; %#ok<AGROW>
    end

    s_true = svd(Y_true_std_all, 'econ');
    true_rank = sum(s_true > 1e-8 * s_true(1));

    fprintf('\n========= Table B. Effective output rank =========\n');
    fprintf('Top-6 singular values of Y_true (standardized): \n');
    fprintf('  %.3e ', s_true(1:min(6,numel(s_true))));
    fprintf('\n');
    fprintf('Effective rank at 1e-8 threshold: %d\n', true_rank);
    fprintf('Latent A_lv: %d\n', A);

    %% 表 C：输入方差贡献
    fprintf('\n========= Table C. Input reconstruction =========\n');
    fprintf('Input channel         VAF(%%)    NRMSE\n');
    for i = 1:Params.n_in_opt
        fprintf('%-18s  %7.2f  %8.4f\n', Params.input_names{i}, ...
            Validation.input_VAF(i), Validation.input_NRMSE(i));
    end
    fprintf('Cumulative EV_X = %.2f%%\n', 100*PLS.EV_X(A));
end
