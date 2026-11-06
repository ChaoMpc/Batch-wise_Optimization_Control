function Params = tune_hyperparameters(Ident_data, Val_data, Params)
% =========================================================================
% tune_hyperparameters -- 基于独立验证集 free-run NRMSE 严格选择超参数
%
% 判据: 在独立验证集上计算 free-run (无状态重置, 长时仿真) NRMSE,
%       取 y1, y2 两个输出的平均 NRMSE 作为选择判据。
% 输出: 控制台对齐表 (对齐截图 Table 2 的列布局)
%       以及可直接粘贴到论文的 LaTeX 代码
% =========================================================================
    fprintf('\n=== Hyperparameter Tuning (Cross-Validation) ===\n');
    fprintf('Selection criterion: Mean free-run NRMSE on validation set\n\n');

    % ---------- 搜索空间 ----------
    if isfield(Params, 'A_lv_cv_range') && ~isempty(Params.A_lv_cv_range)
        A_search = Params.A_lv_cv_range;
    else
        A_search = 1:Params.n_out;
    end
    % 避免 A > n_out 导致 Q 秩亏
    A_search = A_search(A_search <= Params.n_out);
    if isempty(A_search), A_search = Params.n_out; end

    if isfield(Params, 'n_arx_cv_range') && ~isempty(Params.n_arx_cv_range)
        n_arx_search = Params.n_arx_cv_range;
    else
        n_arx_search = [4, 6, 8, 10, 12];
    end

    nA = numel(A_search);
    nN = numel(n_arx_search);
    rec_A     = zeros(nA*nN, 1);
    rec_n     = zeros(nA*nN, 1);
    rec_nrmse = nan(nA*nN, 2);    % [NRMSE_y1, NRMSE_y2]
    rec_vaf   = nan(nA*nN, 2);    % [VAF_y1,   VAF_y2]   单位 %
    rec_ok    = false(nA*nN, 1);

    best_nrmse = inf;
    best_n_arx = Params.n_arx;
    best_A     = Params.A_lv;

    orig_warn_state = warning;
    warning('off', 'all');

    % ---------- 嵌套遍历 ----------
    k = 0;
    for A_try = A_search
        for n_arx_try = n_arx_search
            k = k + 1;
            rec_A(k) = A_try;
            rec_n(k) = n_arx_try;

            Params_temp = Params;
            Params_temp.A_lv = A_try;
            Params_temp.n_arx = n_arx_try;
            Params_temp.plot_model = false;

            try
                % 抑制 build_dypls 内部的 fprintf, 避免打乱表格
                if exist('evalc', 'builtin')
                    evalc('Model_temp = build_dypls(Ident_data, Params_temp);');
                else
                    Model_temp = build_dypls(Ident_data, Params_temp);
                end

                % ---- 判据: 独立验证集 free-run NRMSE ----
                % 与 Xtune_hyperparameters.m 一致: 取第 2, 3 个返回值
                [~, nrmses, vafs] = free_run_validation(Model_temp, Val_data, Params_temp);

                rec_nrmse(k, :) = nrmses(:)';
                rec_vaf(k, :)   = vafs(:)';
                rec_ok(k)       = true;

                mean_nrmse = mean(nrmses);
                if mean_nrmse < best_nrmse
                    best_nrmse = mean_nrmse;
                    best_n_arx = n_arx_try;
                    best_A     = A_try;
                end
            catch ME
                rec_ok(k) = false;
                fprintf('   [FAILED] A=%d, n_arx=%d : %s\n', ...
                    A_try, n_arx_try, ME.message);
            end
        end
    end
    warning(orig_warn_state);

    % ================== 控制台对齐表 ==================
    fprintf('  --------------------------------------------------------------------\n');
    fprintf('   A_lv | n_arx |   NRMSE (y1, y2)      |   VAF (y1, y2) %%      | Mean NRMSE\n');
    fprintf('  --------------------------------------------------------------------\n');
    for k = 1:nA*nN
        if rec_ok(k)
            fprintf('    %2d  |   %2d  | %7.4f, %7.4f     | %7.2f, %7.2f      | %8.4f\n', ...
                rec_A(k), rec_n(k), ...
                rec_nrmse(k,1), rec_nrmse(k,2), ...
                rec_vaf(k,1),   rec_vaf(k,2), ...
                mean(rec_nrmse(k, :)));
        else
            fprintf('    %2d  |   %2d  |   [FAILED]\n', rec_A(k), rec_n(k));
        end
    end
    fprintf('  --------------------------------------------------------------------\n');
    fprintf('  Selected best: A_lv = %d, n_arx = %d (Mean NRMSE = %.4f)\n\n', ...
        best_A, best_n_arx, best_nrmse);

    % ================== LaTeX 表格 (与截图 Table 2 一致) ==================
    fprintf('=== LaTeX Table (paste into manuscript) ===\n\n');
    fprintf('\\begin{table}[htbp]\n');
    fprintf('    \\caption{Cross-validation results for DyPLS-ARX model selection.}\n');
    fprintf('    \\label{tab:cross_validation}\n');
    fprintf('    \\centering\n');
    fprintf('    \\begin{tabular}{ccccc}\n');
    fprintf('        \\toprule\n');
    fprintf('        $A_{lv}$ & $n_{arx}$ & NRMSE ($y_1, y_2$) & VAF ($y_1, y_2$) \\%% & Mean NRMSE \\\\\n');
    fprintf('        \\midrule\n');
    for k = 1:nA*nN
        if ~rec_ok(k), continue; end
        mean_nrmse = mean(rec_nrmse(k, :));
        if rec_A(k) == best_A && rec_n(k) == best_n_arx
            fprintf('        \\textbf{%d} & \\textbf{%d} & \\textbf{%.4f, %.4f} & \\textbf{%.2f, %.2f} & \\textbf{%.4f} \\\\\n', ...
                rec_A(k), rec_n(k), rec_nrmse(k,1), rec_nrmse(k,2), ...
                rec_vaf(k,1), rec_vaf(k,2), mean_nrmse);
        else
            fprintf('        %d & %d & %.4f, %.4f & %.2f, %.2f & %.4f \\\\\n', ...
                rec_A(k), rec_n(k), rec_nrmse(k,1), rec_nrmse(k,2), ...
                rec_vaf(k,1), rec_vaf(k,2), mean_nrmse);
        end
    end
    fprintf('        \\bottomrule\n');
    fprintf('    \\end{tabular}\n');
    fprintf('\\end{table}\n\n');

    % ---------- 写回最优参数 ----------
    Params.A_lv  = best_A;
    Params.n_arx = best_n_arx;
end
