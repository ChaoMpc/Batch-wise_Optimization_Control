function results = evaluate_test_batches(Model, Test_batches, Params)
% =========================================================================
% evaluate_test_batches -- 在独立测试集的多个 batch 上评估模型
%
% 输入:
%   Model        - 训练好的 DyPLS-ARX 模型 (build_dypls 输出)
%   Test_batches - cell 数组, 每个元素是一个 batch 的数据结构体,
%                  需包含与 Val_data 相同的字段 (如 .u, .y)
%   Params       - 参数结构体
%
% 输出:
%   results      - 结构体数组, 包含每个 batch 的 NRMSE 和 VAF
%
% 功能:
%   对每个 batch 调用 free_run_validation, 计算 NRMSE 和 VAF,
%   并打印控制台对齐表格及可直接粘贴的 LaTeX 表格。
% =========================================================================

    n_batches = numel(Test_batches);
    n_out = Params.n_out;   % 输出维数 (通常为 2)

    rec_nrmse = nan(n_batches, n_out);
    rec_vaf   = nan(n_batches, n_out);
    rec_ok    = false(n_batches, 1);

    fprintf('\n=== Independent Test Set Evaluation (Batch-wise) ===\n');
    fprintf('Criterion: free-run NRMSE and VAF on each test batch\n\n');

    % 抑制可能出现的警告
    orig_warn_state = warning;
    warning('off', 'all');

    for b = 1:n_batches
        try
            % 对单个 batch 进行 free-run 验证
            [~, nrmses, vafs] = free_run_validation(Model, Test_batches{b}, Params);
            rec_nrmse(b, :) = nrmses(:)';
            rec_vaf(b, :)   = vafs(:)';
            rec_ok(b)       = true;
        catch ME
            rec_ok(b) = false;
            fprintf('   [FAILED] Batch %d : %s\n', b, ME.message);
        end
    end
    warning(orig_warn_state);

    % ================== 控制台对齐表 ==================
    fprintf('  --------------------------------------------------------------------\n');
    fprintf('   Batch |   NRMSE (y1, y2)      |   VAF (y1, y2) %%      | Mean NRMSE\n');
    fprintf('  --------------------------------------------------------------------\n');
    for b = 1:n_batches
        if rec_ok(b)
            fprintf('    %2d   | %7.4f, %7.4f     | %7.2f, %7.2f      | %8.4f\n', ...
                b, rec_nrmse(b,1), rec_nrmse(b,2), ...
                rec_vaf(b,1),   rec_vaf(b,2), ...
                mean(rec_nrmse(b, :)));
        else
            fprintf('    %2d   |   [FAILED]\n', b);
        end
    end
    fprintf('  --------------------------------------------------------------------\n\n');

    % ================== LaTeX 表格 ==================
    fprintf('=== LaTeX Table (paste into manuscript) ===\n\n');
    fprintf('\\begin{table}[htbp]\n');
    fprintf('    \\caption{Free-run validation results on the independent test set for batches 1--5.}\n');
    fprintf('    \\label{tab:test_batches}\n');
    fprintf('    \\centering\n');
    fprintf('    \\begin{tabular}{ccccc}\n');
    fprintf('        \\toprule\n');
    fprintf('        Batch & NRMSE ($y_1, y_2$) & VAF ($y_1, y_2$) \\%% & Mean NRMSE \\\\\n');
    fprintf('        \\midrule\n');
    for b = 1:n_batches
        if ~rec_ok(b), continue; end
        mean_nrmse = mean(rec_nrmse(b, :));
        fprintf('        %d & %.4f, %.4f & %.2f, %.2f & %.4f \\\\\n', ...
            b, rec_nrmse(b,1), rec_nrmse(b,2), ...
            rec_vaf(b,1),   rec_vaf(b,2), mean_nrmse);
    end
    fprintf('        \\bottomrule\n');
    fprintf('    \\end{tabular}\n');
    fprintf('\\end{table}\n\n');

    % 整理输出
    results = struct('batch', num2cell(1:n_batches), ...
                     'nrmse', num2cell(rec_nrmse, 2), ...
                     'vaf', num2cell(rec_vaf, 2), ...
                     'ok', num2cell(rec_ok));
end
