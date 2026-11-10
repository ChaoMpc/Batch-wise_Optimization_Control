function results = evaluate_on_test_set(Model, Params, Ident_data, n_batches, N_test_override, seed_override)
% =========================================================================
% evaluate_on_test_set -- Evaluate the trained DyPLS-ARX model on the
%                         INDEPENDENT TEST SET, batch by batch.
%
% 输入:
%   Model           - build_dypls 输出的最终模型 (用 Ident_data 训练)
%   Params          - 参数结构体 (init_params 的输出或已调优版本)
%   Ident_data      - 训练集数据 (提供训练统计量, 用于测试集标准化)
%   n_batches       - (可选) 切分的 batch 数目, 默认 5
%   N_test_override - (可选) 测试集总长度, 默认 Params.N_test
%   seed_override   - (可选) 测试集随机种子, 默认 Params.rng_seed_test
%
% 输出:
%   results         - 结构体, 含每个 batch 的 NRMSE / VAF / 切分索引等
%
% 说明:
%   1) 测试集用与识别集 (rng_seed)、验证集 (rng_seed_val) 都不同的
%      rng_seed_test 生成, 且用训练集统计量重新标准化, 属于 untouched
%      test set, 不参与任何超参数选择。
%   2) 每个 batch 单独调用 free_run_validation, 即每个 batch 从
%      前 n_arx 步真实值初始化后 free-run, 得到该 batch 的 NRMSE/VAF。
%   3) 控制台表格和 LaTeX 表格的第一列为 Batch, 其余列与截图一致:
%      Batch | NRMSE (y1, y2) | VAF (y1, y2) % | Mean NRMSE
%
% 典型调用:
%   results = evaluate_on_test_set(Model, Params, Ident_data, 5);
% =========================================================================

    if nargin < 4 || isempty(n_batches),        n_batches       = 5;                  end
    if nargin < 5 || isempty(N_test_override),  N_test_override = Params.N_test;      end
    if nargin < 6 || isempty(seed_override),    seed_override   = Params.rng_seed_test; end

    fprintf('\n=== Independent Test Set Evaluation (batch-wise) ===\n');
    fprintf('  n_batches = %d, N_test = %d, seed = %d\n', ...
        n_batches, N_test_override, seed_override);

    % ---------- 1. 生成独立测试集 ----------
    Params_test = Params;
    Params_test.N_test        = N_test_override;
    Params_test.rng_seed_test = seed_override;

    [Test_data, ~] = gen_validation_data(Params_test, 'test', Ident_data);
    N_total = size(Test_data.X0, 1);
    n_out   = Params.n_out;
    n_arx   = Model.n_arx;

    if N_total < n_batches * (n_arx + 2)
        error(['evaluate_on_test_set: 测试集 N=%d 太短, 不能分成 %d 个 batch ' ...
               '(每个 batch 至少需要 n_arx+2 = %d 个样本)。'], ...
               N_total, n_batches, n_arx + 2);
    end

    % ---------- 2. 切分 (余数均匀分配到前几个 batch) ----------
    N_batch     = floor(N_total / n_batches);
    rem_samples = N_total - N_batch * n_batches;
    batch_idx   = cell(n_batches, 1);
    start       = 1;
    for b = 1:n_batches
        len_b = N_batch + (b <= rem_samples);
        batch_idx{b} = start : (start + len_b - 1);
        start        = start + len_b;
    end

    % ---------- 3. 逐 batch free-run ----------
    rec_nrmse = nan(n_batches, n_out);
    rec_vaf   = nan(n_batches, n_out);
    rec_ok    = false(n_batches, 1);
    Y_sim_all = cell(n_batches, 1);
    Y_true_all = cell(n_batches, 1);

    orig_warn_state = warning;
    warning('off', 'all');
    for b = 1:n_batches
        idx = batch_idx{b};
        Sub_data = subset_validation_data(Test_data, idx);
        try
            % 抑制 free_run_validation 内部的 fprintf, 避免打乱表格
            if exist('evalc', 'builtin')
                evalc('[Y0_sim_b, nrmses_b, vafs_b] = free_run_validation(Model, Sub_data, Params);');
            else
                [Y0_sim_b, nrmses_b, vafs_b] = free_run_validation(Model, Sub_data, Params);
            end
            rec_nrmse(b, :) = nrmses_b(:)';
            rec_vaf(b, :)   = vafs_b(:)';
            Y_sim_all{b}    = Y0_sim_b;
            Y_true_all{b}   = Sub_data.Y0;
            rec_ok(b)       = true;
        catch ME
            fprintf('   [FAILED] Batch %d : %s\n', b, ME.message);
        end
    end
    warning(orig_warn_state);

    % ---------- 4. 控制台对齐表 ----------
    fprintf('\n  --------------------------------------------------------------------\n');
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
    ok_mask = rec_ok;
    if any(ok_mask)
        fprintf('   ------+------------------------+------------------------+----------\n');
        fprintf('   Mean  | %7.4f, %7.4f     | %7.2f, %7.2f      | %8.4f\n', ...
            mean(rec_nrmse(ok_mask,1)), mean(rec_nrmse(ok_mask,2)), ...
            mean(rec_vaf(ok_mask,1)),   mean(rec_vaf(ok_mask,2)), ...
            mean(mean(rec_nrmse(ok_mask, :))));
    end
    fprintf('  --------------------------------------------------------------------\n\n');

    % ---------- 5. LaTeX 表格 ----------
    fprintf('=== LaTeX Table (paste into manuscript) ===\n\n');
    fprintf('\\begin{table}[htbp]\n');
    fprintf('    \\caption{Free-run validation on the independent test set (batch-wise).}\n');
    fprintf('    \\label{tab:test_batches}\n');
    fprintf('    \\centering\n');
    fprintf('    \\begin{tabular}{ccccc}\n');
    fprintf('        \\toprule\n');
    fprintf('        Batch & NRMSE ($y_1, y_2$) & VAF ($y_1, y_2$) \\%% & Mean NRMSE \\\\\n');
    fprintf('        \\midrule\n');
    for b = 1:n_batches
        if ~rec_ok(b), continue; end
        fprintf('        %d & %.4f, %.4f & %.2f, %.2f & %.4f \\\\\n', ...
            b, rec_nrmse(b,1), rec_nrmse(b,2), ...
            rec_vaf(b,1),   rec_vaf(b,2), ...
            mean(rec_nrmse(b, :)));
    end
    if any(ok_mask)
        fprintf('        \\midrule\n');
        fprintf('        Mean & %.4f, %.4f & %.2f, %.2f & %.4f \\\\\n', ...
            mean(rec_nrmse(ok_mask,1)), mean(rec_nrmse(ok_mask,2)), ...
            mean(rec_vaf(ok_mask,1)),   mean(rec_vaf(ok_mask,2)), ...
            mean(mean(rec_nrmse(ok_mask, :))));
    end
    fprintf('        \\bottomrule\n');
    fprintf('    \\end{tabular}\n');
    fprintf('\\end{table}\n\n');

    % ---------- 6. 打包结果 ----------
    results.Test_data  = Test_data;
    results.batch_idx  = batch_idx;
    results.nrmse      = rec_nrmse;
    results.vaf        = rec_vaf;
    results.ok         = rec_ok;
    results.Y0_sim_all = Y_sim_all;
    results.Y0_true_all = Y_true_all;
    results.N_total    = N_total;
    results.n_batches  = n_batches;
    results.n_arx      = n_arx;
end

% ========================================================================
function Sub = subset_validation_data(Data, idx)
% 按行索引截取 Val_data/Test_data, 保持所有字段一致
    Sub = Data;
    row_fields = {'X','Y','X_abs','Y_abs','X0','Y0','time'};
    N_ref = size(Data.X, 1);
    for f = 1:numel(row_fields)
        fname = row_fields{f};
        if isfield(Data, fname) && size(Data.(fname), 1) == N_ref
            Sub.(fname) = Data.(fname)(idx, :);
        end
    end
end
