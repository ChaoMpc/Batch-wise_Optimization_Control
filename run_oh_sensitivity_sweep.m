% =========================================================================
% run_oh_sensitivity_sweep.m
% 针对审稿人意见 3b，对 Oh 方法的 Q_trk 和 R_upd 进行 3x3 局部扫描
% =========================================================================
clear; clc; close all;
if exist('OCTAVE_VERSION','builtin'), pkg load control optim; end
addpath(pwd); rehash();

% 1. 加载所需数据
res_dir = fullfile(pwd, 'results');
S4 = load(fullfile(res_dir,'stage04_nominal.mat'));
S3 = load(fullfile(res_dir,'stage03_tube.mat'));
S5 = load(fullfile(res_dir,'stage05_simulation.mat'));

Params = S4.Params; Model = S4.Model; Tube = S4.Tube;
Problem = S4.Problem; Initial = S4.Initial; Optimal = S4.Optimal;
UncSets = S3.UncSets; Scenario = S5.Scenario;

% 2. 定义扫描网格
Q_list = [10, 100, 1000];
R_list = [1e-5, 1e-4, 1e-3];
nQ = length(Q_list); nR = length(R_list);

Params.N_batches = 5;
Params.plot_results = false;

% 3. 扩展存储矩阵
res_err    = nan(nR, nQ);
res_hits   = nan(nR, nQ);
res_upd    = nan(nR, nQ);
res_tv     = nan(nR, nQ);
res_gap    = nan(nR, nQ);
res_dY     = nan(nR, nQ);
res_dT     = nan(nR, nQ);
res_qp     = nan(nR, nQ);

res_err_batch = nan(nR, nQ, Params.N_batches);
Data_oh_all = cell(nR, nQ);

fprintf('=== 开始 Oh 方法参数敏感度扫描 (共 %d 组组合) ===\n', nR * nQ);
total_timer = tic();

for i = 1:nR
    for j = 1:nQ
        fprintf('\n======================================================================\n');
        fprintf('>>> [进度 %d/%d] 正在运行: Q_trk=%4d, R_upd=%.0e\n', ...
            (i-1)*nQ + j, nR*nQ, Q_list(j), R_list(i));
        fprintf('======================================================================\n');

        Params_try = Params;
        Params_try.Q_weight = Q_list(j);
        Params_try.R_weight = R_list(i);

        timer_inner = tic;

        % 使用 try-catch 防止某组极端参数导致 QP 求解器崩溃中断整个实验
        try
            Data_oh = run_oh_ilmpc(Problem, Initial, Optimal, Scenario, Model, Tube, UncSets, Params_try);

            % 安全存储完整轨迹数据
            Data_oh_all{i, j} = Data_oh;

            if isfield(Data_oh, 'QP_times') && ~isempty(Data_oh.QP_times)
                % 计算每个批次的平均耗时: sum(QP_times, 2) 算出每批耗时，再求均值
                res_qp(i, j) = mean(sum(Data_oh.QP_times, 2));
            end

            if exist('evalc', 'builtin')
                evalc('Metrics_oh = compute_metrics(Data_oh, Params_try);');
            else
                Metrics_oh = compute_metrics(Data_oh, Params_try);
            end

            % --- 提取基础指标 ---
            res_err(i, j)  = Metrics_oh.final_actual_reference_error;
            res_hits(i, j) = Metrics_oh.input_boundary_hits;

            % 提取每个批次的参考点跟踪误差
            if isfield(Data_oh, 'Err_true_act')
                res_err_batch(i, j, :) = Data_oh.Err_true_act;
            end

            % --- 提取控制行为指标 ---
            B_actual = size(Data_oh.V_hist, 3);
            if isfield(Data_oh, 'input_change_norm') && ~isempty(Data_oh.input_change_norm)
                res_upd(i, j) = max(Data_oh.input_change_norm(:));
            end

            tv_batch = zeros(1, B_actual);
            for b = 1:B_actual
                dV_t = diff(squeeze(Data_oh.V_hist(:,:,b)), 1, 2);
                tv_batch(b) = sum(abs(dV_t(:)));
            end
            res_tv(i, j) = mean(tv_batch);

            if isfield(Data_oh, 'V_request_hist') && ~isempty(Data_oh.V_request_hist)
                v_req = Data_oh.V_request_hist(:);
                v_app = Data_oh.V_hist(:);
                res_gap(i, j) = max(abs(v_req - v_app));
            end

            % --- 提取批次发散度指标 ---
            dY_vec = reshape(Data_oh.Y_hist(:,:,B_actual) - Data_oh.Y_hist(:,:,1), [], 1);
            res_dY(i, j) = norm(dY_vec, 2);

            dT_vec = reshape(Data_oh.T_hist(:,:,B_actual) - Data_oh.T_hist(:,:,1), [], 1);
            res_dT(i, j) = norm(dT_vec, 2);

            elapsed = toc(timer_inner);
            fprintf('>>> 本组完成! 耗时: %.1f s | 最终误差: %.4e | 撞线: %d次\n', elapsed, res_err(i, j), res_hits(i, j));
            fprintf('    时域指标: TV=%.2f | MaxUpd=%.2f | Gap=%.2f\n', res_tv(i,j), res_upd(i,j), res_gap(i,j));
            fprintf('    批次指标: ||Y_5-Y_1||=%.2e | ||T_5-T_1||=%.2e\n', res_dY(i,j), res_dT(i,j));

        catch ME
            fprintf(2, '>>> [警告] 该组参数求解失败，已跳过。错误信息: %s\n', ME.message);
        end
    end
end

fprintf('\n\n=== 所有扫描已完成，总耗时: %.1f 分钟 ===\n', toc(total_timer)/60);

% 5. 打印结果汇总表
fprintf('\n=== 扫描结果汇总 ===\n');
for i = 1:nR
    for j = 1:nQ
        fprintf('Q=%4d, R=%.0e : Err=%.4e | Hits=%3d | TV=%.1f | dY=%.1e\n', ...
            Q_list(j), R_list(i), res_err(i, j), res_hits(i, j), res_tv(i,j), res_dY(i,j));
    end
end

% 6. 保存到 .mat
save('-mat7-binary', 'oh_sensitivity_workspace.mat', 'Q_list', 'R_list', ...
    'res_err', 'res_hits', 'res_upd', 'res_tv', 'res_gap', 'res_dY', 'res_dT', 'res_qp', ...
    'res_err_batch', 'Data_oh_all');

% 7. 导出 TXT
fid = fopen('oh_sensitivity_results.txt', 'w');
if fid > 0
    fprintf(fid, '=== Oh IL-MPC Parameter Sensitivity Analysis Results ===\n\n');

    fprintf(fid, '--- 1. Raw Matrices ---\n');
    fprintf(fid, 'Q_list = [%s]\n', num2str(Q_list));
    fprintf(fid, 'R_list = [%s]\n\n', num2str(R_list));
    fprintf(fid, 'res_err (Final Error) =\n'); disp(res_err);
    fprintf(fid, 'res_hits (Constraint Hits) =\n'); disp(res_hits);
    fprintf(fid, 'res_tv (Total Variation of Input) =\n'); disp(res_tv);
    fprintf(fid, 'res_dY (Batch trajectory divergence ||Y_5 - Y_1||) =\n'); disp(res_dY);

    fprintf(fid, '\n--- 2. Batch Tracking Error Evolution ---\n');
    for i = 1:nR
        for j = 1:nQ
            fprintf(fid, 'Q=%4d, R=%.0e : ', Q_list(j), R_list(i));
            fprintf(fid, '%.4e  ', squeeze(res_err_batch(i,j,:)));
            fprintf(fid, '\n');
        end
    end

    fprintf(fid, '\n--- 3. Human-Readable Table ---\n');
    fprintf(fid, '%-8s | %-6s | %-12s | %-6s | %-8s | %-8s | %-8s | %-10s | %-10s\n', ...
        'R_upd', 'Q_trk', 'Error', 'Hits', 'TV', 'Gap', 'QP_Time', 'dY_Norm', 'dT_Norm');
    fprintf(fid, repmat('-', 1, 105)); fprintf(fid, '\n');
    for i = 1:nR
        for j = 1:nQ
            fprintf(fid, '%-8.0e | %-6d | %-12.4e | %-6d | %-8.2f | %-8.2f | %-8.4f | %-10.4e | %-10.4e\n', ...
                R_list(i), Q_list(j), res_err(i,j), res_hits(i,j), ...
                res_tv(i,j), res_gap(i,j), res_qp(i,j), res_dY(i,j), res_dT(i,j));
        end
    end
    fclose(fid);
    fprintf('数据已成功保存到 oh_sensitivity_results.txt 和 oh_sensitivity_workspace.mat 中。\n');
end
