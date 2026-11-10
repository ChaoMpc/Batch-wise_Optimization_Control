% =========================================================================
% run_sensitivity_A.m
% 扫频分析并格式化输出数据表格（彻底修复可行性 Bug 并增强展示）
% =========================================================================

clear; clc; close all;

try
    pkg load control;
catch
    warning('无法加载 control 包，可能会导致后续 dlqr 计算失败。');
end

% 1. 加载基础配置
load('results/stage01_setup.mat', 'Params', 'Reference');

% === 恢复 300 批训练数据 ===
Params.N_train_batches = 300;
Params.N_val_batches = 5;
Params.N_unc_batches = 5;

% 2. 生成一次数据（复用）
fprintf('Generating data for sensitivity analysis (300 batches)...\n');
Train_data = generate_experiment_data(Params, Params.N_train_batches, Params.rng_seed, 'train');
Val_data   = generate_experiment_data(Params, Params.N_val_batches, Params.rng_seed_val, 'val');
Cal_data   = generate_experiment_data(Params, Params.N_unc_batches, Params.rng_seed_unc, 'cal');

% 3. 扫描组合
A_list = [1, 2, 3, 4, 5];
n_list = [2, 3, 4];

results = [];

fprintf('\nStarting sensitivity scan. This may take a few minutes...\n');

for A = A_list
    for n = n_list
        fprintf('Running A_lv=%d, n_arx=%d...\n', A, n);

        Params_iter = Params;
        Params_iter.A_lv = A;
        Params_iter.n_arx = n;

        res = struct('A_lv', A, 'n_arx', n, ...
                     'vaf', NaN, 'nrmse', NaN, ...
                     'c_atr_lo', NaN, 'c_atr_hi', NaN, ...
                     'reach_err', NaN, 'feasible', false);

        try
            % 预处理
            Prep = fit_preprocessing(Train_data, Params_iter);
            Std_train = apply_preprocessing(Train_data, Prep, Params_iter);
            Std_val   = apply_preprocessing(Val_data, Prep, Params_iter);
            Std_cal   = apply_preprocessing(Cal_data, Prep, Params_iter);

            % 辨识模型
            Model = build_identified_model(Prep, Std_train, A, n, Params_iter);
            Params_iter.A_lv = Model.A_lv;
            Params_iter.n_arx = Model.n_arx;
            res.A_lv = Model.A_lv;

            % 验证
            Val = validate_model(Model, Std_val, Val_data, Params_iter);
            res.nrmse = mean(Val.NRMSE);
            res.vaf   = mean(Val.VAF);

            % 估算不确定性
            Cal_resid = collect_model_residuals(Cal_data, Model, Prep, Params_iter);
            UncSets = estimate_uncertainty(Cal_resid, Cal_resid, Cal_resid, Model, Params_iter);

            % 计算紧化约束
            StateError = compute_state_error_set(Model, Params_iter);
            Feedback = design_tube_feedback(Model, Params_iter);
            Invariant = compute_invariant_ellipsoid(Model, Feedback, UncSets, StateError, Params_iter);
            Tube = compute_tightening(Model, Feedback, Invariant, UncSets, StateError, Params_iter);

            % 提取 C_ATR 的紧化区间
            trk = Params_iter.track_idx; % C_ATR 是第4个输出
            res.c_atr_lo = Params_iter.y_phys_min(trk) + Tube.delta_y(Params_iter.n_out + trk);
            res.c_atr_hi = Params_iter.y_phys_max(trk) - Tube.delta_y(trk);

            % 检查是否可行
            Feasibility = check_tightened_feasibility(Model, Tube, Params_iter);

            if Feasibility.feasible
                % === 核心修复：补充这一行，否则 compute_initial_feasible 会报错 ===
                Tube.nominal_feasible = true;
                Tube.feasible_seed = Feasibility.t_feasible;

                Problem = build_nominal_problem(Model, Tube, Params_iter);
                Initial = compute_initial_feasible(Problem, Tube, Params_iter);
                Optimal = compute_optimal_reachable(Problem, Initial.t, Params_iter);

                reach_err = abs(Optimal.selected_outputs - Params_iter.y_ref_phys);
                res.reach_err = reach_err(2); % 记录 120 min 处的误差
                res.feasible = true;
            else
                res.feasible = false;
            end

        catch ME
            fprintf('  Error: %s\n', ME.message);
        end

        results = [results; res];
    end
end

% =========================================================================
% 4. 打印极度直观的表格
% =========================================================================
fprintf('\n\n');
fprintf('========================================================================================================================\n');
fprintf('%-6s | %-6s | %-10s | %-10s | %-22s | %-12s | %-10s\n', ...
    'A_lv', 'n_arx', 'Out VAF(%)', 'Out NRMSE', 'C_ATR Tightened Interval', 'Reach Err(120)', 'Status');
fprintf('------------------------------------------------------------------------------------------------------------------------\n');

for i = 1:numel(results)
    r = results(i);

    if isnan(r.c_atr_lo)
        interval_str = 'N/A';
    else
        interval_str = sprintf('[%.5f, %.5f]', r.c_atr_lo, r.c_atr_hi);
    end

    if r.feasible
        status_str = 'Feasible';
        err_str = sprintf('%12.2e', r.reach_err);
    else
        % 如果不可行，判断是否是因为约束太紧导致参考点被排除
        if r.c_atr_lo > 0.015
            status_str = 'INFEASIBLE (0.015 < Lower Bound)';
        else
            status_str = 'INFEASIBLE (LP failed)';
        end
        err_str = '         N/A';
    end

    fprintf('%-6d | %-6d | %-10.2f | %-10.4f | %-22s | %s | %-10s\n', ...
        r.A_lv, r.n_arx, r.vaf, r.nrmse, interval_str, err_str, status_str);
end
fprintf('========================================================================================================================\n');

% 5. 保存为 CSV
csv_file = 'sensitivity_A_data.csv';
fid = fopen(csv_file, 'w');
fprintf(fid, 'A_lv,n_arx,Out_VAF,Out_NRMSE,C_ATR_Lower,C_ATR_Upper,Reach_Err_120,Feasible\n');
for i = 1:numel(results)
    r = results(i);
    fprintf(fid, '%d,%d,%.4f,%.4f,%.5f,%.5f,%.4e,%d\n', ...
        r.A_lv, r.n_arx, r.vaf, r.nrmse, r.c_atr_lo, r.c_atr_hi, r.reach_err, r.feasible);
end
fclose(fid);
fprintf('\nData saved to: %s\n', csv_file);
