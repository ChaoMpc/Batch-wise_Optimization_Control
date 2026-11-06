% =========================================================================
% run_fast_pipeline.m
%
% 快捷流水线：当仅修改 init_params.m 中的控制参数时，跳过低效的 Stage 2 辨识。
% 流程：跑 Stage 1 (新参数) -> 拼装旧模型 -> 跑 Stage 3, 4, 5, 6。
% =========================================================================

clear; clc; close all;
addpath(pwd); rehash;

fprintf('=== Step 1: 运行 Stage 1 获取最新基线 ===\n');
run('run_stage01_setup.m');

fprintf('\n=== Step 2: 拼装 Stage 2 (复用旧模型，注入新参数) ===\n');
S1 = load('results/stage01_setup.mat');
S2 = load('results/stage02_identification.mat');

% 1. 强制用最新的 Params 和 Reference 覆盖旧的
Params = S1.Params;
Reference = S1.Reference;

% 2. 提取旧的模型和数据集
Model       = S2.Model;
Prep        = S2.Prep;
Train_data  = S2.Train_data;
Val_data    = S2.Val_data;
Cal_data    = S2.Cal_data;
CV          = S2.CV;
Validation  = S2.Validation;
Calibration = S2.Calibration;

% 3. 检查模型维度是否一致（如果 A_lv 或 n_arx 改变，说明模型结构变了，必须重跑 Stage 2）
if Params.A_lv ~= Model.A_lv || Params.n_arx ~= Model.n_arx
    error('【致命错误】模型维度发生变化 (A_lv 或 n_arx 改动)，模型失效，必须重跑完整的 Stage 2！');
end

% 4. 补充并更新 Artifact (确保后续 Stage 的完整性校验通过)
Artifact = S2.Artifact;
Artifact.created_at = datestr(now,31);
Artifact.note = 'Artifact remixed by run_fast_pipeline using new Params.';

% 5. 存回 stage02 文件，供后续读取
save('-mat7-binary', 'results/stage02_identification.mat', ...
     'Params', 'Reference', 'Prep', 'Model', 'CV', 'Validation', 'Calibration', ...
     'Train_data', 'Val_data', 'Cal_data', 'Artifact');

fprintf('Stage 2 拼装完成 (新参数已注入)。\n');

fprintf('\n=== Step 3: 运行 Stage 3 (更新 Tube 紧化) ===\n');
run('run_stage03_tube.m');

fprintf('\n=== Step 4: 运行 Stage 4 (更新 Problem 及最优解) ===\n');
run('run_stage04_nominal.m');

fprintf('\n=== Step 5: 运行 Stage 5 (在线仿真对比) ===\n');
run('run_stage05_simulation.m');

fprintf('\n=== Step 6: 运行 Stage 6 (生成图表) ===\n');
run('run_stage06_reporting.m');

fprintf('\n============================================================\n');
fprintf('全部完成！新参数已生效，请去 results 文件夹查看日志和图表。\n');
fprintf('============================================================\n');
