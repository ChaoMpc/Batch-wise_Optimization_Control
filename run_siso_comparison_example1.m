% =========================================================================
% run_siso_comparison_example1.m
%
% 严格对比 MIMO-ARX 与 SISO-ARX 在 Example 1 新数据上的性能。
% 注意：请确保同目录下存在 build_strict_siso.m 文件。
% =========================================================================

clear; clc; close all;
if exist('OCTAVE_VERSION', 'builtin'), pkg load control optim; end

fprintf('=== Strict SISO vs MIMO Comparison for Example 1 ===\n');

% 1. 初始化参数
Params = init_params();
Params.A_lv = 2;
if exist('rng', 'builtin'), rng(Params.rng_seed);
else, rand('seed', Params.rng_seed); randn('seed', Params.rng_seed); end

% 2. 生成新配置下的数据
fprintf('Generating identification and validation data...\n');
[Ident_data, Params] = gen_prbs_excitation(Params);
[Val_data, Params]   = gen_validation_data(Params, 'val', Ident_data);

% 3. 遍历 n_arx
n_arx_list = [2, 4, 6, 8, 10];
n_cases = length(n_arx_list);
res = struct('n_arx', [], 'siso_nrmse', [], 'siso_vaf', [], 'mimo_nrmse', [], 'mimo_vaf', []);

fprintf('\nEvaluating models (this may take a few minutes, note that n_arx=2 performs poorly due to the 7-step delay)...\n');
for i = 1:n_cases
    n = n_arx_list(i);
    fprintf('  Running n_arx = %d...\n', n);

    Params_temp = Params;
    Params_temp.n_arx = n;

    % --- 1. MIMO-ARX (严格调用原函数) ---
    Model_mimo = build_dypls(Ident_data, Params_temp);
    [~, nrmse_mimo, vaf_mimo] = free_run_validation(Model_mimo, Val_data, Params_temp);

    % --- 2. SISO-ARX (独立文件调用) ---
    Model_siso = build_strict_siso(Ident_data, Params_temp);
    [~, nrmse_siso, vaf_siso] = free_run_validation(Model_siso, Val_data, Params_temp);

    % 记录数据
    res(i).n_arx = n;
    res(i).siso_nrmse = nrmse_siso;
    res(i).siso_vaf   = vaf_siso;
    res(i).mimo_nrmse = nrmse_mimo;
    res(i).mimo_vaf   = vaf_mimo;
end

% 4. 打印表格
fprintf('\n\n');
fprintf('========================================================================================\n');
fprintf('  Comparison of SISO-ARX and MIMO-ARX on Example 1 (New Data, A_lv = 2)\n');
fprintf('========================================================================================\n');
fprintf(' %-6s | %-22s | %-22s \n', 'n_arx', 'SISO-ARX (NRMSE | VAF%)', 'MIMO-ARX (NRMSE | VAF%)');
fprintf('        |   y1          y2        |   y1          y2        \n');
fprintf('----------------------------------------------------------------------------------------\n');
for i = 1:n_cases
    fprintf('  %2d    | %.4f / %5.2f%%   %.4f / %5.2f%% | %.4f / %5.2f%%   %.4f / %5.2f%%\n', ...
        res(i).n_arx, ...
        res(i).siso_nrmse(1), res(i).siso_vaf(1), res(i).siso_nrmse(2), res(i).siso_vaf(2), ...
        res(i).mimo_nrmse(1), res(i).mimo_vaf(1), res(i).mimo_nrmse(2), res(i).mimo_vaf(2));
end
fprintf('========================================================================================\n\n');

% 5. LaTeX 表格输出
fprintf('=== LaTeX Table ===\n\n');
fprintf('\\begin{table}[htbp]\n');
fprintf('    \\caption{Comparison of SISO-ARX and MIMO-ARX under identical new data and hyperparameters ($A_{\\mathrm{lv}} = 2$).}\n');
fprintf('    \\label{tab:siso_mimo_comparison}\n');
fprintf('    \\centering\n');
fprintf('    \\begin{tabular}{c|cc|cc}\n');
fprintf('        \\toprule\n');
fprintf('        \\multirow{2}{*}{$n_{arx}$} & \\multicolumn{2}{c|}{SISO-ARX} & \\multicolumn{2}{c}{MIMO-ARX (Proposed)} \\\\\n');
fprintf('        \\cmidrule(lr){2-3} \\cmidrule(lr){4-5}\n');
fprintf('        & NRMSE ($y_1, y_2$) & VAF (\\%%) & NRMSE ($y_1, y_2$) & VAF (\\%%) \\\\\n');
fprintf('        \\midrule\n');
for i = 1:n_cases
    if res(i).n_arx == 8
        fprintf('        \\textbf{%d} & \\textbf{%.4f, %.4f} & \\textbf{%.2f, %.2f} & \\textbf{%.4f, %.4f} & \\textbf{%.2f, %.2f} \\\\\n', ...
            res(i).n_arx, ...
            res(i).siso_nrmse(1), res(i).siso_nrmse(2), res(i).siso_vaf(1), res(i).siso_vaf(2), ...
            res(i).mimo_nrmse(1), res(i).mimo_nrmse(2), res(i).mimo_vaf(1), res(i).mimo_vaf(2));
    else
        fprintf('        %d & %.4f, %.4f & %.2f, %.2f & %.4f, %.4f & %.2f, %.2f \\\\\n', ...
            res(i).n_arx, ...
            res(i).siso_nrmse(1), res(i).siso_nrmse(2), res(i).siso_vaf(1), res(i).siso_vaf(2), ...
            res(i).mimo_nrmse(1), res(i).mimo_nrmse(2), res(i).mimo_vaf(1), res(i).mimo_vaf(2));
    end
end
fprintf('        \\bottomrule\n');
fprintf('    \\end{tabular}\n');
fprintf('\\end{table}\n\n');
