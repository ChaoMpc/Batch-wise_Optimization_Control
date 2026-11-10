function run_mc_initial_states(Ident_data, Model, Tube, UncSets, Params)
% =========================================================================
% run_mc_initial_states -- MC 统计不同初始状态下的跟踪误差与约束违反
%
% 针对三类初始状态 (small, mid, large)，分别跑 N_mc 次闭环。
% 每次运行:
%   1. 改变随机种子（改变噪声）
%   2. 在固定二次型 e^T P_e^{-1} e = s^2 下，重新随机生成方向向量 z
%   3. 记录第 15 批的 3 个参考点误差，以及约束违反次数
% =========================================================================

    if exist('OCTAVE_VERSION', 'builtin'), pkg load control optim; end

    N_mc = 10; % MC 运行次数，审稿人要求提供，这里设为 10，可自行修改

    % 初始状态类别定义
    s_vals = [0.3, 0.7, 0.99];
    names  = {'Small', 'Mid', 'Large'};

    % 为存储结果预分配内存
    % 维度: [类别数, MC次数, 参考点数]
    err_all = zeros(3, N_mc, numel(Params.Kp));
    viol_v_all = zeros(3, N_mc);
    viol_y_all = zeros(3, N_mc);

    % 计算 P_e 的 Cholesky 分解，用于生成落在椭球内的扰动
    L_Pe = chol(Tube.P_e, 'lower');
    nx   = size(Tube.P_e, 1);

    fprintf('\n=== Monte Carlo Simulation for Initial States (N_mc = %d) ===\n', N_mc);

    for c = 1:3 % 遍历 3 种初始状态类别
        s_val = s_vals(c);
        fprintf('\n>>> Testing %s Perturbation (s = %.2f, s^2 = %.4f) ...\n', ...
            names{c}, s_val, s_val^2);

        for mc = 1:N_mc
            % 1. 设定随机种子 (保证每次噪声不同)
            if exist('rng', 'builtin'), rng(Params.rng_seed + 10000 * c + mc);
            else, rand('seed', Params.rng_seed + 10000 * c + mc); end

            % 2. 重新生成随机方向 z，并计算扰动 e(0)
            z_rand = randn(nx, 1);
            z_rand = z_rand / norm(z_rand);    % 归一化，使其范数为1
            e_pert = s_val * (L_Pe * z_rand);  % 此时 e_pert' * P_e^{-1} * e_pert = s_val^2

            % 3. 配置 Params 并运行闭环
            Params_mc = Params;
            Params_mc.plot_results = false;      % 关闭画图，加速运行
            Params_mc.init_state_vector = e_pert; % 注入初始扰动

            % 为了让 evaluate 不受主流程干扰，这里直接跑
            Data_mc = run_ilmpc(Ident_data, Model, Tube, UncSets, Params_mc);

            % 4. 提取最后一批（第 N_batches 批）的跟踪误差
            b_end = Params.N_batches;
            for j = 1:numel(Params.Kp)
                k_j = Params.Kp(j);
                % 实际误差 = 实际输出偏差 - 参考点偏差
                err_all(c, mc, j) = Data_mc.Y_hist(1, k_j, b_end) - Params.Ref_points(j);
            end

            % 5. 提取约束违反次数
            m_mc = compute_metrics(Data_mc, Params, Model, UncSets);
            viol_v_all(c, mc) = m_mc.n_violations_v;
            viol_y_all(c, mc) = m_mc.n_violations_y;

            fprintf('  MC %2d/%2d done. Viol(v,y)=(%d,%d)\n', mc, N_mc, viol_v_all(c,mc), viol_y_all(c,mc));
        end
    end

    % ================== 统计与打印结果 ==================
    fprintf('\n\n========================================================================================\n');
    fprintf('  Table: Monte Carlo Statistics for Final Batch (Batch %d) Tracking Errors\n', Params.N_batches);
    fprintf('========================================================================================\n');
    fprintf(' %-8s | %-4s | %-20s | %-20s | %-20s | %-10s\n', ...
        'Init', 'MC', 'Error Kp1 (Mean±Std)', 'Error Kp2 (Mean±Std)', 'Error Kp3 (Mean±Std)', 'Viol(v,y)');
    fprintf('----------------------------------------------------------------------------------------\n');

    mean_err = zeros(3, numel(Params.Kp));
    std_err  = zeros(3, numel(Params.Kp));
    total_viol_v = zeros(3,1);
    total_viol_y = zeros(3,1);

    for c = 1:3
        % 计算均值、标准差 (跨 MC 维度)
        mean_err(c, :) = mean(err_all(c, :, :), 2);
        std_err(c, :)  = std(err_all(c, :, :), 0, 2);
        % 计算约束违反总数
        total_viol_v(c) = sum(viol_v_all(c, :));
        total_viol_y(c) = sum(viol_y_all(c, :));

        fprintf(' %-8s | %-4d | %7.4f ± %.4f    | %7.4f ± %.4f    | %7.4f ± %.4f    | %d, %d\n', ...
            names{c}, N_mc, ...
            mean_err(c, 1), std_err(c, 1), ...
            mean_err(c, 2), std_err(c, 2), ...
            mean_err(c, 3), std_err(c, 3), ...
            total_viol_v(c), total_viol_y(c));
    end
    fprintf('========================================================================================\n');

    % ================== LaTeX 表格输出 ==================
    fprintf('\n=== LaTeX Table (paste into manuscript) ===\n\n');
    fprintf('\\begin{table}[htbp]\n');
    fprintf('    \\caption{Monte Carlo simulation results (N_{mc}=%d) for the final batch under different initial-state perturbations.}\n', N_mc);
    fprintf('    \\label{tab:mc_init_stats}\n');
    fprintf('    \\centering\n');
    fprintf('    \\begin{tabular}{c|ccc|c}\n');
    fprintf('        \\toprule\n');
    fprintf('        Initial State & $e_{20}$ (Mean $\\pm$ Std) & $e_{40}$ (Mean $\\pm$ Std) & $e_{60}$ (Mean $\\pm$ Std) & Constraint Violations \\\\\n');
    fprintf('        \\midrule\n');
    for c = 1:3
        fprintf('        %s ($s=%.2f$) & $%7.4f \\pm %.4f$ & $%7.4f \\pm %.4f$ & $%7.4f \\pm %.4f$ & %d (v), %d (y) \\\\\n', ...
            names{c}, s_vals(c), ...
            mean_err(c, 1), std_err(c, 1), ...
            mean_err(c, 2), std_err(c, 2), ...
            mean_err(c, 3), std_err(c, 3), ...
            total_viol_v(c), total_viol_y(c));
    end
    fprintf('        \\bottomrule\n');
    fprintf('    \\end{tabular}\n');
    fprintf('    \\par\\smallskip\n');
    fprintf('    \\footnotesize Note: Errors are calculated at the final batch (%d). Mean and standard deviation are computed over $N_{mc}=%d$ independent runs with different noise realizations and random perturbation directions. Constraint violations are summed over all $N_{mc}$ runs.\n', ...
        Params.N_batches, N_mc);
    fprintf('\\end{table}\n\n');

    % 保存结果到工作区，方便后续画图或检查
    assignin('base', 'mc_mean_err', mean_err);
    assignin('base', 'mc_std_err', std_err);
    assignin('base', 'mc_total_viol_v', total_viol_v);
    assignin('base', 'mc_total_viol_y', total_viol_y);
end
