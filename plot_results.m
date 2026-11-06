function plot_results(Data, Params, varargin)
% =========================================================================
% plot_results -- Visualize DyPLS-IL-MPC results.
%
% Updated for the N-step lifting:
%   * Y_hist / Y_nominal_hist : y(1..N)         -> time axis t_y = (1:N)'*Ts
%   * V_hist / V_nominal_hist : v(0..N-1)       -> time axis t_v = (1:N)'*Ts
%     (v(k-1) is drawn at x = k, aligned with y(k))
%   * Free-run validation reversal uses the TRAINING statistics
%     (Model.Y_std / Model.Y_mean_dev / Model.y_ss), not Val_data's own.
%
% Legend updates:
%   * Figure 1-4: add "Nominal (dashed)" placeholder and rename batch
%     labels to "Batch b (actual)"; enlarged legend font.
%   * Figure 5: 'best' is not implemented in Octave, use 'northeast'.
%   * Figure 7-8: add "Batch 1 (solid)" and "Batch N (dotted)" placeholders.
% =========================================================================

%% 0. Parse inputs
Model = []; Tube = []; VerifyResult = []; metrics = []; ValResults = []; summary_metrics = [];
if numel(varargin) >= 1 && ~isempty(varargin{1}), Model        = varargin{1}; end
if numel(varargin) >= 2 && ~isempty(varargin{2}), Tube         = varargin{2}; end
if numel(varargin) >= 3 && ~isempty(varargin{3}), VerifyResult = varargin{3}; end
if numel(varargin) >= 4 && ~isempty(varargin{4}), metrics      = varargin{4}; end
if numel(varargin) >= 5 && ~isempty(varargin{5}), ValResults   = varargin{5}; end
if numel(varargin) >= 6 && ~isempty(varargin{6}), summary_metrics = varargin{6}; end

%% 1. Basic data
N         = Params.Ble;
N_batches = Params.N_batches;
n_out     = Params.n_out;
n_in      = Params.n_in;
Ts        = Params.Ts;

if isfield(Data,'y_ss') && ~isempty(Data.y_ss), y_ss = Data.y_ss(:); else, y_ss = Params.y_ss(:); end
if isfield(Data,'u_ss') && ~isempty(Data.u_ss), u_ss = Data.u_ss(:); else, u_ss = Params.u_ss(:); end

Y_dev_hist         = Data.Y_hist;              % y(1..N)
V_dev_hist         = Data.V_hist;              % v(0..N-1)
Y_dev_nominal_hist = Data.Y_nominal_hist;      % y(1..N), 3rd dim = 1..N_batches+1
V_dev_nominal_hist = Data.V_nominal_hist;      % v(0..N-1)

Err_actual  = Data.Error_norm_actual(:);
Err_nominal = Data.Error_norm_nominal(:);

y_phys_max = Params.y_phys_max(:); y_phys_min = Params.y_phys_min(:);
u_phys_max = Params.u_phys_max(:); u_phys_min = Params.u_phys_min(:);

Kp             = Params.Kp(:);
y_ref_phys     = Params.y_ref_phys(:);

Y_hist_phys         = bsxfun(@plus, Y_dev_hist,         reshape(y_ss, [n_out,1,1]));
Y_nominal_hist_phys = bsxfun(@plus, Y_dev_nominal_hist, reshape(y_ss, [n_out,1,1]));
V_hist_phys         = bsxfun(@plus, V_dev_hist,         reshape(u_ss, [n_in,1,1]));
V_nominal_hist_phys = bsxfun(@plus, V_dev_nominal_hist, reshape(u_ss, [n_in,1,1]));

%% 2. Tightened bounds (deviation -> physical)
y_tight_max = []; y_tight_min = []; u_tight_max = []; u_tight_min = [];
if ~isempty(Tube) && isfield(Tube,'h_y_tight') && ~isempty(Tube.h_y_tight)
    y_tight_max = Tube.h_y_tight(1:n_out)          + y_ss;
    y_tight_min = -Tube.h_y_tight(n_out+1:end)     + y_ss;
end
if ~isempty(Tube) && isfield(Tube,'h_v_tight') && ~isempty(Tube.h_v_tight)
    u_tight_max = Tube.h_v_tight(1:n_in)           + u_ss;
    u_tight_min = -Tube.h_v_tight(n_in+1:end)      + u_ss;
end

%% 3. Time axes and style
% y(1..N): t_y = (1:N)'*Ts
% v(0..N-1): t_v = (1:N)'*Ts, so that v(k-1) is drawn at x = k
t_y = (1:N)' * Ts;
t_v = (1:N)' * Ts;

colors = lines(6);
lw_nom = 1.0; lw_act = 1.5; lw_bnd = 0.8; lw_tight = 1.0;
ms_ref = 6; ms_star = 8; fs_leg = 10;   % enlarged legend font

desired_batches = [1, 2, 3, 5];   % 收敛快时展示前几批
batch_list = desired_batches(desired_batches <= N_batches);
if isempty(batch_list), batch_list = 1; end
n_show = numel(batch_list);

%% Figure 1 -- y1
figure('Name','Output Tracking y1','Color','w'); hold on; grid on;
i = 1;
h_bnd_hi = plot(t_y, y_phys_max(i)*ones(N,1), 'k--', 'LineWidth', lw_bnd);
plot(t_y, y_phys_min(i)*ones(N,1), 'k--', 'LineWidth', lw_bnd, 'HandleVisibility','off');
h_tight_hi = [];
if ~isempty(y_tight_max)
    h_tight_hi = plot(t_y, y_tight_max(i)*ones(N,1), 'm:', 'LineWidth', lw_tight);
    plot(t_y, y_tight_min(i)*ones(N,1), 'm:', 'LineWidth', lw_tight, 'HandleVisibility','off');
end
h_batch = zeros(1, n_show);
for j = 1:n_show
    b   = batch_list(j);
    clr = colors(mod(j-1, size(colors,1))+1, :);
    plot(t_y, Y_nominal_hist_phys(i,:,b+1), '--', 'Color', clr, ...
         'LineWidth', lw_nom, 'HandleVisibility','off');
    h_batch(j) = plot(t_y, Y_hist_phys(i,:,b), '-', 'Color', clr, 'LineWidth', lw_act);
end

% gray dashed placeholder for "Nominal (dashed)"
h_nom_repr = plot(nan, nan, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', lw_nom);

h_ref = []; h_star = [];
if ~isempty(Kp)
    h_ref = plot(t_y(Kp), y_ref_phys, 'ro', 'MarkerSize', ms_ref, 'MarkerFaceColor','r');
    if ~isempty(VerifyResult) && isfield(VerifyResult,'y_s_at_Kp')
        y_s_phys = VerifyResult.y_s_at_Kp(:, 1) + y_ss(1);
        h_star = plot(t_y(Kp), y_s_phys, 'bo', 'MarkerSize', ms_star, ...
            'MarkerFaceColor','none', 'MarkerEdgeColor','b', 'LineWidth', 1.5);
    end
end
xlabel('Sample k'); ylabel('y_1 (wt%)');
leg_handles = h_batch;
leg_labels  = arrayfun(@(bb) sprintf('Batch %d (actual)', bb), batch_list, 'UniformOutput', false);
leg_handles(end+1) = h_nom_repr;
leg_labels{end+1}  = 'Nominal (dashed)';
if ~isempty(h_ref),  leg_handles(end+1) = h_ref;  leg_labels{end+1} = 'Reference'; end
if ~isempty(h_star), leg_handles(end+1) = h_star; leg_labels{end+1} = 'y_s^*'; end
leg_handles = [leg_handles, h_bnd_hi]; leg_labels = [leg_labels, 'Original Constraint'];
if ~isempty(h_tight_hi), leg_handles(end+1) = h_tight_hi; leg_labels{end+1} = 'Tightened Constraint'; end
legend(leg_handles, leg_labels, 'Location','southeast', 'FontSize',fs_leg);
xlim([1, N]); ylim([y_phys_min(i)-0.1, y_phys_max(i)+0.1]);

%% Figure 2 -- y2
figure('Name','Output Tracking y2','Color','w'); hold on; grid on;
i = 2;
h_bnd_hi = plot(t_y, y_phys_max(i)*ones(N,1), 'k--', 'LineWidth', lw_bnd);
plot(t_y, y_phys_min(i)*ones(N,1), 'k--', 'LineWidth', lw_bnd, 'HandleVisibility','off');
h_tight_hi = [];
if ~isempty(y_tight_max)
    h_tight_hi = plot(t_y, y_tight_max(i)*ones(N,1), 'm:', 'LineWidth', lw_tight);
    plot(t_y, y_tight_min(i)*ones(N,1), 'm:', 'LineWidth', lw_tight, 'HandleVisibility','off');
end
h_batch = zeros(1, n_show);
for j = 1:n_show
    b   = batch_list(j);
    clr = colors(mod(j-1, size(colors,1))+1, :);
    plot(t_y, Y_nominal_hist_phys(i,:,b+1), '--', 'Color', clr, ...
         'LineWidth', lw_nom, 'HandleVisibility','off');
    h_batch(j) = plot(t_y, Y_hist_phys(i,:,b), '-', 'Color', clr, 'LineWidth', lw_act);
end
h_nom_repr = plot(nan, nan, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', lw_nom);

xlabel('Sample k'); ylabel('y_2 (wt%)');
leg_handles = [h_batch, h_nom_repr, h_bnd_hi];
leg_labels  = [arrayfun(@(bb) sprintf('Batch %d (actual)', bb), batch_list, ...
                        'UniformOutput', false), {'Nominal (dashed)', 'Original Constraint'}];
if ~isempty(h_tight_hi), leg_handles(end+1) = h_tight_hi; leg_labels{end+1} = 'Tightened Constraint'; end
legend(leg_handles, leg_labels, 'Location','southeast', 'FontSize',fs_leg);
xlim([1, N]);

%% Figure 3 -- v1
figure('Name','Control Inputs v1','Color','w'); hold on; grid on;
i = 1;
h_bnd_hi = plot(t_v, u_phys_max(i)*ones(N,1), 'k--','LineWidth',lw_bnd);
plot(t_v, u_phys_min(i)*ones(N,1), 'k--','LineWidth',lw_bnd, 'HandleVisibility','off');
h_tight_hi = [];
if ~isempty(u_tight_max)
    h_tight_hi = plot(t_v, u_tight_max(i)*ones(N,1),'m:','LineWidth',lw_tight);
    plot(t_v, u_tight_min(i)*ones(N,1),'m:','LineWidth',lw_tight, 'HandleVisibility','off');
end
h_batch = zeros(1, n_show);
for j = 1:n_show
    b   = batch_list(j);
    clr = colors(mod(j-1, size(colors,1))+1, :);
    plot(t_v, V_nominal_hist_phys(i,:,b+1), '--', 'Color', clr, ...
         'LineWidth', lw_nom, 'HandleVisibility','off');
    h_batch(j) = plot(t_v, V_hist_phys(i,:,b), '-', 'Color', clr, 'LineWidth', lw_act);
end
h_nom_repr = plot(nan, nan, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', lw_nom);

xlabel('Sample k'); ylabel('v_1 = R (lb/min)');
leg_handles = [h_batch, h_nom_repr, h_bnd_hi];
leg_labels  = [arrayfun(@(bb) sprintf('Batch %d (actual)', bb), batch_list, ...
                        'UniformOutput', false), {'Nominal (dashed)', 'Original Constraint'}];
if ~isempty(h_tight_hi), leg_handles(end+1) = h_tight_hi; leg_labels{end+1} = 'Tightened Constraint'; end
legend(leg_handles, leg_labels, 'Location','southeast', 'FontSize',fs_leg);
xlim([1, N]);

%% Figure 4 -- v2
figure('Name','Control Inputs v2','Color','w'); hold on; grid on;
i = 2;
h_bnd_hi = plot(t_v, u_phys_max(i)*ones(N,1), 'k--','LineWidth',lw_bnd);
plot(t_v, u_phys_min(i)*ones(N,1), 'k--','LineWidth',lw_bnd, 'HandleVisibility','off');
h_tight_hi = [];
if ~isempty(u_tight_max)
    h_tight_hi = plot(t_v, u_tight_max(i)*ones(N,1),'m:','LineWidth',lw_tight);
    plot(t_v, u_tight_min(i)*ones(N,1),'m:','LineWidth',lw_tight, 'HandleVisibility','off');
end
h_batch = zeros(1, n_show);
for j = 1:n_show
    b   = batch_list(j);
    clr = colors(mod(j-1, size(colors,1))+1, :);
    plot(t_v, V_nominal_hist_phys(i,:,b+1), '--', 'Color', clr, ...
         'LineWidth', lw_nom, 'HandleVisibility','off');
    h_batch(j) = plot(t_v, V_hist_phys(i,:,b), '-', 'Color', clr, 'LineWidth', lw_act);
end
h_nom_repr = plot(nan, nan, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', lw_nom);

xlabel('Sample k'); ylabel('v_2 = S (lb/min)');
leg_handles = [h_batch, h_nom_repr, h_bnd_hi];
leg_labels  = [arrayfun(@(bb) sprintf('Batch %d (actual)', bb), batch_list, ...
                        'UniformOutput', false), {'Nominal (dashed)', 'Original Constraint'}];
if ~isempty(h_tight_hi), leg_handles(end+1) = h_tight_hi; leg_labels{end+1} = 'Tightened Constraint'; end
legend(leg_handles, leg_labels, 'Location','southeast', 'FontSize',fs_leg);
xlim([1, N]);

%% Figure 5 -- error convergence (Per-reference-point weighted error, Log scale)
figure('Name','Error Convergence','Color','w'); hold on; grid on;

n_ref = numel(Params.Kp);
q_trk = Params.Q_weight(1, 1); % 提取 y1 的追踪权重 Q_trk(1,1)
sqrt_q = sqrt(q_trk);

% 预分配误差数组
Err_act_pt = zeros(N_batches, n_ref);
Err_nom_pt = zeros(N_batches, n_ref);

% 计算每个参考点、每个批次的加权误差
for b = 1:N_batches
    for i = 1:n_ref
        k_i = Params.Kp(i);
        % 实际误差 (偏差量 - 参考偏差量) * sqrt(Q_trk)
        Err_act_pt(b, i) = sqrt_q * (Data.Y_hist(1, k_i, b) - Params.Ref_points(i));
        % 名义误差 (偏差量 - 参考偏差量) * sqrt(Q_trk)
        Err_nom_pt(b, i) = sqrt_q * (Data.Y_nominal_hist(1, k_i, b+1) - Params.Ref_points(i));
    end
end

markers = {'o', '+', 's'};

% 绘制实际误差 (蓝色实线，空心标记)
h_act = zeros(1, n_ref);
for i = 1:n_ref
    h_act(i) = semilogy(1:N_batches, max(abs(Err_act_pt(:,i)), 1e-15), ['b:', markers{i}], ...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b');
end

% 绘制名义误差 (红色实线，空心标记)
h_nom = zeros(1, n_ref);
for i = 1:n_ref
    h_nom(i) = semilogy(1:N_batches, max(abs(Err_nom_pt(:,i)), 1e-15), ['r:', markers{i}], ...
        'LineWidth', 1.5, 'MarkerSize', 5, 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'r');
end

% 计算 Theorem 3 的误差上界 (黑色实线)
if ~isempty(Tube) && isfield(Tube, 'delta_y')
    eps_bound = sqrt_q * Tube.delta_y(1);
    h_bound = semilogy([1, N_batches], [eps_bound, eps_bound], 'k-', 'LineWidth', 2.0);
    fprintf('plot_results: Theorem 3 bound (per point) = %.4e\n', eps_bound);
else
    eps_bound = 1e-5; % 占位防错
    h_bound = semilogy(nan, nan, 'k-', 'LineWidth', 2.0);
end

% 设置图例 (使用 TeX 语法)
leg_handles = [h_act, h_nom, h_bound];
leg_labels  = cell(1, length(leg_handles));

for i = 1:n_ref
    leg_labels{i} = sprintf('Act weighted e_{%d}', Params.Kp(i));
end
for i = 1:n_ref
    leg_labels{n_ref + i} = sprintf('Nom weighted e_{%d}', Params.Kp(i));
end
leg_labels{end} = 'Thm 3 Bound';

% 调用 legend 并尝试设置解释器
lgd = legend(leg_handles, leg_labels, 'Location', 'east', 'FontSize', fs_leg);
try
    set(lgd, 'Interpreter', 'tex');
catch
    % 若当前环境不支持则使用默认
end

xlim([1, N_batches]); xticks(1:N_batches);
xlabel('Batch number'); ylabel('Weighted Error (Log Scale)');

% ================= 动态设置 Y 轴范围 =================
% 收集所有要绘制的有效数据
all_vals = [abs(Err_act_pt(:)); abs(Err_nom_pt(:)); eps_bound];
all_vals = all_vals(all_vals > 0);

if ~isempty(all_vals)
    min_y = min(all_vals);
    max_y = max(all_vals);
    % 上下各留出一个数量级 (10倍) 的余量
    ylim([min_y / 10, max_y * 10]);
else
    ylim([1e-6, 1e2]);
end
% ======================================================

fprintf('plot_results: error norm summary (Weighted by Q_trk)\n');
for i = 1:n_ref
    fprintf('  Kp%d (Batch 1) : nominal = %.4e, actual = %.4e\n', i, Err_nom_pt(1,i), Err_act_pt(1,i));
    fprintf('  Kp%d (Batch %d) : nominal = %.4e, actual = %.4e\n', i, N_batches, Err_nom_pt(end,i), Err_act_pt(end,i));
end

%% Figure 6 -- free-run validation
if ~isempty(ValResults) && isfield(ValResults, 'Y0_sim')
    Val_data  = ValResults.Val_data;
    Y0_sim    = ValResults.Y0_sim;                 % windowed (start_idx .. end_idx)
    nrmses    = ValResults.nrmses;
    vafs      = ValResults.vafs;
    start_idx = ValResults.start_idx;
    window_len = ValResults.window_len;
    end_idx   = min(start_idx + window_len - 1, size(Val_data.Y0,1));
    t_window  = start_idx : end_idx;

    % Use TRAINING statistics for inverse standardization
    if isfield(Model, 'Y_std') && isfield(Model, 'Y_mean_dev') && isfield(Model, 'y_ss')
        Y_std_use      = Model.Y_std(:);
        Y_mean_dev_use = Model.Y_mean_dev(:);
        y_ss_use       = Model.y_ss(:);
    else
        Y_std_use      = Val_data.Y_std(:);
        Y_mean_dev_use = Val_data.Y_mean_dev(:);
        y_ss_use       = Val_data.y_ss(:);
    end

    Y_true_phys = bsxfun(@plus, ...
        bsxfun(@times, Val_data.Y0(t_window,:), Y_std_use'), ...
        (Y_mean_dev_use + y_ss_use)');
    Y_sim_phys  = bsxfun(@plus, ...
        bsxfun(@times, Y0_sim, Y_std_use'), ...
        (Y_mean_dev_use + y_ss_use)');

    t_plot = t_window - 1;   % sample indices 0..window_len-1

    figure('Name', 'Free-run Validation', 'Color', 'w');
    for i = 1:n_out
        subplot(n_out, 1, i); hold on; grid on;
        plot(t_plot, Y_true_phys(:,i), 'b-', 'LineWidth', 1.5);
        plot(t_plot, Y_sim_phys(:,i),  'r--', 'LineWidth', 1.5);
        xlabel('Sample k'); ylabel(sprintf('y_%d (wt%%)', i));
        title(sprintf('y_%d | NRMSE = %.4f, VAF = %.2f%%', i, nrmses(i), vafs(i)), ...
              'FontSize', fs_leg);
        legend('True (validation)', 'Free-run Simulation', 'Location', 'best', 'FontSize', fs_leg);
        xlim([min(t_plot), max(t_plot)]);
    end
end

%% Figure 7 -- initial state y1 (Multi-perturbation test)
if isfield(Data, 'init_test')
    figure('Name', 'Initial State Test y1', 'Color', 'w');
    ax_main = axes(); hold(ax_main, 'on'); grid(ax_main, 'on');

    i = 1;
    % 约束线（不放入 legend，保持画面整洁）
    plot(ax_main, t_y, y_phys_max(i)*ones(N,1), 'k--', 'LineWidth', lw_bnd);
    plot(ax_main, t_y, y_phys_min(i)*ones(N,1), 'k--', 'LineWidth', lw_bnd, 'HandleVisibility','off');
    if ~isempty(y_tight_max)
        plot(ax_main, t_y, y_tight_max(i)*ones(N,1), 'm:', 'LineWidth', lw_tight);
        plot(ax_main, t_y, y_tight_min(i)*ones(N,1), 'm:', 'LineWidth', lw_tight, 'HandleVisibility','off');
    end

    D_nom   = Data.init_test.Nominal;
    D_small = Data.init_test.Small;
    D_mid   = Data.init_test.Mid;
    D_edge  = Data.init_test.Edge;

    Y_nom_phys   = bsxfun(@plus, D_nom.Y_hist,   reshape(y_ss, [n_out,1,1]));
    Y_small_phys = bsxfun(@plus, D_small.Y_hist, reshape(y_ss, [n_out,1,1]));
    Y_mid_phys   = bsxfun(@plus, D_mid.Y_hist,   reshape(y_ss, [n_out,1,1]));
    Y_edge_phys  = bsxfun(@plus, D_edge.Y_hist,  reshape(y_ss, [n_out,1,1]));

    % 绘制主图曲线
    h_nom   = plot(ax_main, t_y, Y_nom_phys(i,:,1),   'k-',  'LineWidth', 1.5);
    h_small1 = plot(ax_main, t_y, Y_small_phys(i,:,1), 'b-',  'LineWidth', 1.5);
    h_mid1   = plot(ax_main, t_y, Y_mid_phys(i,:,1),   'g-',  'LineWidth', 1.5);
    h_edge1  = plot(ax_main, t_y, Y_edge_phys(i,:,1),  'r-',  'LineWidth', 1.5);
    h_smallN = plot(ax_main, t_y, Y_small_phys(i,:,end), 'b:', 'LineWidth', 2.0);
    h_midN   = plot(ax_main, t_y, Y_mid_phys(i,:,end),   'g:', 'LineWidth', 2.0);
    h_edgeN  = plot(ax_main, t_y, Y_edge_phys(i,:,end),  'r:', 'LineWidth', 2.0);

    % 绘制实际参考点（红色空心圆）
    h_ref = plot(ax_main, t_y(Params.Kp), Params.y_ref_phys, 'ro', ...
                 'MarkerSize', 8, 'MarkerFaceColor', 'none', 'LineWidth', 1.5);

    xlabel(ax_main, 'Sample k'); ylabel(ax_main, 'y_1 (wt%)'); xlim(ax_main, [1, N]);

    % 两列 legend
    handles = [h_nom, h_small1, h_mid1, h_edge1, h_smallN, h_midN, h_edgeN, h_ref];
    labels  = {'Nominal Init', 'Small (Batch 1)', 'Mid (Batch 1)', 'Edge (Batch 1)', ...
               'Small (Batch 15)', 'Mid (Batch 15)', 'Edge (Batch 15)', 'Reference'};
    try
        legend(ax_main, handles, labels, 'Location', 'northeast', 'FontSize', fs_leg, 'NumColumns', 2);
    catch
        legend(ax_main, handles, labels, 'Location', 'northeast', 'FontSize', fs_leg);
    end

    % ===== 添加 inset 子图（放大初始区域） =====
    ax_inset = axes('Position', [0.2, 0.2, 0.28, 0.28]); % 左下角
    hold(ax_inset, 'on'); grid(ax_inset, 'on');
    % 绘制 batch 1 的曲线（solid）
    plot(ax_inset, t_y, Y_nom_phys(i,:,1),   'k-', 'LineWidth', 1.0);
    plot(ax_inset, t_y, Y_small_phys(i,:,1), 'b-', 'LineWidth', 1.0);
    plot(ax_inset, t_y, Y_mid_phys(i,:,1),   'g-', 'LineWidth', 1.0);
    plot(ax_inset, t_y, Y_edge_phys(i,:,1),  'r-', 'LineWidth', 1.0);
    xlim(ax_inset, [1, 2]); ylim(ax_inset, [96.2, 96.3]);
    xlabel(ax_inset, 'k'); ylabel(ax_inset, 'y_1');
    title(ax_inset, 'Zoom [1,2]x[96.2,9.3]', 'FontSize', 8);
    set(ax_inset, 'FontSize', 8);
end

%% Figure 8 -- initial state y2 (Multi-perturbation test)
if isfield(Data, 'init_test')
    figure('Name', 'Initial State Test y2', 'Color', 'w');
    ax_main = axes(); hold(ax_main, 'on'); grid(ax_main, 'on');

    i = 2;
    % 约束线
    plot(ax_main, t_y, y_phys_max(i)*ones(N,1), 'k--', 'LineWidth', lw_bnd);
    plot(ax_main, t_y, y_phys_min(i)*ones(N,1), 'k--', 'LineWidth', lw_bnd, 'HandleVisibility','off');
    if ~isempty(y_tight_max)
        plot(ax_main, t_y, y_tight_max(i)*ones(N,1), 'm:', 'LineWidth', lw_tight);
        plot(ax_main, t_y, y_tight_min(i)*ones(N,1), 'm:', 'LineWidth', lw_tight, 'HandleVisibility','off');
    end

    % 绘制主图曲线
    h_nom   = plot(ax_main, t_y, Y_nom_phys(i,:,1),   'k-',  'LineWidth', 1.5);
    h_small1 = plot(ax_main, t_y, Y_small_phys(i,:,1), 'b-',  'LineWidth', 1.5);
    h_mid1   = plot(ax_main, t_y, Y_mid_phys(i,:,1),   'g-',  'LineWidth', 1.5);
    h_edge1  = plot(ax_main, t_y, Y_edge_phys(i,:,1),  'r-',  'LineWidth', 1.5);
    h_smallN = plot(ax_main, t_y, Y_small_phys(i,:,end), 'b:', 'LineWidth', 2.0);
    h_midN   = plot(ax_main, t_y, Y_mid_phys(i,:,end),   'g:', 'LineWidth', 2.0);
    h_edgeN  = plot(ax_main, t_y, Y_edge_phys(i,:,end),  'r:', 'LineWidth', 2.0);

    xlabel(ax_main, 'Sample k'); ylabel(ax_main, 'y_2 (wt%)'); xlim(ax_main, [1, N]);

    % 单列 legend（不包含参考点）
    handles = [h_nom, h_small1, h_mid1, h_edge1, h_smallN, h_midN, h_edgeN];
    labels  = {'Nominal Init', 'Small (Batch 1)', 'Mid (Batch 1)', 'Edge (Batch 1)', ...
               'Small (Batch 15)', 'Mid (Batch 15)', 'Edge (Batch 15)'};
    legend(ax_main, handles, labels, 'Location', 'southeast', 'FontSize', fs_leg);

    % ===== 添加 inset 子图（放大初始区域） =====
    ax_inset = axes('Position', [0.2, 0.2, 0.28, 0.28]); % 左下角
    hold(ax_inset, 'on'); grid(ax_inset, 'on');
    % 绘制 batch 1 的曲线（solid）
    plot(ax_inset, t_y, Y_nom_phys(i,:,1),   'k-', 'LineWidth', 1.0);
    plot(ax_inset, t_y, Y_small_phys(i,:,1), 'b-', 'LineWidth', 1.0);
    plot(ax_inset, t_y, Y_mid_phys(i,:,1),   'g-', 'LineWidth', 1.0);
    plot(ax_inset, t_y, Y_edge_phys(i,:,1),  'r-', 'LineWidth', 1.0);
    xlim(ax_inset, [1, 2]); ylim(ax_inset, [0.45, 0.55]);
    xlabel(ax_inset, 'k'); ylabel(ax_inset, 'y_2');
    title(ax_inset, 'Zoom [1,2]x[0.45,0.55]', 'FontSize', 8);
    set(ax_inset, 'FontSize', 8);
end

%% Validation table (coverage may be N/A)
if ~isempty(metrics)
    fprintf('\n\n========================================================================================\n');
    fprintf('  Table: Monte Carlo / Robustness Validation Summary\n');
    fprintf('========================================================================================\n');
    fprintf(' %-24s | %-8s | %-8s | %-10s | %-10s | %-10s\n', ...
            'Scenario', 'InViol', 'OutViol', 'Cov_Eta', 'Cov_R', 'QP Time');
    fprintf('----------------------------------------------------------------------------------------\n');

    cov_eta_str = ternary_str(~isnan(metrics.coverage_eta), ...
        sprintf('%.2f%%', 100*metrics.coverage_eta), 'N/A');
    cov_r_str   = ternary_str(~isnan(metrics.coverage_r), ...
        sprintf('%.2f%%', 100*metrics.coverage_r), 'N/A');
    fprintf(' %-24s | %-8d | %-8d | %-10s | %-10s | %-8.4fs\n', ...
        'Nominal Init (Main)', metrics.n_violations_v, metrics.n_violations_y, ...
        cov_eta_str, cov_r_str, mean(metrics.qp_times));

    if isfield(summary_metrics, 'pert_small')
        m_small = summary_metrics.pert_small;
        cov_eta_str = ternary_str(~isnan(m_small.coverage_eta), ...
            sprintf('%.2f%%', 100*m_small.coverage_eta), 'N/A');
        cov_r_str   = ternary_str(~isnan(m_small.coverage_r), ...
            sprintf('%.2f%%', 100*m_small.coverage_r), 'N/A');
        fprintf(' %-24s | %-8d | %-8d | %-10s | %-10s | %-8.4fs\n', ...
            'Small Perturbation', m_small.n_violations_v, m_small.n_violations_y, ...
            cov_eta_str, cov_r_str, mean(m_small.qp_times));
    end

    if isfield(summary_metrics, 'pert_mid')
        m_mid = summary_metrics.pert_mid;
        cov_eta_str = ternary_str(~isnan(m_mid.coverage_eta), ...
            sprintf('%.2f%%', 100*m_mid.coverage_eta), 'N/A');
        cov_r_str   = ternary_str(~isnan(m_mid.coverage_r), ...
            sprintf('%.2f%%', 100*m_mid.coverage_r), 'N/A');
        fprintf(' %-24s | %-8d | %-8d | %-10s | %-10s | %-8.4fs\n', ...
            'Mid Perturbation', m_mid.n_violations_v, m_mid.n_violations_y, ...
            cov_eta_str, cov_r_str, mean(m_mid.qp_times));
    end

    if isfield(summary_metrics, 'pert_edge')
        m_edge = summary_metrics.pert_edge;
        cov_eta_str = ternary_str(~isnan(m_edge.coverage_eta), ...
            sprintf('%.2f%%', 100*m_edge.coverage_eta), 'N/A');
        cov_r_str   = ternary_str(~isnan(m_edge.coverage_r), ...
            sprintf('%.2f%%', 100*m_edge.coverage_r), 'N/A');
        fprintf(' %-24s | %-8d | %-8d | %-10s | %-10s | %-8.4fs\n', ...
            'Edge Perturbation', m_edge.n_violations_v, m_edge.n_violations_y, ...
            cov_eta_str, cov_r_str, mean(m_edge.qp_times));
    end

    fprintf('========================================================================================\n');
    fprintf('  * Note: InViol/OutViol denote the number of samples violating physical constraints.\n');
    fprintf('  * Cov_Eta/Cov_R denote the empirical coverage of the uncertainty sets.\n');
    fprintf('  * N/A means no samples were evaluated (e.g., r when A_lv == n_out).\n\n');

    if ~exist('figures','dir'), mkdir('figures'); end
    fid = fopen(fullfile('figures', 'validation_table.csv'), 'w');
    if fid > 0
        fprintf(fid, 'Scenario,InViol,OutViol,Cov_Eta(%%),Cov_R(%%),Mean_QP_Time(s)\n');
        fprintf(fid, 'Nominal Init,%d,%d,%s,%s,%.4f\n', ...
            metrics.n_violations_v, metrics.n_violations_y, ...
            cov_eta_str, cov_r_str, mean(metrics.qp_times));
        if isfield(summary_metrics, 'pert_small')
            fprintf(fid, 'Small Perturbation,%d,%d,%.2f,%.2f,%.4f\n', ...
                m_small.n_violations_v, m_small.n_violations_y, ...
                100*m_small.coverage_eta, 100*m_small.coverage_r, mean(m_small.qp_times));
        end
        if isfield(summary_metrics, 'pert_mid')
            fprintf(fid, 'Mid Perturbation,%d,%d,%.2f,%.2f,%.4f\n', ...
                m_mid.n_violations_v, m_mid.n_violations_y, ...
                100*m_mid.coverage_eta, 100*m_mid.coverage_r, mean(m_mid.qp_times));
        end
        if isfield(summary_metrics, 'pert_edge')
            fprintf(fid, 'Edge Perturbation,%d,%d,%.2f,%.2f,%.4f\n', ...
                m_edge.n_violations_v, m_edge.n_violations_y, ...
                100*m_edge.coverage_eta, 100*m_edge.coverage_r, mean(m_edge.qp_times));
        end
        fclose(fid);
        fprintf('  [Saved] Validation table exported to figures/validation_table.csv\n');
    end
end

%% Tube diagnostics
if ~isempty(Tube)
    fprintf('\nplot_results: Tube diagnostics\n');
    if isfield(Tube,'spectral_radius'), fprintf('  Closed-loop spectral radius : %.6f\n', Tube.spectral_radius); end
    if isfield(Tube,'delta_v'), fprintf('  Input  delta_v              : [%s]\n', num2str(Tube.delta_v(:)','%.3e ')); end
    if isfield(Tube,'delta_y'), fprintf('  Output delta_y              : [%s]\n', num2str(Tube.delta_y(:)','%.3e ')); end
end

%% J_batch summary
if isfield(Data, 'J_batch') && ~isempty(Data.J_batch)
    fprintf('\nplot_results: J_batch summary (Theorem 2 consistency)\n');
    fprintf('  J_batch(1)   = %.4e\n', Data.J_batch(1));
    fprintf('  J_batch(end) = %.4e\n', Data.J_batch(end));
    dJ = diff(Data.J_batch);
    if all(dJ <= 1e-8 * max(1, abs(Data.J_batch(1))))
        fprintf('  J_batch is non-increasing.\n');
    else
        fprintf('  WARNING: J_batch is not strictly non-increasing.\n');
    end
end

%% Tracking error table
n_ref = numel(Params.Kp);
E_act_pt = zeros(N_batches, n_ref);
E_nom_pt = zeros(N_batches, n_ref);
for b = 1:N_batches
    for j = 1:n_ref
        k_j = Params.Kp(j);
        E_act_pt(b,j) = Y_dev_hist(1, k_j, b)           - Params.Ref_points(j);
        E_nom_pt(b,j) = Y_dev_nominal_hist(1, k_j, b+1) - Params.Ref_points(j);
    end
end
E_act_total = sqrt(sum(E_act_pt.^2, 2));
E_nom_total = sqrt(sum(E_nom_pt.^2, 2));

fprintf('\nplot_results: tracking error table (wt%%)\n');
fprintf('%6s', 'batch');
for j = 1:n_ref, fprintf('%12s', sprintf('act Kp%d', j)); end
for j = 1:n_ref, fprintf('%12s', sprintf('nom Kp%d', j)); end
fprintf('%14s%14s\n', 'act total', 'nom total');

for b = 1:N_batches
    fprintf('%6d', b);
    for j = 1:n_ref, fprintf('%12.4f', E_act_pt(b,j)); end
    for j = 1:n_ref, fprintf('%12.4f', E_nom_pt(b,j)); end
    fprintf('%14.4f%14.4f\n', E_act_total(b), E_nom_total(b));
end

if ~exist('figures','dir'), mkdir('figures'); end
save('-mat', fullfile('figures','tracking_error_table.mat'), ...
    'E_act_pt', 'E_nom_pt', 'E_act_total', 'E_nom_total');
fprintf('plot_results: plotting completed.\n');

end

% ========================================================================
function s = ternary_str(cond, a, b)
    if cond, s = a; else, s = b; end
end
