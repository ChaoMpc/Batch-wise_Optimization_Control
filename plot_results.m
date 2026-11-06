function Figures = plot_results(Results)
% =========================================================================
% plot_results -- 7 Comparison Figures (Adapted Style & Units)
% =========================================================================
    Params      = Results.Params;
    Model       = Results.Model;
    Tube        = Results.Tube;
    Data_prop   = Results.Data_prop;
    Data_xue    = Results.Data_xue;

    Data_oh = [];
    if isfield(Results, 'Data_oh')
        Data_oh = Results.Data_oh;
    end

    Val_data = []; Validation = [];
    if isfield(Results, 'Validation') && isstruct(Results.Validation)
        Validation = Results.Validation;
    end
    if isfield(Results, 'Val_data'), Val_data = Results.Val_data; end

    y_opt_track = [];
    if isfield(Results, 'Optimal') && isfield(Results.Optimal, 'selected_outputs')
        y_opt_track = Results.Optimal.selected_outputs';
    end

    N         = Params.Ble;
    N_batches = Params.N_batches;
    n_in_phys = Params.n_in;
    n_in_opt  = Params.n_in_opt;
    n_out     = Params.n_out;
    trk       = Params.track_idx;
    if isempty(trk) || trk < 1 || trk > n_out, trk = 0; end

    y_base = Params.y_base;
    t_y = (0:N)'   * Params.Ts;
    t_v = (1:N)'   * Params.Ts;

    all_batches = 1:N_batches;

    % 严格补全所有输入输出单位
    in_names  = {'T_{j,sp} (K)','F_c (mL/s)','F_{anti} (mL/s)', ...
                 'F_{sol} (mL/s)','N_{stir} (rpm)','F_{seed} (g/min)'};
    out_names = {'T_1 (K)','T_2 (K)','T_3 (K)','C_{ATR} (kg/kg)','C_{NIR} (kg/kg)','C_{\rho} (kg/kg)'};

    if trk > 0
        if isfield(Params, 'y_ref_phys')
            y_ref_phys = Params.y_ref_phys(:)';
        else
            y_ref_phys = Params.y_base(trk, Params.Kp + 1) + Params.Ref_points(:)';
        end
        t_ref = Params.Kp * Params.Ts;
    else
        y_ref_phys = []; t_ref = [];
    end

    has_tube  = ~isempty(Tube);
    has_model = ~isempty(Model);

    Figures = struct('handle',{},'name',{});
    fcount = 0;

%% =======================================================================
% FIG 1 : 比较与追踪 (同时绘制三种方法的 y1)
%% =======================================================================
    if trk > 0
        fh = figure('Name','Fig 1: Proposed vs Xue vs Oh tracking', 'Color','w', 'Position',[80, 80, 1000, 600]);
        hold on; grid on;

        % 1. 物理约束
        plot(t_y, Params.y_phys_max(trk)*ones(N+1,1), 'k--', 'LineWidth', 1.0, 'DisplayName', 'Constraints');
        plot(t_y, Params.y_phys_min(trk)*ones(N+1,1), 'k--', 'LineWidth', 1.0, 'HandleVisibility','off');

        % 2. 紧化约束
        if has_tube
            y_hi = Params.y_phys_max(trk) - Tube.delta_y(trk);
            y_lo = Params.y_phys_min(trk) + Tube.delta_y(n_out + trk);
            plot(t_y, y_hi*ones(N+1,1), 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Tightened Constraints');
            plot(t_y, y_lo*ones(N+1,1), 'r-.', 'LineWidth', 1.5, 'HandleVisibility','off');
        end

        batches_compare = [15, 16, 17, 20];
        batches_compare = batches_compare(batches_compare <= N_batches);
        styles = {'-', '--', '-.', ':'};

        % 3. 三种方法曲线
        for jb = 1:numel(batches_compare)
            b = batches_compare(jb);
            st = styles{mod(jb-1, numel(styles)) + 1};
            Y_act_prop = squeeze(Data_prop.Y_hist(trk,:,b));
            Y_act_xue  = squeeze(Data_xue.Y_hist(trk,:,b));

            % Proposed 使用蓝色 (b)
            plot(t_y, Y_act_prop, ['b', st], 'LineWidth', 1.5, 'DisplayName', sprintf('Prop B%d', b));
            % Xue 使用黑色 (k)
            plot(t_y, Y_act_xue,  ['k', st], 'LineWidth', 1.2, 'DisplayName', sprintf('Xue B%d', b));

            % Oh 使用洋红色 (m)
            if ~isempty(Data_oh)
                Y_act_oh = squeeze(Data_oh.Y_hist(trk,:,b));
                plot(t_y, Y_act_oh, ['m', st], 'LineWidth', 1.2, 'DisplayName', sprintf('Oh B%d', b));
            end
        end

        % 4. 参考点与最优可达点
        plot(t_ref, y_ref_phys, 'ro', 'MarkerSize', 8, 'MarkerFaceColor','r', 'DisplayName', 'True Ref');
        if ~isempty(y_opt_track)
            plot(t_ref, y_opt_track, 'gp', 'MarkerSize', 10, 'MarkerFaceColor','g', 'DisplayName', 'Opt Reachable');
            for jj = 1:length(t_ref)
                plot([t_ref(jj), t_ref(jj)], [y_ref_phys(jj), y_opt_track(jj)], 'k:', 'LineWidth', 1.2, 'HandleVisibility','off');
            end
        end

        xlabel('Time (min)', 'FontSize', 11);
        ylabel(out_names{trk}, 'FontSize', 11);
        title(sprintf('Comparison: Tracking %s', out_names{trk}), 'FontSize', 12);
        xlim([t_y(1), t_y(end)]);
        legend('Location', 'best', 'NumColumns', 2, 'FontSize', 9);

        fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig1_tracking_comparison';
    end

%% =======================================================================
% FIG 2 : PROPOSED ALL OUTPUTS
%% =======================================================================
    fh = figure('Name','Fig 2: Proposed all outputs', 'Color','w', 'Position',[100, 80, 1400, 800]);
    for j = 1:n_out
        subplot(2, 3, j); hold on; grid on;
        if j == trk
            plot(t_y, Params.y_phys_max(j)*ones(N+1,1), 'k--', 'LineWidth', 0.8);
            plot(t_y, Params.y_phys_min(j)*ones(N+1,1), 'k--', 'LineWidth', 0.8);
            if has_tube
                plot(t_y, (Params.y_phys_max(j) - Tube.delta_y(j))*ones(N+1,1), 'r-.', 'LineWidth', 1.0);
                plot(t_y, (Params.y_phys_min(j) + Tube.delta_y(n_out + j))*ones(N+1,1), 'r-.', 'LineWidth', 1.0);
            end
        end
        for b = all_batches
            Y_act = squeeze(Data_prop.Y_hist(j,:,b));
            if b == 15 || b == 16, plot(t_y, Y_act, '-', 'Color', [1, 0, 0], 'LineWidth', 1.2);
            elseif b == 17, plot(t_y, Y_act, '-', 'Color', [0, 0, 1], 'LineWidth', 1.2);
            else, plot(t_y, Y_act, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); end
        end
        if j == trk && ~isempty(y_ref_phys)
            plot(t_ref, y_ref_phys, 'ro', 'MarkerSize', 6, 'MarkerFaceColor','r');
            if ~isempty(y_opt_track), plot(t_ref, y_opt_track, 'gp', 'MarkerSize', 8, 'MarkerFaceColor','g'); end
        end
        xlabel('Time (min)', 'FontSize', 11); ylabel(out_names{j}, 'FontSize', 11);
        title(sprintf('%s', out_names{j}), 'FontSize', 11);
        xlim([t_y(1), t_y(end)]);
    end
    try, sgtitle('Proposed method - all outputs (Batches 15,16 Red, 17 Blue)', 'FontWeight', 'bold', 'FontSize', 14); catch, end
    fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig2_proposed_outputs';

%% =======================================================================
% FIG 3 : PROPOSED ALL INPUTS
%% =======================================================================
    fh = figure('Name','Fig 3: Proposed all inputs', 'Color','w', 'Position',[120, 60, 1400, 800]);
    for i = 1:n_in_phys
        subplot(2, 3, i); hold on; grid on;
        plot(t_v, Params.u_phys_max(i)*ones(N,1), 'k--', 'LineWidth', 0.8);
        plot(t_v, Params.u_phys_min(i)*ones(N,1), 'k--', 'LineWidth', 0.8);
        if i <= n_in_opt && has_tube
            plot(t_v, (Params.u_phys_max(i) - Tube.delta_v(i))*ones(N,1), 'r-.', 'LineWidth', 1.0);
            plot(t_v, (Params.u_phys_min(i) + Tube.delta_v(n_in_opt + i))*ones(N,1), 'r-.', 'LineWidth', 1.0);
        end
        for b = all_batches
            if i <= n_in_opt
                V_act = squeeze(Data_prop.V_hist(i,:,b));
                if b == 15 || b == 16, plot(t_v, V_act, '-', 'Color', [1, 0, 0], 'LineWidth', 1.2);
                elseif b == 17, plot(t_v, V_act, '-', 'Color', [0, 0, 1], 'LineWidth', 1.2);
                else, plot(t_v, V_act, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); end
            else
                plot(t_v, Params.u_base(i, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
            end
        end
        xlabel('Time (min)', 'FontSize', 11); ylabel(in_names{i}, 'FontSize', 11);
        title(sprintf('%s', in_names{i}), 'FontSize', 11);
        xlim([t_v(1), t_v(end)]);
    end
    try, sgtitle('Proposed method - all inputs (Batches 15,16 Red, 17 Blue)', 'FontWeight', 'bold', 'FontSize', 14); catch, end
    fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig3_proposed_inputs';

%% =======================================================================
% FIG 4 : XUE ALL OUTPUTS
%% =======================================================================
    fh = figure('Name','Fig 4: Xue all outputs', 'Color','w', 'Position',[160, 80, 1400, 800]);
    for j = 1:n_out
        subplot(2, 3, j); hold on; grid on;
        if j == trk
            plot(t_y, Params.y_phys_max(j)*ones(N+1,1), 'k--', 'LineWidth', 0.8);
            plot(t_y, Params.y_phys_min(j)*ones(N+1,1), 'k--', 'LineWidth', 0.8);
            if has_tube
                plot(t_y, (Params.y_phys_max(j) - Tube.delta_y(j))*ones(N+1,1), 'r-.', 'LineWidth', 1.0);
                plot(t_y, (Params.y_phys_min(j) + Tube.delta_y(n_out + j))*ones(N+1,1), 'r-.', 'LineWidth', 1.0);
            end
        end
        for b = all_batches
            Y_xue = squeeze(Data_xue.Y_hist(j,:,b));
            if b == 15 || b == 16, plot(t_y, Y_xue, '-', 'Color', [1, 0, 0], 'LineWidth', 1.2);
            elseif b == 17, plot(t_y, Y_xue, '-', 'Color', [0, 0, 1], 'LineWidth', 1.2);
            else, plot(t_y, Y_xue, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); end
        end
        if j == trk && ~isempty(y_ref_phys)
            plot(t_ref, y_ref_phys, 'ro', 'MarkerSize', 6, 'MarkerFaceColor','r');
            if ~isempty(y_opt_track), plot(t_ref, y_opt_track, 'gp', 'MarkerSize', 8, 'MarkerFaceColor','g'); end
        end
        xlabel('Time (min)', 'FontSize', 11); ylabel(out_names{j}, 'FontSize', 11);
        title(sprintf('%s', out_names{j}), 'FontSize', 11);
        xlim([t_y(1), t_y(end)]);
    end
    try, sgtitle('Xue method - all outputs (Batches 15,16 Red, 17 Blue)', 'FontWeight', 'bold', 'FontSize', 14); catch, end
    fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig4_xue_outputs';

%% =======================================================================
% FIG 5 : XUE ALL INPUTS
%% =======================================================================
    fh = figure('Name','Fig 5: Xue all inputs', 'Color','w', 'Position',[180, 60, 1400, 800]);
    for i = 1:n_in_phys
        subplot(2, 3, i); hold on; grid on;
        plot(t_v, Params.u_phys_max(i)*ones(N,1), 'k--', 'LineWidth', 0.8);
        plot(t_v, Params.u_phys_min(i)*ones(N,1), 'k--', 'LineWidth', 0.8);
        for b = all_batches
            if i <= n_in_opt
                V_xue = squeeze(Data_xue.V_hist(i,:,b));
                if b == 15 || b == 16, plot(t_v, V_xue, '-', 'Color', [1, 0, 0], 'LineWidth', 1.2);
                elseif b == 17, plot(t_v, V_xue, '-', 'Color', [0, 0, 1], 'LineWidth', 1.2);
                else, plot(t_v, V_xue, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); end
            else
                plot(t_v, Params.u_base(i, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
            end
        end
        xlabel('Time (min)', 'FontSize', 11); ylabel(in_names{i}, 'FontSize', 11);
        title(sprintf('%s (Xue)', in_names{i}), 'FontSize', 11);
        xlim([t_v(1), t_v(end)]);
    end
    try, sgtitle('Xue method - all inputs (Batches 15,16 Red, 17 Blue)', 'FontWeight', 'bold', 'FontSize', 14); catch, end
    fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig5_xue_inputs';

%% =======================================================================
% FIG 6 : OH ALL OUTPUTS (新增：完全仿照 Fig 2)
%% =======================================================================
    if ~isempty(Data_oh)
        fh = figure('Name','Fig 6: Oh all outputs', 'Color','w', 'Position',[200, 100, 1400, 800]);
        for j = 1:n_out
            subplot(2, 3, j); hold on; grid on;
            if j == trk
                plot(t_y, Params.y_phys_max(j)*ones(N+1,1), 'k--', 'LineWidth', 0.8);
                plot(t_y, Params.y_phys_min(j)*ones(N+1,1), 'k--', 'LineWidth', 0.8);
                if has_tube
                    plot(t_y, (Params.y_phys_max(j) - Tube.delta_y(j))*ones(N+1,1), 'r-.', 'LineWidth', 1.0);
                    plot(t_y, (Params.y_phys_min(j) + Tube.delta_y(n_out + j))*ones(N+1,1), 'r-.', 'LineWidth', 1.0);
                end
            end
            for b = all_batches
                Y_oh = squeeze(Data_oh.Y_hist(j,:,b));
                if b == 15 || b == 16, plot(t_y, Y_oh, '-', 'Color', [1, 0, 0], 'LineWidth', 1.2);
                elseif b == 17, plot(t_y, Y_oh, '-', 'Color', [0, 0, 1], 'LineWidth', 1.2);
                else, plot(t_y, Y_oh, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); end
            end
            if j == trk && ~isempty(y_ref_phys)
                plot(t_ref, y_ref_phys, 'ro', 'MarkerSize', 6, 'MarkerFaceColor','r');
                if ~isempty(y_opt_track), plot(t_ref, y_opt_track, 'gp', 'MarkerSize', 8, 'MarkerFaceColor','g'); end
            end
            xlabel('Time (min)', 'FontSize', 11); ylabel(out_names{j}, 'FontSize', 11);
            title(sprintf('%s', out_names{j}), 'FontSize', 11);
            xlim([t_y(1), t_y(end)]);
        end
        try, sgtitle('Oh method (Ref [12]) - all outputs (Batches 15,16 Red, 17 Blue)', 'FontWeight', 'bold', 'FontSize', 14); catch, end
        fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig6_oh_outputs';
    end

%% =======================================================================
% FIG 7 : OH ALL INPUTS (新增：完全仿照 Fig 3)
%% =======================================================================
    if ~isempty(Data_oh)
        fh = figure('Name','Fig 7: Oh all inputs', 'Color','w', 'Position',[220, 80, 1400, 800]);
        for i = 1:n_in_phys
            subplot(2, 3, i); hold on; grid on;
            plot(t_v, Params.u_phys_max(i)*ones(N,1), 'k--', 'LineWidth', 0.8);
            plot(t_v, Params.u_phys_min(i)*ones(N,1), 'k--', 'LineWidth', 0.8);
            if i <= n_in_opt && has_tube
                plot(t_v, (Params.u_phys_max(i) - Tube.delta_v(i))*ones(N,1), 'r-.', 'LineWidth', 1.0);
                plot(t_v, (Params.u_phys_min(i) + Tube.delta_v(n_in_opt + i))*ones(N,1), 'r-.', 'LineWidth', 1.0);
            end
            for b = all_batches
                if i <= n_in_opt
                    V_oh = squeeze(Data_oh.V_hist(i,:,b));
                    if b == 15 || b == 16, plot(t_v, V_oh, '-', 'Color', [1, 0, 0], 'LineWidth', 1.2);
                    elseif b == 17, plot(t_v, V_oh, '-', 'Color', [0, 0, 1], 'LineWidth', 1.2);
                    else, plot(t_v, V_oh, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8); end
                else
                    plot(t_v, Params.u_base(i, :), '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8);
                end
            end
            xlabel('Time (min)', 'FontSize', 11); ylabel(in_names{i}, 'FontSize', 11);
            title(sprintf('%s', in_names{i}), 'FontSize', 11);
            xlim([t_v(1), t_v(end)]);
        end
        try, sgtitle('Oh method (Ref [12]) - all inputs (Batches 15,16 Red, 17 Blue)', 'FontWeight', 'bold', 'FontSize', 14); catch, end
        fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig7_oh_inputs';
    end

%% =======================================================================
% FIG 8 : FREE-RUN VALIDATION (原图6)
%% =======================================================================
    if ~isempty(Validation) && isfield(Validation, 'Y_free') && has_model && ~isempty(Val_data)
        fh = figure('Name','Fig 8: Free-run validation','Color','w', 'Position',[100, 80, 1400, 800]);
        t_full = (0:N)' * Params.Ts;
        cols = Validation.output_columns; % 63:121
        for j = 1:n_out
            subplot(2, 3, j); hold on; grid on;
            Y_true_full = Val_data.Y(j, 1:N+1, 1);
            Y_pred_full = Y_true_full;
            Y_pred_full(cols) = Validation.Y_free(j, :, 1);
            plot(t_full, Y_true_full, 'k-',  'LineWidth', 1.5);
            plot(t_full, Y_pred_full, 'r--', 'LineWidth', 1.8);
            yl = ylim();
            plot([Params.ctrl_start * Params.Ts, Params.ctrl_start * Params.Ts], yl, 'b-.', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            ylim(yl);
            if isfield(Validation, 'NRMSE') && isfield(Validation, 'VAF')
                title(sprintf('%s  | NRMSE=%.3f  VAF=%.1f%%', out_names{j}, Validation.NRMSE(j), Validation.VAF(j)), 'FontSize', 11);
            else
                title(out_names{j}, 'FontSize', 11);
            end
            xlabel('Time (min)', 'FontSize', 11); ylabel(out_names{j}, 'FontSize', 11);
            xlim([0, Params.Ble * Params.Ts]);
            if j == 1, legend('Measured','Free-run (60-120)','Location','northeast','FontSize', 9); end
        end
        try, sgtitle('Free-run validation (Evaluation on 60-120 min)', 'FontWeight','bold', 'FontSize', 14); catch, end
        fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig8_free_run_validation';
    end

%% =======================================================================
% FIG 9 : ERROR CONVERGENCE (包含 Proposed, Xue, Oh 三种方法)
%% =======================================================================
    fh = figure('Name','Fig 9: Error Convergence','Color','w','Position',[240 200 900 550]);
    hold on; grid on;

    % 防御性检查，并将数据强制转为列向量
    e1_prop = Data_prop.Err_true_act(:);
    if isfield(Data_prop, 'Err_opt_nom')
        e3_prop = Data_prop.Err_opt_nom(:);
    else
        e3_prop = e1_prop; % 旧数据保底
    end
    b_prop = (1:length(e1_prop))';

    % ==========================================
    % 我们的方法 (蓝色)
    % ==========================================
    h1 = plot(b_prop, e1_prop, 'b-d', ...
        'LineWidth', 1.8, 'MarkerSize', 7);
    h2 = plot(b_prop, e3_prop, 'b:+', ...
        'LineWidth', 1.5, 'MarkerSize', 8);

    % ==========================================
    % Xue 方法 (黑色)
    % ==========================================
    e1_xue = [];
    if ~isempty(Data_xue)
        e1_xue = Data_xue.Err_true_act(:);
        b_xue  = (1:length(e1_xue))';
        h3 = plot(b_xue, e1_xue, 'k--s', ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    else
        h3 = plot(nan, nan, 'k--s', 'LineWidth', 1.5, 'MarkerSize', 6);
        b_xue = [];
    end

    % ==========================================
    % Oh 方法 (洋红色, 文献[12])
    % ==========================================
    e1_oh = [];
    h4 = plot(nan, nan, 'm-o', 'LineWidth', 1.5, 'MarkerSize', 6); % 占位
    h5 = plot(nan, nan, 'm:x', 'LineWidth', 1.5, 'MarkerSize', 8); % 占位

    if ~isempty(Data_oh)
        e1_oh = Data_oh.Err_true_act(:);
        if isfield(Data_oh, 'Err_opt_nom')
            e3_oh = Data_oh.Err_opt_nom(:);
        else
            e3_oh = e1_oh;
        end
        b_oh = (1:length(e1_oh))';

        % 覆盖原本的占位线
        h4 = plot(b_oh, e1_oh, 'm-o', 'LineWidth', 1.5, 'MarkerSize', 6);
    end

    xlabel('Batch Number', 'FontSize', 11);
    ylabel('Error Norm', 'FontSize', 11);
    title('Tracking Error Convergence Comparison', 'FontSize', 12);

    % 动态生成图例
    leg_handles = [h1, h2, h3, h4];
    leg_labels  = {'Proposed: Act vs True Ref', ...
                   'Proposed: Nom vs Opt Reachable', ...
                   'Xue: Act vs True Ref', ...
                   'Oh: Act vs True Ref'};

    legend(leg_handles, leg_labels, 'Location', 'northeast', 'FontSize', 9);

    % 计算 X 轴的上限
    max_batches = length(e1_prop);
    if ~isempty(Data_xue), max_batches = max(max_batches, length(Data_xue.Err_true_act)); end
    if ~isempty(Data_oh),  max_batches = max(max_batches, length(Data_oh.Err_true_act));  end
    xlim([1, max_batches]);

    % === 动态计算 Y 轴上限，包含所有曲线 ===
    max_err = max([e1_prop; e3_prop]);
    if ~isempty(e1_xue), max_err = max(max_err, max(e1_xue)); end
    if ~isempty(e1_oh),  max_err = max(max_err, max(e1_oh)); end

    if isempty(max_err) || max_err <= 0 || isnan(max_err)
        max_err = 0.01; % 保底防错
    end

    % 留 10% 的顶部余量
    ylim([0, max_err * 1.1]);

    fcount = fcount + 1; Figures(fcount).handle = fh; Figures(fcount).name = 'fig9_error_convergence';
    fprintf('plot_results: generated %d figures.\n', numel(Figures));
end
