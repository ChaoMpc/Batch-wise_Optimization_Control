function plot_figure16_redraw(Results, out_dir)
% Plot Figure 16: Tracking error convergence comparison of the three methods,
% with the theoretical error bound \bar{\epsilon}_{k_i} superimposed.

    Data_prop = Results.Data_prop;
    Params    = Results.Params;
    Model     = Results.Model;
    Tube      = Results.Tube;
    UncSets   = Results.UncSets;

    Kp       = Params.Kp(:);
    trk      = Params.track_idx;
    y_ref    = Params.y_ref_phys(:);
    N_batch  = Params.N_batches;

    % Get the optimal reachable points
    y_opt = Results.Optimal.selected_outputs(:);

    % Compute the theoretical bound \bar{\epsilon}_{k_i}
    eps_bar = compute_tube_epsilon_bar(trk, Tube, Model, UncSets, Params);
    if numel(eps_bar) == 1
        eps_bar = repmat(eps_bar(1), size(Kp));
    end

    % Handle weights in scalar or matrix form
    if isscalar(Params.Q_weight)
        w_trk = Params.Q_weight;
    else
        w_trk = Params.Q_weight(trk, trk);
    end
    sqrt_w = sqrt(w_trk);

    %% =====================================================================
    % Figure 1: Comparison of three methods (2 subplots, 5 lines each)
    %% =====================================================================
    fig = figure('Name','Figure 16 Redrawn: Comparison','Color','w', ...
        'Position',[100 100 1000 450*numel(Kp)]);

    % Fix Octave cropping issue
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPosition', [0 0 10 4.5*numel(Kp)]);
    set(fig, 'PaperSize', [10 4.5*numel(Kp)]);

    for j = 1:numel(Kp)
        k_abs = Kp(j);
        col   = k_abs + 1;

        subplot(numel(Kp),1,j); hold on; grid on;

        % 1. Weighted error between nominal output and optimal reachable point (Proposed)
        if isfield(Data_prop,'Y_nominal_hist')
            e_nom = zeros(N_batch,1);
            for b = 1:N_batch
                y_n = Data_prop.Y_nominal_hist(trk, col, b);
                e_nom(b) = sqrt_w * abs(y_n - y_opt(j));
            end
            plot(1:N_batch, max(e_nom,1e-14), 'r--+', 'LineWidth', 1.5, 'MarkerSize', 8, ...
                'DisplayName', 'Proposed: Nom vs Opt');
        end

        % 2. Weighted error between actual output and actual reference point (Proposed)
        e_act_prop = zeros(N_batch,1);
        for b = 1:N_batch
            y_a = Data_prop.Y_hist(trk, col, b);
            e_act_prop(b) = sqrt_w * abs(y_a - y_ref(j));
        end
        plot(1:N_batch, max(e_act_prop,1e-14), 'r--o', 'LineWidth', 1.5, 'MarkerSize', 7, ...
            'MarkerFaceColor', 'none', 'DisplayName', 'Proposed: Act vs Ref');

        % 3. Xue method
        if isfield(Results, 'Data_xue') && ~isempty(Results.Data_xue)
            e_act_xue = zeros(N_batch,1);
            for b = 1:N_batch
                y_x = Results.Data_xue.Y_hist(trk, col, b);
                e_act_xue(b) = sqrt_w * abs(y_x - y_ref(j));
            end
            plot(1:N_batch, max(e_act_xue,1e-14), 'k--+', 'LineWidth', 1.2, 'MarkerSize', 7, ...
                'DisplayName', 'Xue: Act vs Ref');
        end

        % 4. Oh method
        if isfield(Results, 'Data_oh') && ~isempty(Results.Data_oh)
            e_act_oh = zeros(N_batch,1);
            for b = 1:N_batch
                y_o = Results.Data_oh.Y_hist(trk, col, b);
                e_act_oh(b) = sqrt_w * abs(y_o - y_ref(j));
            end
            plot(1:N_batch, max(e_act_oh,1e-14), 'b--+', 'LineWidth', 1.2, 'MarkerSize', 7, ...
                'DisplayName', 'Oh: Act vs Ref');
        end

        % 5. Error bound
        plot([1, N_batch], [eps_bar(j), eps_bar(j)], 'k-', 'LineWidth', 2, ...
            'DisplayName', 'Thm 3 Bound');

        set(gca,'YScale','log');
        xlabel('Batch number', 'FontSize', 11);
        ylabel('Weighted tracking error', 'FontSize', 11);
        title(sprintf('Reference point k_i = %d (t = %.0f min), track C_{ATR}', ...
            k_abs, k_abs*Params.Ts), 'FontSize', 12);
        legend('Location','east', 'FontSize', 9);
        xlim([1, N_batch]);
    end

    % Save comparison figure
    if nargin >= 2 && ~isempty(out_dir)
        try
            print(fig, fullfile(out_dir,'fig16_comparison.png'), '-dpng', '-r300');
            print(fig, fullfile(out_dir,'fig16_comparison.pdf'), '-dpdf');
            print(fig, fullfile(out_dir,'fig16_comparison.eps'), '-depsc');
            fprintf('Figure 16 (comparison) successfully saved to:\n  %s\n', out_dir);
        catch ME
            warning('Failed to save comparison figure: %s', ME.message);
        end
    end % Close handle to prevent conflicts

end
