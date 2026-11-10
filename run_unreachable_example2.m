% =========================================================================
% run_unreachable_example2.m
%
% Unreachable reference scenario for Example 2.
%
% This script directly loads the pre-computed artifacts (Stage 1-3) to save
% time. It overrides ONLY the 120-min reference point to 0.013, making it
% physically unreachable for the identified model. It then builds the new
% limiting nominal problem, runs the proposed method, and checks the bound
% \bar{\epsilon}_{k_i} (Theorem 3) against the actual tracking error.
% =========================================================================

clear; clc; close all;
if ~exist('OCTAVE_VERSION','builtin')
    error('This script requires GNU Octave.');
end
try, pkg load control; catch, end
addpath(pwd); rehash();

%% 1. Load pre-computed artifacts
fprintf('\n=== Loading Stage 1-3 artifacts ===\n');
res_dir = fullfile(pwd, 'results');

if exist(fullfile(res_dir, 'stage01_setup.mat'), 'file') ~= 2 || ...
   exist(fullfile(res_dir, 'stage02_identification.mat'), 'file') ~= 2 || ...
   exist(fullfile(res_dir, 'stage03_tube.mat'), 'file') ~= 2
    error('Missing Stage 1-3 files. Please run run_stage01, 02, 03 first.');
end

S1 = load(fullfile(res_dir, 'stage01_setup.mat'));
S2 = load(fullfile(res_dir, 'stage02_identification.mat'));
S3 = load(fullfile(res_dir, 'stage03_tube.mat'));

Params   = S1.Params;
Reference= S1.Reference;
Model    = S2.Model;
Tube     = S3.Tube;
UncSets  = S3.UncSets;

% Override the reference point (120 min)
Params.y_ref_phys = [0.022; 0.013];

[Reference, Params] = generate_reference(Params);

trk = Params.track_idx;   % C_ATR channel (4)
ni  = Params.n_in_opt;    % 5 optimized inputs
no  = Params.n_out;

C_ATR_lo = Params.y_phys_min(trk) + Tube.delta_y(no + trk);
C_ATR_hi = Params.y_phys_max(trk) - Tube.delta_y(trk);
fprintf('\nC_ATR tightened interval: [%.6f, %.6f]\n', C_ATR_lo, C_ATR_hi);
fprintf('Reference [90, 120]     : [%.6f, %.6f]\n', ...
        Params.y_ref_phys(1), Params.y_ref_phys(2));
if Params.y_ref_phys(2) >= C_ATR_lo
    error('Reference at 120 min is still reachable.');
end
fprintf('Confirmed: ref(120) = %.6f < tightened lower = %.6f  -> UNREACHABLE\n', ...
        Params.y_ref_phys(2), C_ATR_lo);

%% 2. Rebuild Problem and Optimal for the new reference
fprintf('\n=== Rebuilding nominal problem for unreachable target ===\n');
Problem = build_nominal_problem(Model, Tube, Params);

% Ensure feasibility seed exists (it should from Stage 3, but we refresh it)
Feasibility = check_tightened_feasibility(Model, Tube, Params);
Tube.nominal_feasible = Feasibility.feasible;
Tube.feasible_seed    = Feasibility.t_feasible;

Initial = compute_initial_feasible(Problem, Tube, Params);
Optimal = compute_optimal_reachable(Problem, Initial.t, Params);

y_star_k = Optimal.selected_outputs(:);       % (nKp x 1)
fprintf('\nOptimal reachable points (y*_s):\n');
for j = 1:numel(Params.Kp)
    fprintf('  k=%3d min: ref=%.6f, y*_s=%.6f, gap=%+.3e\n', ...
            Params.Kp(j), Params.y_ref_phys(j), y_star_k(j), ...
            Params.y_ref_phys(j) - y_star_k(j));
end

%% 3. Run proposed method
fprintf('\n=== Running proposed Tube method (unreachable scenario) ===\n');
Scenario  = generate_simulation_scenario(Params);
Data_prop = run_ilmpc(Problem, Initial, Optimal, Scenario, ...
                      Model, Tube, UncSets, Params);

%% 4. Compute bar-epsilon_{k_i} (Theorem 3)
Sy   = diag(Model.Y_std);
Qout = Model.Q;
Pe   = Tube.P;
Wr   = UncSets.W_r;
Cy   = Sy * Qout * Model.barC;         % Physical output error map

bar_epsilon = zeros(numel(Params.Kp), 1);

if isscalar(Params.Q_weight)
    q = Params.Q_weight;
else
    q = Params.Q_weight(trk, trk);
end
sqrt_q = sqrt(q);

for j = 1:numel(Params.Kp)
    row = Cy(trk, :);
    s1 = sqrt(max(0, row * Pe * row'));                 % e in E(Pe)
    s2 = sqrt(max(0, Sy(trk,trk)^2 * Wr(trk,trk)));     % r in E(Wr)
    bar_epsilon(j) = sqrt_q * (s1 + s2);
end

fprintf('\n======================================================\n');
fprintf('【Verification】bar-epsilon_{k_i} (Theorem 3)\n');
fprintf('  k_i = %d : %.12e\n', Params.Kp(1), bar_epsilon(1));
fprintf('  k_i = %d : %.12e\n', Params.Kp(2), bar_epsilon(2));
fprintf('  (Compare this with the black line in previous Figure 18)\n');
fprintf('======================================================\n');

out_dir = fullfile(pwd, 'unreachable_results');
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end

%% 5. Per-batch, per-k_i Q_trk-weighted errors
B   = Params.N_batches;
nKp = numel(Params.Kp);

Err_nom_ki = zeros(nKp, B);
Err_act_ki = zeros(nKp, B);

for j = 1:nKp
    y_star_j = y_star_k(j);
    col_phys = Params.Kp(j) + 1;    % 1-based column
    for b = 1:B
        y_nom_b = Data_prop.Y_nominal_hist(trk, col_phys, b);
        y_act_b = Data_prop.Y_hist(trk,         col_phys, b);
        Err_nom_ki(j, b) = sqrt_q * abs(y_nom_b - y_star_j);
        Err_act_ki(j, b) = sqrt_q * abs(y_act_b - y_star_j);
    end
end


%% ========================================================================
%% 6. PLOTTING SECTION - Figure 1: Figure 18-style (Vertical Subplots)
%% ========================================================================
fig1 = figure('Name','Unreachable reference scenario (Errors)','Color','w', ...
             'Position',[100, 100, 800, 800]);

% 统一尺寸设置，防止 EPS/PDF 导出变形
set(fig1, 'PaperUnits', 'inches');
set(fig1, 'PaperPosition', [0 0 8 8]);
set(fig1, 'PaperSize', [8 8]);

batch_idx = (1:B)';

for j = 1:nKp
    subplot(nKp, 1, j); hold on; grid on; box on;
    set(gca, 'LineWidth', 1.2);

    % 使用 1e-12 防止 log 坐标下出现 -inf
    % 【颜色修正】Nominal 用红色圆圈，Actual 用蓝色方块（与轨迹图保持一致）
    h_nom = semilogy(batch_idx, max(Err_nom_ki(j,:), 1e-12), ...
                     'r-o', 'LineWidth', 1.5, 'MarkerSize', 7);
    h_act = semilogy(batch_idx, max(Err_act_ki(j,:), 1e-12), ...
                     'b-s', 'LineWidth', 1.5, 'MarkerSize', 7);
    h_eps = semilogy(batch_idx, bar_epsilon(j) * ones(B,1), ...
                     'k-', 'LineWidth', 2);

    xlabel('Batch number', 'FontSize', 11);
    ylabel('Q_{trk}-weighted tracking error', 'FontSize', 11);
    title(sprintf('k_i = %d min (y*_s = %.6f)', ...
                  Params.Kp(j), y_star_k(j)), 'FontSize', 12);

    legend([h_nom, h_act, h_eps], ...
           {'Nominal error to y*_s', 'Actual error to y*_s', ...
            'Thm 3 Bound'}, ...
           'Location', 'southwest', 'FontSize', 11);

    xlim([0.5, B + 0.5]);

    finite_vals = [bar_epsilon(j), Err_nom_ki(j,:), Err_act_ki(j,:)];
    finite_vals = finite_vals(isfinite(finite_vals) & finite_vals > 0);
    if isempty(finite_vals)
        ymax = 1;
    else
        ymax = max(finite_vals) * 10;
    end
    ylim([1e-12, ymax]);
end

saveas(fig1, fullfile(out_dir, 'fig_unreachable_tracking.png'));
try, print(fig1, fullfile(out_dir, 'fig_unreachable_tracking.eps'), '-depsc2'); catch, end


%% ========================================================================
%% 7. PLOTTING SECTION - Figure 2: Figure 11-style (Trajectories + Constraints + Zoom-in)
%% ========================================================================
fig2 = figure('Name','Unreachable trajectory comparison', 'Color','w', ...
              'Position',[100, 150, 1100, 600]);
hold on; grid on; box on; % Added box on
set(gca, 'LineWidth', 1.2); % Thicker axes border

t_y = (0:Params.Ble) * Params.Ts;
batches_to_plot = [1, 3, 5, 20];
styles = {'-', '--', '-.', ':'};

% Extract bounds
y_phys_max = Params.y_phys_max(trk);
y_phys_min = Params.y_phys_min(trk);
y_tight_max = Params.y_phys_max(trk) - Tube.delta_y(trk);
y_tight_min = Params.y_phys_min(trk) + Tube.delta_y(no + trk);

% 准备图例句柄和字符串
h_leg = [];
leg_str = {};

% --- Main Plot ---
% 1. Plot Physical & Tightened Constraints (仅上界加图例，下界与之对称)
h_phy = plot([0, Params.Ble*Params.Ts], [y_phys_max, y_phys_max], 'k--', 'LineWidth', 1.5);
plot([0, Params.Ble*Params.Ts], [y_phys_min, y_phys_min], 'k--', 'LineWidth', 1.5);
h_tight = plot([0, Params.Ble*Params.Ts], [y_tight_max, y_tight_max], 'm-.', 'LineWidth', 1.5);
plot([0, Params.Ble*Params.Ts], [y_tight_min, y_tight_min], 'm-.', 'LineWidth', 1.5);

% 2. Plot Trajectories and record handles
for idx = 1:length(batches_to_plot)
    b = batches_to_plot(idx);
    st = styles{idx};

    h_nom_i = plot(t_y, squeeze(Data_prop.Y_nominal_hist(trk, :, b)), ['r', st], 'LineWidth', 1.5);
    h_act_i = plot(t_y, squeeze(Data_prop.Y_hist(trk, :, b)), ['b', st], 'LineWidth', 1.5);

    h_leg = [h_leg, h_nom_i, h_act_i];
    leg_str{end+1} = sprintf('Nominal B%d', b);
    leg_str{end+1} = sprintf('Actual B%d', b);
end

% 3. Plot Markers
h_ref = plot(Params.Kp * Params.Ts, Params.y_ref_phys, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
h_opt = plot(Params.Kp * Params.Ts, y_star_k, 'gp', 'MarkerSize', 12, 'MarkerFaceColor', 'g');

% 4. Complete the legend handles & strings
h_leg = [h_leg, h_ref, h_opt, h_phy, h_tight];
leg_str{end+1} = 'Ref Target';
leg_str{end+1} = 'y*_s (Reachable)';
leg_str{end+1} = 'Physical Bounds';
leg_str{end+1} = 'Tightened Bounds';

xlabel('Time (min)', 'FontSize', 11);
ylabel('C_{ATR} (kg/kg)', 'FontSize', 11);
title('Tracking Trajectories with Constraints (Unreachable)', 'FontSize', 12);
xlim([0, Params.Ble * Params.Ts]);

% 图2图例设置为两列，放在右上角，字号增大到 10
legend(h_leg, leg_str, 'Location', 'northeast', 'NumColumns', 2, 'FontSize', 10);

% --- Inset Zoom (110-120 min) ---
% 将图中图放到左下角 [x, y, width, height]
ax2 = axes('Position', [0.15, 0.15, 0.35, 0.35]);
box on; hold on; grid on;
set(ax2, 'LineWidth', 1.2); % Thicker axes border for Zoom

% Plot Constraints in Zoom
plot([110, 120], [y_phys_max, y_phys_max], 'k--', 'LineWidth', 1.5);
plot([110, 120], [y_phys_min, y_phys_min], 'k--', 'LineWidth', 1.5);
plot([110, 120], [y_tight_max, y_tight_max], 'm-.', 'LineWidth', 1.5);
plot([110, 120], [y_tight_min, y_tight_min], 'm-.', 'LineWidth', 1.5);

% Plot Trajectories in Zoom
for idx = 1:length(batches_to_plot)
    b = batches_to_plot(idx);
    st = styles{idx};
    plot(t_y, squeeze(Data_prop.Y_nominal_hist(trk, :, b)), ['r', st], 'LineWidth', 1.5);
    plot(t_y, squeeze(Data_prop.Y_hist(trk, :, b)), ['b', st], 'LineWidth', 1.5);
end

% Plot Markers in Zoom
plot(Params.Kp * Params.Ts, Params.y_ref_phys, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(Params.Kp * Params.Ts, y_star_k, 'gp', 'MarkerSize', 12, 'MarkerFaceColor', 'g');

xlim([115, 120]);
ylim([0.013, 0.018]);
title('Zoom: 110-120 min', 'FontSize', 10);
set(ax2, 'FontSize', 8);

saveas(fig2, fullfile(out_dir, 'fig_unreachable_trajectory_zoom.png'));
try, print(fig2, fullfile(out_dir, 'fig_unreachable_trajectory_zoom.eps'), '-depsc2'); catch, end


%% 8. Active-constraint detection at nominal solutions
batch_list = [1, 3, 5, 10, 18, 19, 20];
batch_list = batch_list(batch_list <= B);

u_max_tight = Params.u_phys_max(Params.opt_input_idx) - Tube.delta_v(1:ni);
u_min_tight = Params.u_phys_min(Params.opt_input_idx) + Tube.delta_v(ni+1:2*ni);
y_max_tight = Params.y_phys_max - Tube.delta_y(1:no);
y_min_tight = Params.y_phys_min + Tube.delta_y(no+1:2*no);

tol_act_in  = 1e-2;
tol_act_out = 1e-6;

input_names  = Params.input_names(Params.opt_input_idx);
in_intervals = cell(1, ni);
for i = 1:ni
    in_intervals{i} = sprintf('[%.4f, %.4f]', u_min_tight(i), u_max_tight(i));
end
out_interval_str = sprintf('[%.6f, %.6f]', ...
                           y_min_tight(trk), y_max_tight(trk));

Mv = diag(Model.X_std) * Model.P;   % (ni x A_lv)

input_status  = cell(numel(batch_list), ni);
output_status = cell(numel(batch_list), 1);

for bi = 1:numel(batch_list)
    b   = batch_list(bi);
    t_b = Data_prop.T_hist(:,:,b);   % nominal latent input

    % Reconstruct physical nominal inputs
    V_nom = zeros(ni, Problem.H);
    for k = 1:Problem.H
        jj     = Params.ctrl_start + k;
        v_off  = Params.u_base(Model.input_idx, jj) + Model.X_mean;
        V_nom(:, k) = v_off + Mv * t_b(:, k);
    end

    for i = 1:ni
        if any(abs(V_nom(i,:) - u_max_tight(i)) <= tol_act_in) || ...
           any(abs(V_nom(i,:) - u_min_tight(i)) <= tol_act_in)
            input_status{bi, i} = 'Active';
        else
            input_status{bi, i} = 'Inactive';
        end
    end

    % Nominal C_ATR at tracked instants
    cols_track  = Params.Kp + 1;
    Y_nom_track = zeros(numel(cols_track), 1);
    for jj = 1:numel(cols_track)
        Y_nom_track(jj) = Data_prop.Y_nominal_hist(trk, cols_track(jj), b);
    end
    if any(abs(Y_nom_track - y_max_tight(trk)) <= tol_act_out) || ...
       any(abs(Y_nom_track - y_min_tight(trk)) <= tol_act_out)
        output_status{bi} = 'Active';
    else
        output_status{bi} = 'Inactive';
    end
end

%% 9. Print and Save Table
table_file = fullfile(out_dir, 'table_unreachable_active.csv');
fid = fopen(table_file, 'w');

% --- 1st Row: Table Header ---
fprintf(fid, 'Batch,y*_s(90),y*_s(120)');
for i = 1:ni
    fprintf(fid, ',%s', input_names{i});
end
fprintf(fid, ',C_ATR\n');

% --- 2nd Row: Tightened Intervals ---
fprintf(fid, 'Intervals,-,-');
for i = 1:ni
    fprintf(fid, ',"%s"', in_intervals{i}); % Add double quotes for CSV safety
end
fprintf(fid, ',"%s"\n', out_interval_str);

% --- 3rd Row onwards: Data ---
for bi = 1:numel(batch_list)
    b = batch_list(bi);
    fprintf(fid, '%d,%.8f,%.8f', b, y_star_k(1), y_star_k(2));
    for i = 1:ni
        fprintf(fid, ',%s', input_status{bi, i});
    end
    fprintf(fid, ',%s\n', output_status{bi});
end
fclose(fid);
fprintf('\nSaved table: %s\n', table_file);

% --- Print to Console (Formatted) ---
fprintf('\n=== Summary: Table ===\n');

% Header
fprintf('%-10s | %-10s | %-10s', 'Batch', 'y*_s(90)', 'y*_s(120)');
for i = 1:ni
    fprintf(' | %-28s', input_names{i});
end
fprintf(' | %-25s\n', 'C_ATR');

% Intervals (Second Row)
fprintf('%-10s | %-10s | %-10s', 'Intervals', '-', '-');
for i = 1:ni
    fprintf(' | %-28s', in_intervals{i});
end
fprintf(' | %-25s\n', out_interval_str);
fprintf(repmat('-', 1, 150)); fprintf('\n');

% Batches
for bi = 1:numel(batch_list)
    b = batch_list(bi);
    fprintf('%-10d | %-10.6f | %-10.6f', b, y_star_k(1), y_star_k(2));
    for i = 1:ni
        fprintf(' | %-28s', input_status{bi, i});
    end
    fprintf(' | %-25s\n', output_status{bi});
end

%% 10. Save full state
save('-mat7-binary', fullfile(out_dir, 'unreachable_example2_results.mat'), ...
     'Params','Model','Tube','UncSets','Problem','Initial','Optimal', ...
     'Data_prop','Err_nom_ki','Err_act_ki','bar_epsilon', ...
     'C_ATR_lo','C_ATR_hi','y_star_k', ...
     'batch_list','input_status','output_status', ...
     'u_max_tight','u_min_tight','y_max_tight','y_min_tight');

fprintf('\nAll outputs saved to: %s\n', out_dir);
