% =========================================================================
% run_diagnostics_example2.m
% -------------------------------------------------------------------------
% One-click diagnostics: extract all tables and statistics requested by the reviewer.
% =========================================================================
clear; clc; close all;
if ~exist('OCTAVE_VERSION','builtin'), error('Octave is required'); end

proj  = fileparts(mfilename('fullpath'));
addpath(proj); addpath(fullfile(proj,'diagnostics'));
rehash();

out_dir    = fullfile(proj,'diagnostics','out');
table_dir  = fullfile(out_dir,'tables');
figure_dir = fullfile(out_dir,'figures');
for d = {out_dir, table_dir, figure_dir}
    if exist(d{1},'dir')~=7, mkdir(d{1}); end
end

% =========================================================================
% Start Diary: Capture all terminal output to a single text file
% =========================================================================
master_log = fullfile(out_dir, 'full_diagnostics_report.txt');
if exist(master_log, 'file'), delete(master_log); end
diary(master_log);

res_dir = fullfile(proj,'results');
S5 = load(fullfile(res_dir,'stage05_simulation.mat'));
S2 = load(fullfile(res_dir,'stage02_identification.mat'));
S3 = load(fullfile(res_dir,'stage03_tube.mat'));

Params  = S5.Params;  Model = S5.Model;  Tube = S5.Tube;
Optimal = S5.Optimal; Problem = S5.Problem;
UncSets = S3.UncSets; Prep = S2.Prep;
Validation = S2.Validation; Val_data = S2.Val_data;
Results = S5;  Results.UncSets = UncSets;

%% [A] PLS Evidence Table
report_pls_evidence_example2(Model, Prep, Validation, Val_data, Params);

%% [B] QP Condition Number
QPCond = compute_qp_conditioning(Problem, Model, Params);

%% [C] \epsilon_bar and Redrawn Figure 16
plot_figure16_redraw(Results, figure_dir);

%% [D] Active Constraints at the Nominal Optimum (Reviewer Major 7a)
fprintf('\n[Active Constraints at the Nominal Optimum (Problem 47)]\n');
% Compute slack = b - A*t
slack = Problem.b - Problem.A * Optimal.t;
tol_active = 1e-6;
active_idx = find(slack <= tol_active);

H  = Problem.H;
ni = Params.n_in_opt;
no = Params.n_out;

active_file = fullfile(table_dir, 'active_constraints.txt');
fid_active = fopen(active_file, 'w');
fprintf(fid_active, '=== Active Constraints at the Nominal Optimal Solution ===\n');
fprintf(fid_active, 'Total inequality constraints: %d\n', length(slack));
fprintf(fid_active, 'Number of active constraints (slack <= %g): %d\n\n', tol_active, length(active_idx));

if isempty(active_idx)
    fprintf(fid_active, 'No active constraints found (the solution is strictly interior).\n');
else
    for idx = 1:length(active_idx)
        row = active_idx(idx);
        val = Problem.A(row,:) * Optimal.t;
        bnd = Problem.b(row);

        if row <= 2*ni*H
            % Input constraint
            k = floor((row-1) / (2*ni)) + 1;
            rem = mod(row-1, 2*ni) + 1;
            if rem <= ni
                dir = 'UPPER'; i = rem;
            else
                dir = 'LOWER'; i = rem - ni;
            end
            time_abs = Params.ctrl_start + k;
            in_name = Params.input_names{Model.input_idx(i)};
            fprintf(fid_active, 'Input : %-6s | k=%2d (t=%3d min) | Dir: %s | Val: %10.4e | Bound: %10.4e\n', ...
                in_name, time_abs, time_abs*Params.Ts, dir, val, bnd);
        else
            % Output constraint
            rem_row = row - 2*ni*H;
            k = floor((rem_row-1) / (2*no)) + 1;
            rem = mod(rem_row-1, 2*no) + 1;
            if rem <= no
                dir = 'UPPER'; j = rem;
            else
                dir = 'LOWER'; j = rem - no;
            end
            time_abs = Params.ctrl_start + k + 1; % Output corresponds to time k+1
            out_name = Params.output_names{j};
            fprintf(fid_active, 'Output: %-6s | k=%2d (t=%3d min) | Dir: %s | Val: %10.4e | Bound: %10.4e\n', ...
                out_name, time_abs, time_abs*Params.Ts, dir, val, bnd);
        end
    end
end
fclose(fid_active);
fprintf('  Active constraints saved to: %s\n', active_file);
if ~isempty(active_idx)
    fprintf('  Found %d active constraints (see table for details).\n', length(active_idx));
else
    fprintf('  No active constraints (interior solution).\n');
end

%% [E] Reference vs. Optimal Reachable Points
fprintf('\n[Reference Points vs. Optimal Reachable Points]\n');
Kp = Params.Kp(:);
y_ref = Params.y_ref_phys(:);
y_opt = Optimal.selected_outputs(:);

reach_file = fullfile(table_dir, 'table_reachable_points.txt');
fid_reach = fopen(reach_file, 'w');
fprintf(fid_reach, '=== Reference Points vs. Optimal Reachable Points ===\n');
fprintf(fid_reach, '%-10s | %-18s | %-18s | %-18s\n', 'Time (min)', 'Reference Value', 'Reachable Value', 'Difference');
fprintf(fid_reach, '-------------------------------------------------------------------------\n');

for j = 1:numel(Kp)
    t_min = Kp(j) * Params.Ts;
    diff_val = y_ref(j) - y_opt(j);
    fprintf(fid_reach, '%-10.1f | %-18.8e | %-18.8e | %+.4e\n', ...
        t_min, y_ref(j), y_opt(j), diff_val);

    fprintf('  k_i=%2d (t=%3.0f min): Ref = %.6e, Opt = %.6e, Diff = %+.2e\n', ...
        Kp(j), t_min, y_ref(j), y_opt(j), diff_val);
end
fclose(fid_reach);
fprintf('  Comparison saved to: %s\n', reach_file);

%% [F] Generate Parameter Summary Table (for Reviewer Major 5e)
params_file = fullfile(table_dir, 'table_params.txt');
fid = fopen(params_file, 'w');
fprintf(fid, '=== Simulation Parameters Summary ===\n');
fprintf(fid, 'Ts              : %.2f min\n', Params.Ts);
fprintf(fid, 'Ble             : %d\n', Params.Ble);
fprintf(fid, 'ctrl_start      : %d\n', Params.ctrl_start);
fprintf(fid, 'H (Horizon)     : %d\n', Problem.H);
fprintf(fid, 'N_train_batches : %d\n', Params.N_train_batches);
fprintf(fid, 'A_lv            : %d\n', Params.A_lv);
fprintf(fid, 'n_arx           : %d\n', Params.n_arx);
fprintf(fid, 'lambda_a (ARX)  : %.1e\n', Params.arx_lambda);
fprintf(fid, 'Q_weight        : %.1e\n', Params.Q_weight);
fprintf(fid, 'R_weight        : %.1e\n', Params.R_weight);
fprintf(fid, 'P_weight        : %.1e\n', Params.P_weight);
fprintf(fid, 'delta_t_max     : %.2f\n', Params.delta_t_max);
fprintf(fid, 'ddt_max         : %.2f\n', Params.ddt_max);
fprintf(fid, 'tube_beta       : %.2f\n', Params.tube_beta);
fprintf(fid, 'k_anti          : %.2f\n', Params.k_anti);
fprintf(fid, 'UA0             : %.2f\n', Params.UA0);
fprintf(fid, 'F_c0            : %.2f\n', Params.F_c0);
fprintf(fid, 'N0              : %.2f\n', Params.N0);
fprintf(fid, 'C_feed          : %.4f\n', Params.C_feed);
fprintf(fid, 'k_seed          : %.4f\n', Params.k_seed);
fprintf(fid, '\n[Input Handling]\n');
fprintf(fid, 'opt_input_idx   : [%s]\n', num2str(Params.opt_input_idx'));
fprintf(fid, 'F_seed is excluded from opt_input_idx (Channel 6), hence Sv^-1 for F_seed is bypassed.\n');
fclose(fid);

%% [G] Summary Output (Appended to diary)
fprintf('\n=== Example 2 diagnostics ===\n');
fprintf('A_lv = %d, n_arx = %d\n', Model.A_lv, Model.n_arx);
fprintf('nx   = %d, H = %d\n', Model.nx, Problem.H);
fprintf('sigma_Q = %.4e, cond(Q) = %.4e\n', Model.PLS.sigma_min_Q, cond(Model.PLS.Q));
fprintf('cond(H) latent = %.4e, raw = %.4e\n', QPCond.lat_cond, QPCond.raw_cond);
fprintf('\n--- Reference vs Optimal Reachable ---\n');
for j = 1:numel(Kp)
    t_min = Kp(j) * Params.Ts;
    fprintf('Time %.1f min: Ref = %.6e, Opt = %.6e\n', t_min, y_ref(j), y_opt(j));
end

%% [H] Append Table Contents to the Master Log
fprintf('\n\n==============================================\n');
fprintf('APPENDIX: Content of generated tables\n');
fprintf('==============================================\n');

files_to_append = {active_file, reach_file, params_file};
for f = 1:numel(files_to_append)
    file_path = files_to_append{f};
    if exist(file_path, 'file')
        [~, name, ext] = fileparts(file_path);
        fprintf('\n--- Start of %s%s ---\n', name, ext);
        fid_tmp = fopen(file_path, 'r');
        while ~feof(fid_tmp)
            line = fgetl(fid_tmp);
            if ischar(line)
                fprintf('%s\n', line);
            end
        end
        fclose(fid_tmp);
        fprintf('--- End of %s%s ---\n', name, ext);
    end
end

fprintf('\nAll diagnostics completed. The complete report has been saved to:\n  %s\n', master_log);

% Close diary
diary off;
