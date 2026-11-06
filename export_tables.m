function export_tables(Results, out_dir)
% =========================================================================
% export_tables
% 自动生成并导出 6 张核心数据表 (CSV 和 TXT)
% =========================================================================

    if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end

    txt_file = fullfile(out_dir, 'Metrics_Tables.txt');
    fid = fopen(txt_file, 'w');
    clean_txt = onCleanup(@() fclose(fid));

    fprintf(fid, '========================================================\n');
    fprintf(fid, '                PAPER METRICS AND TABLES\n');
    fprintf(fid, '========================================================\n\n');

    %% 表 1：潜变量交叉验证数据表
    fprintf(fid, 'Table 1: Cross-Validation Results for PLS-ARX\n');
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, '%-6s %-5s | %-12s %-12s | %-12s %-12s\n', ...
        'A_req', 'n', 'Out_NRMSE', 'Out_VAF(%)', 'In_NRMSE', 'In_VAF(%)');
    fprintf(fid, '--------------------------------------------------------\n');
    if isfield(Results, 'CV') && isfield(Results.CV, 'entries')
        entries = Results.CV.entries;
        for i = 1:numel(entries)
            if entries(i).success
                fprintf(fid, '%-6d %-5d | %-12.5f %-12.2f | %-12.5f %-12.2f\n', ...
                    entries(i).A_requested, entries(i).n_arx, ...
                    mean(entries(i).NRMSE), mean(entries(i).VAF), ...
                    mean(entries(i).in_NRMSE), mean(entries(i).in_VAF));
            end
        end
    end
    fprintf(fid, '\n');

    %% 表 2：不确定估计得到的椭圆矩阵参数表
    fprintf(fid, 'Table 2: Uncertainty Sets and Invariant Ellipsoid Parameters\n');
    fprintf(fid, '--------------------------------------------------------\n');
    Unc = Results.UncSets; Inv = Results.Invariant;
    fprintf(fid, 'Innovation Set (eta) trace(W_eta) : %.4e\n', trace(Unc.W_eta));
    fprintf(fid, 'Reconstruction Set (r) trace(W_r) : %.4e\n', trace(Unc.W_r));
    fprintf(fid, 'Initial Error Set trace(W_0)       : %.4e\n', trace(Results.StateError.W_initial));
    fprintf(fid, 'Online Est Error trace(W_eps)      : %.4e\n', trace(Results.StateError.W_epsilon));
    fprintf(fid, 'Disturbance Aggregation trace(W_e): %.4e (beta=%.2f)\n', trace(Inv.W_d), Inv.tube_beta);
    fprintf(fid, 'Invariant Ellipsoid trace(P_e)     : %.4e (alpha=%.2f)\n', trace(Inv.P), Inv.alpha);
    fprintf(fid, '\n');

    %% 表 3：RPI 外估计椭圆矩阵参数敏感度分析 (alpha 和 beta)
    fprintf(fid, 'Table 3: Sensitivity of Ellipsoid Volumes to alpha and beta\n');
    fprintf(fid, '--------------------------------------------------------\n');
    if isfield(Inv, 'Sens')
        S = Inv.Sens;
        % 打印 3(a)：不同 beta 对扰动集 W_e 的影响
        fprintf(fid, 'Table 3a: Effect of beta on Disturbance Set W_e (trace)\n');
        fprintf(fid, '%-8s', 'alpha\beta');
        for ib = 1:numel(S.beta_candidates), fprintf(fid, '| %-8.2f ', S.beta_candidates(ib)); end
        fprintf(fid, '\n');
        for ia = 1:numel(S.alpha_candidates)
            fprintf(fid, '%-8.2f', S.alpha_candidates(ia));
            for ib = 1:numel(S.beta_candidates)
                fprintf(fid, '| %-9.3e', S.trace_We(ia, ib));
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n');

        % 打印 3(b)：不同 alpha 和 beta 对最终包络 P_e 的影响
        fprintf(fid, 'Table 3b: Effect of alpha and beta on Invariant Set P_e (trace)\n');
        fprintf(fid, '%-8s', 'alpha\beta');
        for ib = 1:numel(S.beta_candidates), fprintf(fid, '| %-8.2f ', S.beta_candidates(ib)); end
        fprintf(fid, '\n');
        for ia = 1:numel(S.alpha_candidates)
            fprintf(fid, '%-8.2f', S.alpha_candidates(ia));
            for ib = 1:numel(S.beta_candidates)
                fprintf(fid, '| %-9.3e', S.trace_Pe(ia, ib));
            end
            fprintf(fid, '\n');
        end
    else
        fprintf(fid, ' (Sensitivity data not found. Please rerun Stage 3)\n');
    end
    fprintf(fid, '\n');

    %% 表 4：优化问题的条件数与降维对比
    fprintf(fid, 'Table 4: Optimization Problem Conditioning and Dimensionality Reduction\n');
    fprintf(fid, '--------------------------------------------------------\n');
    Prob = Results.Problem;
    H_opt = 2*(Prob.L'*Prob.L + Prob.Pprox);
    cond_H = cond(H_opt);
    var_before = Prob.H * Prob.n_in;
    var_after  = Prob.nt;
    fprintf(fid, 'Decision Variables (Before/After PLS) : %d -> %d\n', var_before, var_after);
    fprintf(fid, 'Reduction Percentage                  : %.1f%%\n', (1 - var_after/var_before)*100);
    fprintf(fid, 'Hessian Condition Number cond(H)      : %.4e\n', cond_H);
    fprintf(fid, '\n');

    %% 表 5：方法与 Xue 方法的对比
    fprintf(fid, 'Table 5: Performance Comparison (Proposed vs Xue)\n');
    fprintf(fid, '--------------------------------------------------------\n');
    M_p = Results.Metrics_prop; M_x = Results.Metrics_xue;
    fprintf(fid, '%-28s | %-15s | %-15s\n', 'Metric', 'Proposed Tube', 'Xue Method');
    fprintf(fid, '--------------------------------------------------------\n');
    fprintf(fid, 'Closed-loop Coverage (%%)     | %-15.2f | %-15.2f\n', 100*M_p.r_ellipsoid_coverage, 100*M_x.r_ellipsoid_coverage);
    fprintf(fid, 'Input Violations             | %-15d | %-15d\n', M_p.input_applied_violations, M_x.input_applied_violations);
    fprintf(fid, 'Output Violations            | %-15d | %-15d\n', M_p.output_measured_violations, M_x.output_measured_violations);
    fprintf(fid, 'Actuator Saturation (clips)  | %-15d | %-15d\n', M_p.input_clips, M_x.input_clips);
    fprintf(fid, 'Final Error (vs True Ref)    | %-15.4e | %-15.4e\n', M_p.final_actual_reference_error, M_x.final_actual_reference_error);
    fprintf(fid, 'Final Error (vs Opt Reach)   | %-15.4e | %-15.4e\n', M_p.final_actual_optimal_error, M_x.final_actual_optimal_error);
    fprintf(fid, '\n');

    %% 表 6：实际参考点与最优可达点
    fprintf(fid, 'Table 6: True Reference vs Optimal Reachable Points\n');
    fprintf(fid, '--------------------------------------------------------\n');
    if isfield(Results, 'Optimal')
        fprintf(fid, '%-15s | %-15s | %-15s | %-15s\n', 'Tracking Time', 'True Ref', 'Opt Reachable', 'Difference');
        for j = 1:numel(Results.Optimal.selected_references)
            t_k = Results.Params.Kp(j) * Results.Params.Ts;
            ref_val = Results.Optimal.selected_references(j);
            opt_val = Results.Optimal.selected_outputs(j);
            fprintf(fid, 't = %-11.1f | %-15.6f | %-15.6f | %-15.6f\n', ...
                t_k, ref_val, opt_val, ref_val - opt_val);
        end
    end
    fprintf(fid, '========================================================\n');

    fprintf('\nResults Tables (1-6) exported successfully to:\n  %s\n', txt_file);
end
