function report_simulation_comparison(Prop, Xue, Oh)
% =========================================================================
% report_simulation_comparison
%
% Print comparison without inferring unobserved causal conclusions.
% =========================================================================

    if nargin < 3, Oh = []; end

    fprintf('\n================ Simulation comparison ================\n');
    if isempty(Oh)
        fprintf('%-34s %14s %14s\n','Metric','Tube','Xue adapted');
    else
        fprintf('%-34s %14s %14s %14s\n','Metric','Tube','Xue adapted','Oh adapted');
    end

    % ==== 输入/输出约束违背统计 ====
    print_row_3('Applied input violations', ...
        Prop.input_applied_violations, Xue.input_applied_violations, get_m(Oh,'input_applied_violations'));

    print_row_3('Requested input violations', ...
        Prop.input_request_violations, Xue.input_request_violations, get_m(Oh,'input_request_violations'));

    print_row_3('Input clipping channel-samples', ...
        Prop.input_clips, Xue.input_clips, get_m(Oh,'input_clips'));

    print_row_3('Input boundary hits', ...
        Prop.input_boundary_hits, Xue.input_boundary_hits, get_m(Oh,'input_boundary_hits'));

    print_row_3('Measured output violations', ...
        Prop.output_measured_violations, Xue.output_measured_violations, get_m(Oh,'output_measured_violations'));

    print_row_3('Noise-free output violations', ...
        Prop.output_true_violations, Xue.output_true_violations, get_m(Oh,'output_true_violations'));

    print_row_3('Max measured output violation', ...
        Prop.max_output_measured_violation, Xue.max_output_measured_violation, get_m(Oh,'max_output_measured_violation'));

    % ==== 误差统计 ====
    print_row_3('Final actual / target error', ...
        Prop.final_actual_reference_error, Xue.final_actual_reference_error, get_m(Oh,'final_actual_reference_error'));

    print_row_3('Final actual / optimum error', ...
        Prop.final_actual_optimal_error, Xue.final_actual_optimal_error, get_m(Oh,'final_actual_optimal_error'));

    print_row_3('Final nominal / optimum error', ...
        Prop.final_nominal_optimal_error, Xue.final_nominal_optimal_error, get_m(Oh,'final_nominal_optimal_error'));

    print_row_3('Mean solver time (s)', ...
        Prop.mean_solver_time, Xue.mean_solver_time, get_m(Oh,'mean_solver_time'));

    % ==== 不确定集覆盖率 ====
    fprintf('\nOnline uncertainty coverage (joint vectors)\n');
    print_row_3('Innovation box coverage (%)', ...
        100*Prop.eta_box_coverage, 100*Xue.eta_box_coverage, get_m(Oh,'eta_box_coverage')*100);

    print_row_3('Residual box coverage (%)', ...
        100*Prop.r_box_coverage, 100*Xue.r_box_coverage, get_m(Oh,'r_box_coverage')*100);

    print_row_3('Innovation ellipsoid coverage (%)', ...
        100*Prop.eta_ellipsoid_coverage, 100*Xue.eta_ellipsoid_coverage, get_m(Oh,'eta_ellipsoid_coverage')*100);

    print_row_3('Residual ellipsoid coverage (%)', ...
        100*Prop.r_ellipsoid_coverage, 100*Xue.r_ellipsoid_coverage, get_m(Oh,'r_ellipsoid_coverage')*100);

    % ==== Tube 状态覆盖率 ====
    fprintf('\nTube-method state-error coverage\n');
    fprintf('  Before input: %.3f%%; after input: %.3f%%.\n', ...
        100*Prop.tube_pre_coverage,100*Prop.tube_post_coverage);

    fprintf('  Maximum post-input Tube quadratic: %.6g.\n', ...
        Prop.max_tube_post_quadratic);

    fprintf('========================================================\n');
end

function val = get_m(S, field)
    if isempty(S) || ~isfield(S, field) || isnan(S.(field))
        val = NaN;
    else
        val = S.(field);
    end
end

function print_row_3(label, v1, v2, v3)
    if isnan(v3)
        fprintf('%-34s %14.6g %14.6g %14s\n', label, v1, v2, 'N/A');
    else
        fprintf('%-34s %14.6g %14.6g %14.6g\n', label, v1, v2, v3);
    end
end
