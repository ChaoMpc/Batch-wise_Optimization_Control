function report_simulation_comparison(Prop, Xue, Oh)
% =========================================================================
% report_simulation_comparison
% 严格对齐论文 Table 14 的极简格式输出（自动计算耗时）
% =========================================================================

    if nargin < 1, Prop = []; end
    if nargin < 2, Xue = []; end
    if nargin < 3, Oh = []; end

    % 1. 先计算好均值(每批时间)和总和(总时间)
    mean_prop = NaN; total_prop = NaN;
    mean_xue  = NaN; total_xue  = NaN;
    mean_oh   = NaN; total_oh   = NaN;

    if ~isempty(Prop) && isfield(Prop, 'batch_times')
        mean_prop = mean(Prop.batch_times);
        total_prop = sum(Prop.batch_times);
    end
    if ~isempty(Xue) && isfield(Xue, 'batch_times')
        mean_xue = mean(Xue.batch_times);
        total_xue = sum(Xue.batch_times);
    end
    if ~isempty(Oh) && isfield(Oh, 'batch_times')
        mean_oh = mean(Oh.batch_times);
        total_oh = sum(Oh.batch_times);
    end

    % 2. 打印极简表格
    fprintf('\n');
    fprintf('%-34s %16s %16s %16s\n', 'Metric', 'Proposed Method', 'Xue''s Method[30]', 'Oh''s Method[13]');
    fprintf(repmat('-', 1, 84)); fprintf('\n');

    print_row_3('Applied Input Violations', ...
        get_m(Prop,'input_applied_violations'), get_m(Xue,'input_applied_violations'), get_m(Oh,'input_applied_violations'));

    print_row_3('Input Boundary Hits', ...
        get_m(Prop,'input_boundary_hits'), get_m(Xue,'input_boundary_hits'), get_m(Oh,'input_boundary_hits'));

    cov_prop = NaN;
    if ~isempty(Prop) && isfield(Prop, 'tube_post_coverage') && ~isnan(Prop.tube_post_coverage)
        cov_prop = 100 * Prop.tube_post_coverage;
    end
    print_row_3('Tube State-Error Coverage (%)', cov_prop, NaN, NaN);

    print_row_3('Final Actual Tracking Error', ...
        get_m(Prop,'final_actual_reference_error'), get_m(Xue,'final_actual_reference_error'), get_m(Oh,'final_actual_reference_error'));

    print_row_3('Mean Solver Time (s/batch)', mean_prop, mean_xue, mean_oh);
    print_row_3('Total Simulation Time (s)', total_prop, total_xue, total_oh);

    fprintf('\n');
end

% 辅助函数
function val = get_m(S, field)
    if isempty(S) || ~isfield(S, field) || isempty(S.(field)) || any(isnan(S.(field)))
        val = NaN;
    else
        val = S.(field);
    end
end

function print_row_3(label, v1, v2, v3)
    if isnan(v1), v1_str = 'N/A'; else, v1_str = sprintf('%14.4g', v1); end
    if isnan(v2), v2_str = 'N/A'; else, v2_str = sprintf('%14.4g', v2); end
    if isnan(v3), v3_str = 'N/A'; else, v3_str = sprintf('%14.4g', v3); end
    fprintf('%-34s %16s %16s %16s\n', label, v1_str, v2_str, v3_str);
end
