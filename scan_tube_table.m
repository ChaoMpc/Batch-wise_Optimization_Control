function scan_tube_table(Ident_data, Model, UncSets, Params, opts)
% =========================================================================
% scan_tube_table -- Grid sweep over (alpha, beta) and print one combined table.
%
% 设计:
%   固定 alpha 的 3 个值 (opts.alpha_list), 外层循环 alpha, 内层遍历 beta.
%   每个 (alpha, beta) 计算:
%     tr(P_e), cond(P_e), delta_v, delta_y, 紧化物理区间, 未紧化物理区间,
%     盒可行性, 可行时的剩余优化空间比例 (v1, v2).
%
% 输出:
%   控制台逐条打印 + 一张 LaTeX 表 (第一列为 (alpha, beta)), 不同 alpha 组用
%   \midrule 分隔.
%
% 用法:
%   scan_tube_table(Ident_data, Model, UncSets, Params);
%   scan_tube_table(Ident_data, Model, UncSets, Params, opts);
%
% opts 字段 (均可选):
%   .alpha_list  : 默认 [0.03, 0.10, 0.20]
%   .beta_grid   : 默认 [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
%   .latex_dir   : 默认 'figures'
%   .tag         : 默认 ''
% =========================================================================

    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'alpha_list') || isempty(opts.alpha_list)
        opts.alpha_list = [0.03, 0.10, 0.20];
    end
    if ~isfield(opts,'beta_grid') || isempty(opts.beta_grid)
        opts.beta_grid = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0];
    end
    if ~isfield(opts,'latex_dir'), opts.latex_dir = 'figures'; end
    if ~isfield(opts,'tag'),       opts.tag = '';              end

    if exist('OCTAVE_VERSION','builtin'), pkg load control; end

    fprintf('\n=== scan_tube_table (alpha x beta grid) ===\n');

    % 通过 compute_tube 获得 K_bar / A_K / W_x / W_eta / W_r
    if exist('evalc','builtin')
        evalc('Tube0 = compute_tube(Ident_data, Model, UncSets, Params);');
    else
        Tube0 = compute_tube(Ident_data, Model, UncSets, Params);
    end

    ctx = build_ctx(Ident_data, Model, UncSets, Params, Tube0);
    fprintf('  rho_K = %.4f, alpha_max (stability) = %.4f\n', ...
        ctx.rho_K, ctx.alpha_max);

    % 过滤 alpha (稳定性)
    alpha_valid = opts.alpha_list(opts.alpha_list > 0 & ...
        (1 + opts.alpha_list) * ctx.rho_K^2 < 1 - 1e-6);
    removed = setdiff(opts.alpha_list, alpha_valid);
    if ~isempty(removed)
        fprintf('  [info] 以下 alpha 违反稳定性条件, 已跳过: %s\n', ...
            num2str(removed, '%.4f '));
    end
    if isempty(alpha_valid)
        error('scan_tube_table: 没有满足稳定性条件的 alpha 候选。');
    end

    n_alpha = numel(alpha_valid);
    n_beta  = numel(opts.beta_grid);
    fprintf('  Sweeping %d alpha x %d beta = %d combos\n', ...
        n_alpha, n_beta, n_alpha * n_beta);

    % ---- 组合扫描 ----
    rec = cell(n_alpha * n_beta, 1);
    idx = 0;
    for ai = 1:n_alpha
        for bi = 1:n_beta
            idx = idx + 1;
            r = eval_one(ctx, alpha_valid(ai), opts.beta_grid(bi));
            r.alpha_group = ai;
            rec{idx} = r;
        end
    end

    % ---- 控制台打印 ----
    print_records_console(rec, alpha_valid);

    % ---- LaTeX 表 ----
    latex_s = latex_table_string(rec, alpha_valid);
    fprintf('\n--- LaTeX Table ---\n%s\n', latex_s);

    if ~isempty(opts.latex_dir)
        if ~exist(opts.latex_dir,'dir'), mkdir(opts.latex_dir); end
        fn = fullfile(opts.latex_dir, sprintf('tab_sens_grid%s.tex', opts.tag));
        fid = fopen(fn,'w'); fprintf(fid,'%s\n',latex_s); fclose(fid);
        fprintf('  [Saved] %s\n', fn);
    end

    assignin('base','scan_grid', rec);
end

% =========================================================================
function ctx = build_ctx(Ident_data, Model, UncSets, Params, Tube0)
    ctx.A_K    = Tube0.A_K;
    ctx.K_bar  = Tube0.K;
    ctx.W_x    = Tube0.W_x;
    ctx.W_eta  = Tube0.W_eta;
    ctx.W_r    = Tube0.W_r;
    ctx.rho_K  = Tube0.spectral_radius;
    ctx.alpha_max = 1 / ctx.rho_K^2 - 1;

    ctx.A_lv  = Params.A_lv;
    ctx.n_in  = Params.n_in;
    ctx.n_out = Params.n_out;
    ctx.nx    = size(ctx.A_K, 1);

    ctx.B_bar = Model.barB;
    ctx.E_bar = Model.barE;
    ctx.C_bar = Model.barC;
    ctx.P_pls = Model.P;
    ctx.Q_pls = Model.Q;
    ctx.S_v   = diag(Model.X_std(:));
    ctx.S_y   = diag(Model.Y_std(:));

    ctx.H_v = Params.H_v;
    ctx.H_y = Params.H_y;
    ctx.h_v_phys = Params.h_v_phys(:);
    ctx.h_y_phys = Params.h_y_phys(:);

    ctx.u_ss = Params.u_ss(:);
    ctx.y_ss = Params.y_ss(:);
    ctx.u_phys_min = Params.u_phys_min(:);
    ctx.u_phys_max = Params.u_phys_max(:);
    ctx.y_phys_min = Params.y_phys_min(:);
    ctx.y_phys_max = Params.y_phys_max(:);
end

% =========================================================================
function r = eval_one(ctx, alpha, beta)
    r = struct();
    r.alpha = alpha;
    r.beta  = beta;
    r.ok_lyap = false;
    r.trace_Pe = NaN;
    r.cond_Pe  = NaN;
    r.delta_v  = [];
    r.delta_y  = [];
    r.v1_tight = [NaN NaN]; r.v2_tight = [NaN NaN];
    r.y1_tight = [NaN NaN]; r.y2_tight = [NaN NaN];
    r.feasible = false;
    r.ratio_v1 = NaN;
    r.ratio_v2 = NaN;

    nx = ctx.nx;
    BK = ctx.B_bar * ctx.K_bar;

    W_e = (1 + beta) * BK * ctx.W_x * BK' ...
        + (1 + 1/beta) * ctx.E_bar * ctx.W_eta * ctx.E_bar';
    W_e = 0.5 * (W_e + W_e');

    try
        As = sqrt(1 + alpha) * ctx.A_K;
        Qs = (1 + 1/alpha) * W_e + 1e-12 * eye(nx);
        P_e = dlyap(As, Qs);
        P_e = 0.5 * (P_e + P_e');
        [L, pflag] = chol(P_e, 'lower');
        if pflag ~= 0, return; end
        image_shape = (1 + alpha) * ctx.A_K * P_e * ctx.A_K' ...
                    + (1 + 1/alpha) * W_e;
        inv_ratio = max(eig(L \ image_shape / L'));
        if inv_ratio > 1 + 1e-7, return; end
    catch
        return;
    end
    r.ok_lyap = true;
    r.trace_Pe = trace(P_e);
    r.cond_Pe  = cond(P_e);

    M_v = ctx.S_v * ctx.P_pls * ctx.K_bar;
    n_constr_v = size(ctx.H_v, 1);
    delta_v = zeros(n_constr_v, 1);
    for i = 1:n_constr_v
        h_i = ctx.H_v(i, :);
        t1 = h_i * M_v * P_e * M_v' * h_i';
        t2 = h_i * M_v * ctx.W_x * M_v' * h_i';
        delta_v(i) = sqrt(max(t1, 0)) + sqrt(max(t2, 0));
    end
    r.delta_v = delta_v;

    C_y = ctx.S_y * ctx.Q_pls * ctx.C_bar;
    n_constr_y = size(ctx.H_y, 1);
    delta_y = zeros(n_constr_y, 1);
    for i = 1:n_constr_y
        h_i = ctx.H_y(i, :);
        t1 = h_i * C_y * P_e * C_y' * h_i';
        t2 = h_i * ctx.S_y * ctx.W_r * ctx.S_y' * h_i';
        delta_y(i) = sqrt(max(t1, 0)) + sqrt(max(t2, 0));
    end
    r.delta_y = delta_y;

    r.v1_tight = [ctx.u_phys_min(1) + delta_v(ctx.n_in + 1), ...
                  ctx.u_phys_max(1) - delta_v(1)];
    r.v2_tight = [ctx.u_phys_min(2) + delta_v(ctx.n_in + 2), ...
                  ctx.u_phys_max(2) - delta_v(2)];
    r.y1_tight = [ctx.y_phys_min(1) + delta_y(ctx.n_out + 1), ...
                  ctx.y_phys_max(1) - delta_y(1)];
    r.y2_tight = [ctx.y_phys_min(2) + delta_y(ctx.n_out + 2), ...
                  ctx.y_phys_max(2) - delta_y(2)];

    h_v_tight = ctx.h_v_phys - delta_v;
    h_y_tight = ctx.h_y_phys - delta_y;
    r.feasible = all(h_v_tight > 0) && all(h_y_tight > 0);

    if r.feasible
        av1 = ctx.u_phys_max(1) - ctx.u_phys_min(1);
        av2 = ctx.u_phys_max(2) - ctx.u_phys_min(2);
        r.ratio_v1 = (r.v1_tight(2) - r.v1_tight(1)) / av1;
        r.ratio_v2 = (r.v2_tight(2) - r.v2_tight(1)) / av2;
    end
end

% =========================================================================
function print_records_console(rec, alpha_valid)
    fprintf('\n%s\n', repmat('=', 1, 125));
    fprintf('  Row |  alpha   beta  |  tr(P_e)      cond(P_e)     | Feasible | Space(v1,v2)\n');
    fprintf('%s\n', repmat('-', 1, 125));
    for k = 1:numel(rec)
        r = rec{k};
        if ~r.ok_lyap
            fprintf('  %3d |  %.4f  %.2f  : Lyapunov/invariance FAILED\n', k, r.alpha, r.beta);
            continue;
        end
        fprintf('  %3d |  %.4f  %.2f  | %11.3e %11.3e | %-8s | %s\n', ...
            k, r.alpha, r.beta, r.trace_Pe, r.cond_Pe, ...
            ternary_str(r.feasible, 'YES', 'NO'), ratio_str(r.ratio_v1, r.ratio_v2));
        fprintf('       delta_v = [%s]  delta_y = [%s]\n', ...
            num2str(r.delta_v(:)', '%.3e '), num2str(r.delta_y(:)', '%.3e '));
        fprintf('       v1 [%.4f, %.4f]  v2 [%.4f, %.4f]  y1 [%.4f, %.4f]  y2 [%.4f, %.4f]\n', ...
            r.v1_tight(1), r.v1_tight(2), r.v2_tight(1), r.v2_tight(2), ...
            r.y1_tight(1), r.y1_tight(2), r.y2_tight(1), r.y2_tight(2));
    end
    fprintf('%s\n', repmat('=', 1, 125));
    fprintf('  alpha groups: %s\n', num2str(alpha_valid, '%.4f '));
end

% =========================================================================
function s = latex_table_string(rec, alpha_valid) %#ok<INUSD>
% 生成一张 LaTeX 表, 不同 alpha 组之间用 \midrule 分隔.
    lines = {};
    lines{end+1} = '\\begin{table}[htbp]';
    lines{end+1} = '    \\caption{Tube parameter sensitivity to $\\alpha$ and $\\beta$.}';
    lines{end+1} = '    \\label{tab:sens_grid}';
    lines{end+1} = '    \\centering';
    lines{end+1} = '    \\resizebox{\\textwidth}{!}{%';
    lines{end+1} = '    \\begin{tabular}{c|cc|cc|c|c|c|c}';
    lines{end+1} = '        \\toprule';
    lines{end+1} = ['        $(\\alpha, \\beta)$ & $\\mathrm{tr}(P_e)$ & $\\mathrm{cond}(P_e)$ & ' ...
        '$[\\delta_{v,1}, \\delta_{v,2}]$ & $[\\delta_{y,1}, \\delta_{y,2}]$ & ' ...
        '$[c_{v1}, d_{v1}] / [c_{v2}, d_{v2}]$ & $[c_{y1}, d_{y1}] / [c_{y2}, d_{y2}]$ & ' ...
        'Feasible & Space $(v_1, v_2)$ \\\\'];
    lines{end+1} = '        \\midrule';

    prev_group = -1;
    for k = 1:numel(rec)
        r = rec{k};
        if r.alpha_group ~= prev_group
            if prev_group > 0
                lines{end+1} = '        \\midrule'; %#ok<AGROW>
            end
            prev_group = r.alpha_group;
        end
        tag = sprintf('$(%.4f, %.2f)$', r.alpha, r.beta);
        if ~r.ok_lyap
            lines{end+1} = sprintf('        %s & \\multicolumn{8}{c}{Lyapunov/invariance failed} \\\\', tag); %#ok<AGROW>
            continue;
        end
        dv_str = sprintf('$[%.2e, %.2e]$', r.delta_v(1), r.delta_v(2));
        dy_str = sprintf('$[%.2e, %.2e]$', r.delta_y(1), r.delta_y(2));
        v_str  = sprintf('$[%.4f, %.4f]$ / $[%.4f, %.4f]$', ...
            r.v1_tight(1), r.v1_tight(2), r.v2_tight(1), r.v2_tight(2));
        y_str  = sprintf('$[%.4f, %.4f]$ / $[%.4f, %.4f]$', ...
            r.y1_tight(1), r.y1_tight(2), r.y2_tight(1), r.y2_tight(2));
        feas_str = ternary_str(r.feasible, 'yes', 'no');
        if r.feasible
            sp_str = sprintf('$(%.3f, %.3f)$', r.ratio_v1, r.ratio_v2);
        else
            sp_str = 'N/A';
        end
        lines{end+1} = sprintf('        %s & $%.3e$ & $%.3e$ & %s & %s & %s & %s & %s & %s \\\\', ... %#ok<AGROW>
            tag, r.trace_Pe, r.cond_Pe, dv_str, dy_str, v_str, y_str, feas_str, sp_str);
    end

    lines{end+1} = '        \\bottomrule';
    lines{end+1} = '    \\end{tabular}}';
    lines{end+1} = '    \\par\\smallskip';
    lines{end+1} = '    \\footnotesize Lower-bound margins satisfy $\\delta_{v,3}=\\delta_{v,1}$ and $\\delta_{v,4}=\\delta_{v,2}$ (and likewise for $\\delta_y$), due to the symmetric Lyapunov solution.';
    lines{end+1} = '\\end{table}';
    s = strjoin(lines, newline);
end

% =========================================================================
function s = ratio_str(r1, r2)
    if isnan(r1) || isnan(r2)
        s = 'N/A';
    else
        s = sprintf('(%.3f, %.3f)', r1, r2);
    end
end

function s = ternary_str(cond, a, b)
    if cond, s = a; else, s = b; end
end
