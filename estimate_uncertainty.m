function UncSets = estimate_uncertainty( ...
    Calibration_residuals, Training_residuals, ...
    Selection_residuals, Model, Params)
% =========================================================================
% 使用最经典的维度乘子法：将有界 Box 包络为椭球时，乘以维度 d
% =========================================================================

    Eta = Calibration_residuals.eta;
    R = Calibration_residuals.r;

    q = Params.q_level;
    safety = Params.safety_factor;

    % 分位数计算
    eta_bound = safety * column_quantile_local(abs(Eta), q)';
    r_bound = safety * column_quantile_local(abs(R), q)';

    % 防止奇异
    eta_bound = max(eta_bound, 1e-12);
    r_bound = max(r_bound, 1e-12);

    A = Model.A_lv;
    no = size(Model.Q,1);

    UncSets = struct();

    UncSets.state_convention = 'observed_score_history';
    UncSets.fit_source = 'independent calibration dataset';

    UncSets.q_level = q;
    UncSets.safety_factor = safety;

    UncSets.eta_bound = eta_bound;
    UncSets.r_bound = r_bound;

    % ====== 退回：乘以维度 A 和 no ======
    UncSets.W_eta = A * diag(eta_bound.^2);
    UncSets.W_r = no * diag(r_bound.^2);

    % 诊断
    UncSets.calibration = coverage_local(Calibration_residuals, eta_bound, r_bound);
    UncSets.training = coverage_local(Training_residuals, eta_bound, r_bound);
    UncSets.selection = coverage_local(Selection_residuals, eta_bound, r_bound);

    fprintf('\nUncertainty bounds fitted on calibration data (Original Box-to-Ellipsoid)\n');
    fprintf('  q_level=%.4f; safety_factor=%.4f.\n',q,safety);
    fprintf('  eta bound: [%s]\n',num2str(eta_bound','%.5e '));
    fprintf('  r bound:   [%s]\n',num2str(r_bound','%.5e '));

    print_coverage_local('Training',UncSets.training);
    print_coverage_local('Selection',UncSets.selection);
    print_coverage_local('Calibration',UncSets.calibration);
end

function qvalue = column_quantile_local(X,p)
    X = sort(X,1);
    count = size(X,1);
    position = 1+p*(count-1);
    lower = floor(position);
    upper = ceil(position);
    weight = position-lower;
    qvalue = (1-weight)*X(lower,:) + weight*X(upper,:);
end

function Report = coverage_local(Residuals,eta_bound,r_bound)
    eta = Residuals.eta;
    r = Residuals.r;
    eta_ratio = bsxfun(@rdivide,abs(eta),eta_bound');
    r_ratio = bsxfun(@rdivide,abs(r),r_bound');
    Report.n_samples = size(eta,1);
    Report.eta_box_coverage = mean(all(eta_ratio <= 1+1e-10,2));
    Report.r_box_coverage = mean(all(r_ratio <= 1+1e-10,2));
    eta_quad = sum(eta_ratio.^2,2)/numel(eta_bound);
    r_quad = sum(r_ratio.^2,2)/numel(r_bound);
    Report.eta_ellipsoid_coverage = mean(eta_quad <= 1+1e-10);
    Report.r_ellipsoid_coverage = mean(r_quad <= 1+1e-10);
end

function print_coverage_local(label,Report)
    fprintf(['  %-11s: eta box=%6.2f%%, r box=%6.2f%%; ', ...
             'eta ellipsoid=%6.2f%%, r ellipsoid=%6.2f%%.\n'], ...
        label, 100*Report.eta_box_coverage, 100*Report.r_box_coverage, ...
        100*Report.eta_ellipsoid_coverage, 100*Report.r_ellipsoid_coverage);
end
