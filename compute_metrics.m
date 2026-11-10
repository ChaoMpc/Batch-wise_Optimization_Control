function Metrics = compute_metrics(Data, Params)
% =========================================================================
% compute_metrics
%
% 所有 violation 计数按 channel-samples。
% 覆盖率按每个控制步的联合向量判断。
% =========================================================================

    B   = Params.N_batches;
    idx = Params.opt_input_idx;

    tolerance = 1e-8;

    applied_count   = 0;
    requested_count = 0;

    measured_output_count = 0;
    true_output_count     = 0;

    max_applied         = 0;
    max_requested       = 0;
    max_measured_output = 0;
    max_true_output     = 0;

    boundary_hits = 0;
    clips         = 0;

    eta_box = [];
    r_box   = [];
    eta_ellipsoid = [];
    r_ellipsoid   = [];

    tube_pre      = [];
    tube_post     = [];
    tube_pre_obs  = [];
    tube_post_obs = [];

    for b = 1:B
        V        = Data.V_hist(idx,:,b);
        Vrequest = Data.V_request_hist(idx,:,b);

        Y     = Data.Y_hist(:,:,b);
        Ytrue = Data.Y_true_hist(:,:,b);

        vlow  = Params.u_phys_min(idx);
        vhigh = Params.u_phys_max(idx);

        applied_violation = max( ...
            bsxfun(@minus,V,vhigh), ...
            bsxfun(@minus,vlow,V));

        requested_violation = max( ...
            bsxfun(@minus,Vrequest,vhigh), ...
            bsxfun(@minus,vlow,Vrequest));

        measured_violation = max( ...
            bsxfun(@minus,Y,Params.y_phys_max), ...
            bsxfun(@minus,Params.y_phys_min,Y));

        true_violation = max( ...
            bsxfun(@minus,Ytrue,Params.y_phys_max), ...
            bsxfun(@minus,Params.y_phys_min,Ytrue));

        applied_count   = applied_count   + sum(applied_violation(:)   > tolerance);
        requested_count = requested_count + sum(requested_violation(:) > tolerance);

        measured_output_count = measured_output_count + sum(measured_violation(:) > tolerance);
        true_output_count     = true_output_count     + sum(true_violation(:)     > tolerance);

        max_applied         = max([max_applied;         applied_violation(:);   0]);
        max_requested       = max([max_requested;       requested_violation(:); 0]);
        max_measured_output = max([max_measured_output; measured_violation(:);  0]);
        max_true_output     = max([max_true_output;     true_violation(:);      0]);

        at_upper = abs(bsxfun(@minus,V,vhigh)) < 1e-7;
        at_lower = abs(bsxfun(@minus,V,vlow))  < 1e-7;

        boundary_hits = boundary_hits + sum(at_upper(:) | at_lower(:));

        Batch = Data.Batch{b};

        clips = clips + sum(Batch.clip_flags(:));

        eta_box       = [eta_box,       Batch.eta_box_in];
        r_box         = [r_box,         Batch.r_box_in];
        eta_ellipsoid = [eta_ellipsoid, Batch.eta_ellipsoid_in];
        r_ellipsoid   = [r_ellipsoid,   Batch.r_ellipsoid_in];

        tube_pre      = [tube_pre,      Batch.tube_pre_quadratic];
        tube_post     = [tube_post,     Batch.tube_post_quadratic];

        if isfield(Batch,'tube_pre_obs_quadratic')
            tube_pre_obs  = [tube_pre_obs,  Batch.tube_pre_obs_quadratic];
            tube_post_obs = [tube_post_obs, Batch.tube_post_obs_quadratic];
        end
    end

    Metrics = struct();

    Metrics.method = Data.method;

    Metrics.input_applied_violations = applied_count;
    Metrics.input_request_violations = requested_count;

    Metrics.output_measured_violations = measured_output_count;
    Metrics.output_true_violations     = true_output_count;

    Metrics.max_input_applied_violation   = max_applied;
    Metrics.max_input_request_violation   = max_requested;
    Metrics.max_output_measured_violation = max_measured_output;
    Metrics.max_output_true_violation     = max_true_output;

    Metrics.input_boundary_hits = boundary_hits;
    Metrics.input_clips         = clips;

    Metrics.input_samples  = numel(idx)*Params.Ble*B;
    Metrics.output_samples = Params.n_out*(Params.Ble+1)*B;

    Metrics.eta_box_coverage = mean(eta_box);
    Metrics.r_box_coverage   = mean(r_box);

    Metrics.eta_ellipsoid_coverage = mean(eta_ellipsoid);
    Metrics.r_ellipsoid_coverage   = mean(r_ellipsoid);

    Metrics.tube_pre_coverage      = mean(tube_pre  <= 1+1e-7);
    Metrics.tube_post_coverage     = mean(tube_post <= 1+1e-7);
    Metrics.max_tube_pre_quadratic  = max(tube_pre);
    Metrics.max_tube_post_quadratic = max(tube_post);

    if isempty(tube_pre_obs)
        Metrics.tube_pre_obs_coverage  = NaN;
        Metrics.tube_post_obs_coverage = NaN;
        Metrics.max_tube_pre_obs_quadratic  = NaN;
        Metrics.max_tube_post_obs_quadratic = NaN;
    else
        Metrics.tube_pre_obs_coverage  = mean(tube_pre_obs  <= 1+1e-7);
        Metrics.tube_post_obs_coverage = mean(tube_post_obs <= 1+1e-7);
        Metrics.max_tube_pre_obs_quadratic  = max(tube_pre_obs);
        Metrics.max_tube_post_obs_quadratic = max(tube_post_obs);
    end

    Metrics.final_actual_reference_error  = Data.Err_true_act(end);
    Metrics.final_actual_optimal_error    = Data.Err_opt_act(end);
    Metrics.final_nominal_optimal_error   = Data.Err_opt_nom(end);

    % 注意：耗时字段已经从本函数中移除，交由主脚本提取并计算
end
