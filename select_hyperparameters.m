function [Model, CV] = select_hyperparameters( ...
    Prep, Standard_train, Standard_val, Val_data, Params)
% =========================================================================
% select_hyperparameters
%
% 打印格式升级：同时显示输入和输出的平均 VAF 和 NRMSE
% =========================================================================

    A_range = Params.A_lv_cv_range;
    n_range = Params.n_arx_cv_range;

    entries = struct( ...
        'A_requested',{}, 'A_effective',{}, 'n_arx',{}, ...
        'score',{}, 'NRMSE',{}, 'VAF',{}, ...
        'in_NRMSE',{}, 'in_VAF',{}, ...
        'n_parameters',{}, 'rho_open',{}, 'elapsed',{}, ...
        'success',{}, 'message',{});

    best_score = Inf;
    best_parameters = Inf;
    best_model = [];
    best_index = 0;

    index = 0;

    fprintf('\n================ Cross-Validation ================\n');
    fprintf([' A_req A_eff   n   Out NRMSE  Out VAF(%%)  ', ...
             ' In NRMSE  In VAF(%%)  rho(A)   status\n']);
    fprintf('--------------------------------------------------\n');

    for ia = 1:numel(A_range)
        for in = 1:numel(n_range)
            index = index+1;
            timer = tic();

            entry.A_requested = A_range(ia);
            entry.A_effective = NaN;
            entry.n_arx = n_range(in);
            entry.score = Inf;
            entry.NRMSE = NaN(1,Params.n_out);
            entry.VAF = NaN(1,Params.n_out);
            entry.in_NRMSE = NaN(1,Params.n_in_opt);
            entry.in_VAF = NaN(1,Params.n_in_opt);
            entry.n_parameters = Inf;
            entry.rho_open = NaN;
            entry.elapsed = NaN;
            entry.success = false;
            entry.message = '';

            try
                Candidate = build_identified_model( ...
                    Prep,Standard_train, ...
                    A_range(ia),n_range(in),Params);

                Metrics = validate_model( ...
                    Candidate,Standard_val,Val_data,Params);

                % --- 提取输出指标 ---
                out_nrmse_mean = mean(Metrics.NRMSE);
                out_vaf_mean   = mean(Metrics.VAF);

                % --- 提取输入指标 ---
                in_nrmse_mean  = mean(Metrics.input_NRMSE);
                in_vaf_mean    = mean(Metrics.input_VAF);

                n_parameters = Candidate.ARX.n_parameters;

                entry.A_effective = Candidate.A_lv;
                entry.score = out_nrmse_mean; % 仍然用输出 NRMSE 选优
                entry.NRMSE = Metrics.NRMSE;
                entry.VAF = Metrics.VAF;
                entry.in_NRMSE = Metrics.input_NRMSE;
                entry.in_VAF = Metrics.input_VAF;
                entry.n_parameters = n_parameters;
                entry.rho_open = Candidate.open_loop_spectral_radius;
                entry.success = true;

                if ~isfinite(entry.score)
                    error('Nonfinite cross-validation score.');
                end

                is_better = entry.score < best_score-1e-12;
                is_tie_better = ...
                    abs(entry.score-best_score) <= 1e-12 && ...
                    n_parameters < best_parameters;

                if isempty(best_model) || is_better || is_tie_better
                    best_score = entry.score;
                    best_parameters = n_parameters;
                    best_model = Candidate;
                    best_index = index;
                end

                % 打印这一行，格式齐整
                fprintf('%6d %5d %3d %11.5f %10.2f %10.5f %10.2f %9.5f   OK\n', ...
                    entry.A_requested,entry.A_effective,entry.n_arx, ...
                    out_nrmse_mean, out_vaf_mean, in_nrmse_mean, in_vaf_mean, ...
                    entry.rho_open);

            catch ME
                entry.success = false;
                entry.message = ME.message;

                fprintf('%6d %5s %3d %11s %10s %10s %10s %9s   FAILED\n', ...
                    entry.A_requested,'-',entry.n_arx,'-','-','-','-','-');

                fprintf('  Reason: %s\n',ME.message);
            end

            entry.elapsed = toc(timer);
            entries(index) = entry;
        end
    end

    if isempty(best_model)
        error('select_hyperparameters: every candidate failed.');
    end

    Model = best_model;

    CV = struct();
    CV.entries = entries;
    CV.best_index = best_index;
    CV.best_score = best_score;
    CV.criterion = 'mean control-segment free-run output NRMSE';

    fprintf('\nSelected model: A_requested=%d, A_effective=%d, n=%d.\n', ...
        Model.A_requested,Model.A_lv,Model.n_arx);

    best_entry = entries(best_index);
    fprintf('Selected Output NRMSE: %.5f, Output VAF: %.2f%%\n', ...
        mean(best_entry.NRMSE), mean(best_entry.VAF));
    fprintf('Selected Input  NRMSE: %.5f, Input  VAF: %.2f%%\n', ...
        mean(best_entry.in_NRMSE), mean(best_entry.in_VAF));
    fprintf('==================================================\n');
end
