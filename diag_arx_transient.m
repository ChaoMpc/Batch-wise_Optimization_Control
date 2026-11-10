function diag_arx_transient(Ident_data, Val_data, Params)
% Diagnostic: how much transient growth is in A_a vs A_K, and how
% sensitive it is to arx_lambda. This is READ-ONLY; no Params modification.

A_lv   = Params.A_lv;
n_arx  = Params.n_arx;
n_state = 2*n_arx - 1;
lambda_list = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];

fprintf('\n================ ARX transient-growth diagnostics ================\n');
fprintf('%8s | %4s | %6s | %10s | %10s | %8s | %s\n', ...
    'lambda', 'LV', 'rho(A)', '||A||_2', 'max_j||A^j||', 'j_peak', 'free-run NRMSE');
fprintf('------------------------------------------------------------------\n');

for li = 1:numel(lambda_list)
    lam = lambda_list(li);
    Params_try = Params;
    Params_try.arx_lambda = lam;
    Params_try.plot_model = false;

    try
        if exist('evalc','builtin')
            evalc('Model_try = build_dypls(Ident_data, Params_try);');
        else
            Model_try = build_dypls(Ident_data, Params_try);
        end
    catch ME
        fprintf('%8.0e | (build failed: %s)\n', lam, ME.message);
        continue;
    end

    for a = 1:A_lv
        idx = (a-1)*n_state+1 : a*n_state;
        A_i = Model_try.barA(idx, idx);

        % Power-norm curve for j = 1..30
        J_show = 30;
        A_pow  = eye(n_state);
        norms  = zeros(J_show, 1);
        for j = 1:J_show
            A_pow = A_pow * A_i;
            norms(j) = norm(A_pow);
        end
        [maxn, jpeak] = max(norms);

        % Free-run NRMSE (whole-val)
        [~, nrmses, ~, ~, ~] = free_run_validation(Model_try, Val_data, Params_try);
        nrmse_str = num2str(nrmses, '%.3f ');

        fprintf('%8.0e | %4d | %6.3f | %10.3e | %10.3e | %8d | [%s]\n', ...
            lam, a, max(abs(eig(A_i))), norm(A_i), maxn, jpeak, nrmse_str);
    end
end

fprintf('==================================================================\n\n');
end
