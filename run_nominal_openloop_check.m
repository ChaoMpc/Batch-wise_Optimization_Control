function R = run_nominal_openloop_check(Results)
% 仅用标称物理输入（无 tube 反馈）驱动非线性晶化模型，观察
% 参考点处的实际 y4。用于回答"不可达是否来自 PLS 输入截断"。

    Params = Results.Params;
    Optimal = Results.Optimal;
    Kp = Params.Kp(:);
    trk = Params.track_idx;
    N = Params.Ble;

    V_nom = Params.u_base;
    V_nom(Params.opt_input_idx, Optimal.Trajectory.input_columns) = ...
        Optimal.Trajectory.V;

    x = Params.x0_nom;
    y = zeros(Params.n_out, N+1);
    [y(:,1)] = cryst_measure(x, Params, false);

    for j = 1:N
        [~, x, ~] = cryst_plant(x, V_nom(:,j), Params, false);
        y(:, j+1) = cryst_measure(x, Params, false);
    end

    R.V_nom    = V_nom;
    R.y        = y;
    R.y_at_kp  = y(trk, Kp+1);
    R.y_ref    = Params.y_ref_phys(:);
    R.residual = R.y_at_kp(:) - R.y_ref;

    fprintf('\n[Open-loop nominal check, no tube feedback]\n');
    fprintf('  k_i = %d : y_actual=%.6e  y_ref=%.6e  residual=%+.3e\n', ...
        [Kp, R.y_at_kp(:), R.y_ref, R.residual]');
end
