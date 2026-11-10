function [History, Observation] = update_score_history( ...
    History, y_next, v_applied, input_column, Model, Params)
% =========================================================================
% update_score_history
%
% 当前物理输入施加后更新得分历史。
%
% 使用实际施加的物理输入（包括执行器限幅）。
% =========================================================================

    j = input_column;
    output_column = j + 1;

    vdev = v_applied(Model.input_idx) - ...
           Params.u_base(Model.input_idx, j);

    vstd = (vdev - Model.X_mean) ./ Model.X_std;
    t_applied = Model.R_x.' * vstd;

    ydev = y_next - Params.y_base(:, output_column);
    ystd = (ydev - Model.Y_mean) ./ Model.Y_std;
    u_next = Model.R_u.' * ystd;

    r_next = ystd - Model.Q * u_next;

    History.U = [u_next, History.U(:, 1:end-1)];
    History.T = [t_applied, History.T(:, 1:end-1)];

    Observation = struct();

    Observation.t_applied = t_applied;
    Observation.u_next = u_next;
    Observation.r_next = r_next;
end
