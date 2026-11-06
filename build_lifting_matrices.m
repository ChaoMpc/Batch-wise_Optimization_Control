function Lift = build_lifting_matrices(Model, Params)
% =========================================================================
% build_lifting_matrices
%
% Fixed nominal control-segment prediction.
%
% t = [t_1; ...; t_H].
%
% At control step k:
%   x_pre(k)  : state BEFORE applying t_k
%   x_post(k) : state AFTER applying t_k
%
% Physical predictions:
%   v_abs = Phi_v*t + v_offset
%   y_abs = Phi_y*t + y_offset
%
% y_abs contains END-of-interval outputs.
% =========================================================================

    H = Params.Ble-Params.ctrl_start;

    nx = Model.nx;
    nu = Model.A_lv;

    ni = Params.n_in_opt;
    no = Params.n_out;

    nt = nu*H;

    A = Model.barA;
    B = Model.barB;

    Mv = diag(Model.X_std)*Model.P;
    My = diag(Model.Y_std)*Model.Q*Model.barC;

    Phi_pre = zeros(nx*H,nt);
    offset_pre = zeros(nx*H,1);

    Phi_post = zeros(nx*H,nt);
    offset_post = zeros(nx*H,1);

    Phi_v = zeros(ni*H,nt);
    v_offset = zeros(ni*H,1);

    Phi_y = zeros(no*H,nt);
    y_offset = zeros(no*H,1);

    state_map = zeros(nx,nt);
    state_offset = Model.x0_ctrl(:);

    for k = 1:H
        cx = (k-1)*nx+(1:nx);
        cv = (k-1)*ni+(1:ni);
        cy = (k-1)*no+(1:no);
        ct = (k-1)*nu+(1:nu);

        j = Params.ctrl_start+k;

        % Current state, before applying the controlled input.
        Phi_pre(cx,:) = state_map;
        offset_pre(cx) = state_offset;

        % Current physical input.
        Phi_v(cv,ct) = Mv;

        v_offset(cv) = ...
            Params.u_base(Model.input_idx,j)+Model.X_mean;

        % Propagate nominal state.
        state_map = A*state_map;
        state_map(:,ct) = state_map(:,ct)+B;

        state_offset = A*state_offset;

        Phi_post(cx,:) = state_map;
        offset_post(cx) = state_offset;

        % Output at physical time j.
        Phi_y(cy,:) = My*state_map;

        y_offset(cy) = ...
            Params.y_base(:,j+1)+Model.Y_mean+My*state_offset;
    end

    arrays = {Phi_pre,offset_pre,Phi_post,offset_post, ...
              Phi_v,v_offset,Phi_y,y_offset};

    for j = 1:numel(arrays)
        if any(~isfinite(arrays{j}(:)))
            error('build_lifting_matrices: nonfinite prediction matrices.');
        end
    end

    Lift = struct();

    Lift.H = H;
    Lift.nt = nt;

    Lift.Phi_x_pre = Phi_pre;
    Lift.x_pre_offset = offset_pre;

    Lift.Phi_x_post = Phi_post;
    Lift.x_post_offset = offset_post;

    Lift.Phi_v = Phi_v;
    Lift.v_offset = v_offset;

    Lift.Phi_y = Phi_y;
    Lift.y_offset = y_offset;

    Lift.input_columns = Params.ctrl_start+(1:H);
    Lift.output_columns = Lift.input_columns+1;
end
