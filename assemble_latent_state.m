function x = assemble_latent_state(History, Model)
% =========================================================================
% assemble_latent_state
%
% History.U(:,1)：当前观测输出得分。
% History.T(:,1)：先前已施加输入得分。
%
% x_a(k) = [u_a(k),...,u_a(k-n+1),
%           t_a(k-1),...,t_a(k-n+1)]'.
%
% 该状态用于统一 MIMO ARX 的状态空间实现；跨通道耦合由 barA/barB 处理。
% =========================================================================

    n = Model.n_arx;
    A = Model.A_lv;
    block = Model.n_state_per_lv;

    if ~isequal(size(History.U), [A, n]) || ...
       ~isequal(size(History.T), [A, n])
        error('assemble_latent_state: invalid history dimensions.');
    end

    x = zeros(Model.nx, 1);

    for a = 1:A
        idx = (a-1)*block + (1:block);

        x(idx(1:n)) = History.U(a, :).';

        if n > 1
            x(idx(n+1:end)) = History.T(a, 1:n-1).';
        end
    end
end
