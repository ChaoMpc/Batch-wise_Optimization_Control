function SS = realize_latent_model(ARX)
% =========================================================================
% realize_latent_model
%
% 统一 MIMO ARX 的状态空间实现。
%
% 状态：
%   x_a(k) = [u_a(k),...,u_a(k-n+1),
%             t_a(k-1),...,t_a(k-n+1)]'
%
% 模型：
%   x(k+1) = barA*x(k) + barB*t(k) + barE*eta(k)
%   u(k)   = barC*x(k)
%
% 其中 barA 的 off-diagonal 块显式包含全部潜变量通道之间的耦合。
% =========================================================================

    A = ARX.A_lv;
    n = ARX.n_arx;

    block = 2*n - 1;
    nx = A * block;

    if size(ARX.Theta_matrix, 1) ~= 2*A*n || ...
       size(ARX.Theta_matrix, 2) ~= A
        error('realize_latent_model: Theta_matrix has inconsistent dimensions.');
    end

    barA = zeros(nx, nx);
    barB = zeros(nx, A);
    barC = zeros(A, nx);
    barE = zeros(nx, A);

    for a = 1:A
        ia = (a-1)*block + (1:block);

        % 本通道输出历史状态的移位寄存器。
        for ell = 2:n
            barA(ia(ell), ia(ell-1)) = 1;
        end

        % 本通道输入记忆状态的移位寄存器。
        if n > 1
            barB(ia(n+1), a) = 1;

            for ell = n+2:block
                barA(ia(ell), ia(ell-1)) = 1;
            end
        end

        % 输出与 innovation 注入。
        barC(a, ia(1)) = 1;
        barE(ia(1), a) = 1;

        theta = ARX.Theta{a};

        % 统一 MIMO ARX：输出 u_a 依赖全部 b=1,...,A 的潜变量。
        for b = 1:A
            ib = (b-1)*block + (1:block);

            iu = (b-1)*n + (1:n);
            it = A*n + (b-1)*n + (1:n);

            alpha = theta(iu);   % alpha_{ab,1..n}
            beta  = theta(it);   % beta_{ab,1..n}

            % u_b(k),...,u_b(k-n+1) 对 u_a(k+1) 的影响。
            barA(ia(1), ib(1:n)) = alpha.';

            % t_b(k-1),...,t_b(k-n+1) 对 u_a(k+1) 的影响。
            if n > 1
                barA(ia(1), ib(n+1:end)) = beta(2:end).';
            end

            % 当前输入 t_b(k) 对 u_a(k+1) 的影响。
            barB(ia(1), b) = beta(1);
        end
    end

    SS = struct();

    SS.barA = barA;
    SS.barB = barB;
    SS.barC = barC;
    SS.barE = barE;

    SS.nx = nx;
    SS.n_state_per_lv = block;

    SS.open_loop_spectral_radius = max(abs(eig(barA)));
end
