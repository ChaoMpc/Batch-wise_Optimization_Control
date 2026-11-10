% =========================================================================
% 从已保存的 .mat 文件中提取“相邻批次最大改变量”
% =========================================================================
clear; clc;
fprintf('=== 正在从 oh_sensitivity_workspace.mat 提取相邻批次指标 ===\n');

S = load('oh_sensitivity_workspace.mat');
Q_list = S.Q_list;
R_list = S.R_list;
Data_oh_all = S.Data_oh_all;

nR = length(R_list);
nQ = length(Q_list);

% 预分配
max_adj_dY = nan(nR, nQ);
max_adj_dT = nan(nR, nQ);

for i = 1:nR
    for j = 1:nQ
        Data = Data_oh_all{i, j};
        if isempty(Data), continue; end

        B = size(Data.Y_hist, 3); % 批次数
        dY_norms = zeros(B-1, 1);
        dT_norms = zeros(B-1, 1);

        % 计算相邻批次之间的 2-范数改变量
        for b = 2:B
            dY_vec = reshape(Data.Y_hist(:,:,b) - Data.Y_hist(:,:,b-1), [], 1);
            dY_norms(b-1) = norm(dY_vec, 2);

            dT_vec = reshape(Data.T_hist(:,:,b) - Data.T_hist(:,:,b-1), [], 1);
            dT_norms(b-1) = norm(dT_vec, 2);
        end

        % 取所有相邻批次变化中的最大值
        max_adj_dY(i, j) = max(dY_norms);
        max_adj_dT(i, j) = max(dT_norms);
    end
end

% 打印控制台直观表格
fprintf('\n%-10s | %-8s | %-12s | %-6s | %-18s | %-18s\n', ...
    'R_upd', 'Q_trk', 'Error', 'Hits', 'MaxAdjBatch_dY', 'MaxAdjBatch_dT');
fprintf(repmat('-', 1, 85)); fprintf('\n');
for i = 1:nR
    for j = 1:nQ
        fprintf('%-10.0e | %-8d | %-12.4e | %-6d | %-18.4e | %-18.4e\n', ...
            R_list(i), Q_list(j), S.res_err(i,j), S.res_hits(i,j), ...
            max_adj_dY(i,j), max_adj_dT(i,j));
    end
end
