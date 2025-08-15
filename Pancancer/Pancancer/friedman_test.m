function [p, tbl, stats]= friedman_test(hallmark_order_matrix)

num_cancers = size(hallmark_order_matrix,1);
num_hallmarks = size(hallmark_order_matrix,2);

rank_matrix = zeros(num_cancers, num_hallmarks);
for i = 1:num_cancers 
    for k = 1:num_hallmarks 

        hallmark_code = hallmark_order_matrix(i, k);

        rank_matrix(i, hallmark_code) = k;
    end
end

disp('转换后的秩矩阵 (用于Friedman检验):');
disp(rank_matrix(1:5, :));


% 3. 执行弗里德曼检验
% [p, tbl, stats] = friedman(Y, reps, 'off')
% Y: m x n 矩阵, m 是区组数 (癌症数), n 是处理组数 (Hallmark数)
% reps: 每个单元格的重复次数，这里是1
% 'off': 不显示图形输出
% p: 检验的p值
% tbl: 一个包含方差分析表的单元格数组
% stats: 一个包含额外统计信息的结构体
[p, tbl, stats] = friedman(rank_matrix, 1, 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4. 显示并解释结果
fprintf('\n--- 弗里德曼检验结果 ---\n');
fprintf('P-value: %e\n', p);

% 解释 p 值
if p < 0.05
    fprintf('结果解释: P值小于0.05，我们拒绝原假设。\n');
    fprintf('这表明在15种癌症中，10个Hallmark的排序存在统计上显著的一致性。\n');
    fprintf('也就是说，这种排序模式不太可能是由随机巧合产生的。\n');
else
    fprintf('结果解释: P值大于或等于0.05，我们不能拒绝原假设。\n');
    fprintf('这表明没有足够的证据证明在15种癌症中，Hallmark的排序存在一致的模式。\n');
end

% 显示详细的统计表
fprintf('\n详细统计表 (ANOVA-like Table):\n');
disp(tbl);

% 显示每个 Hallmark 的秩和 (Rank Sums)
% 这对于后续分析（post-hoc）很重要，可以帮助看出哪些 Hallmark 倾向于排名更高（秩和更小）或更低（秩和更大）
fprintf('\n每个 Hallmark 的秩和 (列1-10 对应 Hallmark 1-10):\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




