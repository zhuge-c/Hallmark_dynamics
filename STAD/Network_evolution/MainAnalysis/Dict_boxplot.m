function Dict_boxplot(time_dif)

gene=size(time_dif,1);

figure('WindowState', 'maximized')

hold on; % 保持图的当前内容以添加更多的图形

% 循环遍历字典中的每个键值对并绘制箱线图
for i=1:gene
    data = time_dif{i};
    boxplot(data, 'positions', i);
end

hold off; % 结束在图上添加内容
