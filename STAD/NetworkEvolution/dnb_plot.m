function dnb_plot(total_dnb, local_dnb, file_path, figure_type)


figure('WindowState', 'maximized','Color', 'w'); % 创建一个白色背景的图形窗口
hold on; % 允许在同一坐标轴上绘制多个图形

% --- 3. 添加背景高亮区域 ---
% 使用 patch 函数创建一个半透明的矩形区域
% 矩形的四个顶点坐标: (x1,y1), (x2,y1), (x2,y2), (x1,y2)
patch_x = [0 100 100 0];
patch_y = [0.00125 0.00125 0.00205 0.00205];
patch(patch_x, patch_y, [1, 0.9, 0.8], ... % 颜色 (淡橙色)
      'EdgeColor', 'none', ...           % 无边框
      'FaceAlpha', 0.5);                 % 设置透明度

% --- 4. 绘制所有曲线 ---
% 首先绘制 10 条细线 (H1-H10)
num_h_curves = size(local_dnb,2);
x=1:size(local_dnb,1);
h_plots = gobjects(num_h_curves, 1); % 初始化图形句柄数组
colors = lines(num_h_curves); % 获取一组区分度高的颜色
for i = 1:num_h_curves
    h_plots(i) = plot(x', local_dnb(:, i), 'LineWidth', 1, 'Color', colors(i,:));
end

% 最后绘制加粗的 "Total" 曲线，使其显示在最上层
total_plot = plot(x, total_dnb, 'k', 'LineWidth', 2.5); % 'k' 代表黑色

% --- 5. 添加垂直虚线和文本注释 ---
% 在 x=37 处添加一条黑色垂直虚线
xline(37, '--r', 'LineWidth', 1.5);

% 在指定位置添加文本
text(10, 0.0016, 'Critical signal', 'FontSize', 16);


% --- 6. 设置坐标轴和标题 ---
% 设置坐标轴标签和图形标题
xlabel('Time Points', 'FontSize', 16);
ylabel('DIND Score', 'FontSize', 16);
title('Dynamic Evolution of DIND Scores Over Time', 'FontSize', 16, 'FontWeight', 'bold');

% 设置坐标轴的范围和刻度
xlim([0 100]);
ylim([0 0.0021]);
xticks([0 20 37 40 60 80 100]); % 自定义 X 轴刻度点
grid on; % 添加网格线，使其更清晰

ax = gca; % 获取当前坐标轴对象 (Get Current Axes)
ax.FontSize = 16; % <-- 新增：设置刻度数字的字体大小
% --- 7. 创建并美化图例 ---
% 创建图例标签
legend_labels = {'Total'};
for i = 1:num_h_curves
    legend_labels{end+1} = ['H' num2str(i)];
end

% 创建图例对象
% 将所有曲线的句柄一起传入 legend 函数
lgd = legend([total_plot; h_plots], legend_labels);

% 美化图例
lgd.NumColumns = 11;              % 设置图例列数，使其水平排列
lgd.Location = 'southoutside';   % 将图例放置在图的下方外部
lgd.Title.String = 'Hallmark';   % 为图例添加标题
lgd.Box = 'on';                 % 去掉图例的边框
lgd.FontSize = 16; % <-- 新增：设置图例各项文字的字体大小
lgd.Title.FontSize = 16; % <-- 新增：设置图例标题的字体大小
lgd.ItemTokenSize = [20, 18];
% --- 绘图结束 ---
hold off;

print([file_path,'dnb_plot'],figure_type)
