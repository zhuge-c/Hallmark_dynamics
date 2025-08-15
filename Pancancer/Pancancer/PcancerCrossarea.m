clear
clc
result_path='./Result/';

fontsize = 10;

data=readtable([result_path,'cancer_overleap.csv']);
cancer=table2cell(data(:,1));

overleap=table2array(data(:,2:end));
order_cancer=[];
for i=1:size(overleap,1)
    [val,ind]=sort(overleap(i,:));
    order_cancer=[order_cancer;ind];
end

corss_area=order_cancer;


colors = [
    0.0, 0.45, 0.74; % 深蓝色
    0.85, 0.33, 0.10; % 深橙色
    0.93, 0.69, 0.13; % 深黄色
    0.49, 0.18, 0.56; % 深紫色
    0.47, 0.67, 0.19; % 深绿色
    0.30, 0.75, 0.93; % 深青色
    0.64, 0.08, 0.18; % 深红色
    0.58, 0.44, 0.86; % 深蓝紫色
    0.75, 0.75, 0.00; % 深橄榄色
    0.39, 0.60, 0.64; % 深青灰色
];

marks={'o','h','s','d','p','o','h','s','d','p'};
txt_hallmarks={'H1','H2','H3','H4','H5','H6','H7','H8','H9','H10'};
msize=500;

x=1:size(corss_area,2);
x0=x;
x0(1)=x(1)-0.5;
x0(end)=x(end)+0.5;


% figure('WindowState', 'maximized')
figure('Position',[10,10,800,800])
hold on
for cc=1:size(corss_area,1)
    y=ones(1,size(corss_area,2))*cc;
    values = corss_area(cc,:);
    % 创建散点图
    for i=1:length(x)
        scatter(x(i), y(i), msize, colors(values(i),:), 'filled','Marker','s');
        text(x(i), y(i),txt_hallmarks{values(i)},'Color','k',...
            'HorizontalAlignment','center','fontsize',fontsize, 'fontname', 'times new roman');
    end
    plot(x0,y+0.5,'LineStyle','-','color','k')
end

% % 添加手动图例条目
% unique_values = unique(values);
% h = zeros(length(unique_values), 1);
% legend_entries = cell(length(unique_values), 1);
% 
% for i = 1:length(unique_values)
%     % 绘制隐藏的点，只用于图例
%     h(i) = scatter(NaN, NaN, msize, colors(unique_values(i),:), 'filled','Marker',marks{unique_values(i)});
%     legend_entries{i} = ['Hallmark ' num2str(unique_values(i))];
% end
% 
% % 创建图例
% legend(h, legend_entries, 'Location', 'northeastoutside','fontsize',16);

hold off

% 设置图形标题和轴标签

% xlabel('Divergence order of Hallmarks','fontsize',16);
xlabel('Ordered by JS divergence of Hallmarks','fontsize',fontsize, 'FontName','times new roman');
xlim([0,11])
ylabel('Cancers','fontsize',fontsize, 'FontName','times new roman');
yticks(1:size(cancer,1))
yticklabels(cancer)
ax=gca;
ax.FontSize = fontsize;

figure_type='-djpeg';

fig = gcf;

fig.PaperUnits = 'centimeters';
width_cm = 12;
height_cm = 15;
fig.PaperSize = [width_cm, height_cm]; 
fig.PaperPosition = [0, 0, width_cm, height_cm]; 

print([result_path,'Pancer_cross_area'],figure_type, '-r1200')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% [p, tbl, stats]= friedman_test(corss_area);
% 
% outfile = './Result/area_friedman_result.txt';
% 
% test_output(p, tbl, stats, outfile)
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% num_cancers = size(corss_area,1);
% 
% average_ranks = stats.meanranks;
% 
% % 4. 绘制条形图
% % 为了更清晰地展示，我们先按平均排名对Hallmark进行排序
% [sorted_ranks, sort_order] = sort(average_ranks, 'ascend');
% sorted_names = txt_hallmarks(sort_order);
% 
% % 创建图形
% figure;
% bar(sorted_ranks);
% 
% % 添加标题和标签
% title('Average Rank of Hallmarks Across 15 Cancer Types');
% ylabel('Average Rank (Lower is Better)');
% xlabel('Hallmarks of Cancer');
% set(gca, 'xtick', 1:length(txt_hallmarks), 'xticklabel', sorted_names);
% xtickangle(45); % 让x轴的标签倾斜，避免重叠
% grid on; % 添加网格线
% 
% % 在图上标注弗里德曼检验的p值
% text_str = sprintf('Friedman Test, p = %.2e', p);
% % 获取图形的坐标轴范围，以确定文本放置的位置
% ax = gca;
% text(ax.XLim(2)*0.6, ax.YLim(2)*0.9, text_str, 'FontSize', 12, 'FontWeight', 'bold');
% 
% % 调整图形边距
% set(gca, 'LooseInset', get(gca, 'TightInset'));
% 
% 
% 
% 
