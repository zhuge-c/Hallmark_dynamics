clear
clc
result_path='./Result/';

fontsize = 10;

data=readtable([result_path,'cancer_time_t1t2.csv']);
cancer=table2cell(data(:,1));
times_diff=table2array(data(:,2:end));
order_cancer=[];
for i=1:size(times_diff,1)
    [val,ind]=sort(times_diff(i,:),'descend');
    order_cancer=[order_cancer;ind];
end

time_diff=order_cancer;

cancer_id=data(:,2);
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

x=1:size(time_diff,2);
x0=x;
x0(1)=x(1)-0.5;
x0(end)=x(end)+0.5;


% figure('WindowState', 'maximized')
figure('Position',[10,10,800,800])
hold on
for cc=1:size(time_diff,1)
    y=ones(1,size(time_diff,2))*cc;
    values = time_diff(cc,:);
    % 创建散点图
    for i=1:length(x)
        scatter(x(i), y(i), msize, colors(values(i),:), 'filled','Marker','s');
        text(x(i), y(i),txt_hallmarks{values(i)},'Color','k',...
            'HorizontalAlignment','center','fontsize',fontsize, 'fontname', 'times new roman');
    end

    plot(x0,y+0.5,'LineStyle','-','color','k')
end

hold off

% 设置图形标题和轴标签

xlabel('Time order of Hallmrks','fontsize',fontsize);
xlim([0,11])

ylabel('Cancers','fontsize',fontsize);
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

print([result_path,'Pancer_t2t1'],figure_type, '-r1200')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [p, tbl, stats]= friedman_test(time_diff);
% 
% 
% outfile = './REsult/time_friedman_result.txt';
% 
% test_output(p, tbl, stats, outfile)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% num_cancers = size(time_diff,1);
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
