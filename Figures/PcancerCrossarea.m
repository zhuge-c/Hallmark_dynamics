clear
clc
result_path='./Result/';

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
msize=1000;

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
        text(x(i), y(i),txt_hallmarks{values(i)},'Color','k','HorizontalAlignment','center');
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
xlabel('Ordered by JS divergence of Hallmarks','fontsize',16);
xlim([0,11])
ylabel('Cancers','fontsize',16);
yticks(1:size(cancer,1))
yticklabels(cancer)

figure_type='-djpeg';

print([result_path,'Pancer_cross_area'],figure_type)



