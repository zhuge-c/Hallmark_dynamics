% clear
clc
addpath('./MainAnalysis/');
figure_type='-djpeg';
figure_size=[10,10,800,400];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_Hallmarks ={};
for i=1:length(Hallmarks)
    temp_Hallmarks ={temp_Hallmarks{:},['H',num2str(i)]};
end

Hallmarks = temp_Hallmarks;

file_path = './Result/result/';

result_file = "./Result/result/result.xlsx";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rand_All_sample=zeros(group_num,gene,t_num,cell_num);
for i=1:group_num
    rs=group_id(i,:);
%     rs = randperm(sample_num,cell_num);
    Rand_All_sample(i,:,:,:)=All_sample(:,:,rs);
end

reference_num=6;
edge_weights_threshold=0.25;
addpath('./DIND/');
[DNB_score_dind,local_dnb]=Compute_DNB(Rand_All_sample,reference_num,edge_weights_threshold,figure_type,file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1
Plot_Hallmarks_Expression(Hallmarks,All_sample_X,DNB_score_dind,...
    Parameters.dt,Parameters.total_t,figure_type,file_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure 2
[cross_areas, cancer_peak, nor_peak,JS_div]=Plot_Norm_Cancer_Difference_fill_DF(Hallmarks, ...
    All_sample_X,nor_t,cancer_t,figure_type,file_path);

data_frame0 = table(Hallmarks', cross_areas);
writetable(data_frame0, result_file, "Sheet","cross_area")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure JS div
fontSize=16;
% figure('WindowState', 'maximized')
figure('Position',[10,10,1200,800])
b=bar(JS_div,0.5);
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;
% 格式化标签，只显示小数点后两位
labels1 = arrayfun(@(x) sprintf('%.3f', x), b.YData, 'UniformOutput', false);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','FontSize',fontSize)
ylim([0,0.8])
xlabel("Hallmarks",'fontsize',fontSize,'fontname','times new roman')
xticks(1:length(Hallmarks))
xticklabels(Hallmarks)
ylabel("JS divergence",'fontsize',fontSize,'fontname','times new roman')
ax=gca;
ax.FontSize=fontSize;
print([file_path,'JSDiv'],figure_type)

data_frame1 = table(Hallmarks', JS_div,'VariableNames', {'Hallmark_Name', 'JS_Divergence'});

writetable(data_frame1, result_file, "Sheet","JS_div")
writetable(data_frame1, "./Result/jsd.csv")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 3
[maxval,maxindex]=max(JS_div);
[minval,minindex]=min(JS_div);
analysis_gene=[minindex,maxindex];

% Plot_Exp_DNB4(Hallmarks,All_sample_X,DNB_score_dind,analysis_gene, ...
%         Parameters.dt,Parameters.total_t,figure_type,file_path)

Plot_Interval_Exp_DNB(Hallmarks,All_sample_X(:,:,group_id(1,1:100)),...
        DNB_score_dind,analysis_gene, nor_t,cancer_t, time_point, ...
        Parameters.dt,Parameters.total_t,figure_type,file_path)

% time_val_global = Plot_Interval_Exp_DNB2(Hallmarks,All_sample_X(:,:,group_id(1,1:100)),...
%         DNB_score_dind,nor_t,cancer_t, time_point, ...
%         Parameters.dt,Parameters.total_t,figure_type,file_path);
% 
% data_frame6 = table(Hallmarks', time_val_global);
% writetable(data_frame6, result_file, "Sheet","Global_DNB_t_Cancer_t")
% 
% Plot_Interval_Exp_Local_DNB(Hallmarks,All_sample_X(:,:,group_id(1,1:100)),...
%         local_dnb,analysis_gene, nor_t,cancer_t, time_point, ...
%         Parameters.dt,Parameters.total_t,figure_type,file_path)

% time_val_local = Plot_Interval_Exp_Local_DNB2(Hallmarks,All_sample_X(:,:,group_id(1,1:100)),...
%         local_dnb, nor_t,cancer_t, time_point, ...
%         Parameters.dt,Parameters.total_t,figure_type,file_path);
% 
% data_frame7 = table(Hallmarks', time_val_local);
% writetable(data_frame7, result_file, "Sheet","Local_DNB_t_Cancer_t")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure 5

[hist_data, Box_Data]= Plot_Norm_Cancer_Difference_bar_and_box(All_sample_X,figure_type,Hallmarks,file_path);

box_title = {};
for i=1:length(Hallmarks)
box_title = {box_title{:}, ['Normal H',num2str(i)], ['Cancer H',num2str(i)]};
end

data_frame2 = table(box_title', Box_Data);
writetable(data_frame2, result_file, "Sheet","BoxData_N_C_expression")
writetable(data_frame2, "./Result/hallmark_data.csv")

data_frame2_2 = table(Hallmarks', hist_data);

writetable(data_frame2_2, result_file, "Sheet","histData_N_C_expression")

save("BoxData.mat","Box_Data")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表达检验

% 创建一个表格来存储结果
results = table('Size', [10 3], 'VariableTypes', ...
    {'string', 'double', 'string'}, 'VariableNames',...
    {'FeaturePair', 'PValue', 'Significant'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 表达检验
% 设置显著性水平 alpha
alpha = 0.05;

for i = 1:10
    % 确定 'normal' 和 'cancer' 数据所在的列索引
    normal_col_index = 2*i - 1;
    cancer_col_index = 2*i;
    % 提取两组数据
    normal_data = Box_Data(normal_col_index,:)';
    cancer_data = Box_Data(cancer_col_index,:)';
    % 执行双样本t检验 (ttest2): h = 1 表示在 alpha 水平下拒绝原假设 (即均值有显著差异) % h = 0 表示不能拒绝原假设 (即均值无显著差异)
    % p 是 p-value
%     [h, p] = ttest2(normal_data, cancer_data);
%     执行Mean U 检验
    [p,h,stats] = ranksum(normal_data,cancer_data);
    % 获取当前特征对的名称
    feature_pair_name = sprintf('H%d_normal vs H%d_cancer', i, i);
    % 判断是否显著
    if h == 1
        is_significant = '是';
    else
        is_significant = '否';
    end
    
    % 将结果存入表格
    results.FeaturePair(i) = feature_pair_name;
    results.PValue(i) = p;
    results.Significant(i) = is_significant;
    
    % 在命令窗口打印单次结果
    fprintf('特征对: %s\n', feature_pair_name);
    fprintf('  p-value: %.4f\n', p);
    fprintf('  在 alpha=%.2f 水平下是否显著: %s\n', alpha, is_significant);
    fprintf('--------------------------------------------------\n');
end
disp('所有特征对的显著性检验结果汇总:');
disp(results);

writetable(results, result_file, "sheet","BoxData_nor_cancer_test")



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 4
Rand_All_sample_X=zeros(group_num,gene,time_seq,cell_num);

for i=1:group_num
    rs=group_id(i,:);
    Rand_All_sample_X(i,:,:,:)=All_sample_X(:,:,rs);
end
All_Hallmarks_t1t2=Plot_DNB_Expression_Time_Difference(DNB_score_dind,Rand_All_sample_X,nor_t,cancer_t,...
    time_point,Hallmarks,figure_type,file_path);

t1t2_dif =[];
for i=1:size(All_Hallmarks_t1t2,1)
    t1t2_dif = [t1t2_dif; All_Hallmarks_t1t2{i}];
end

data_frame3 = table(Hallmarks', t1t2_dif);
writetable(data_frame3, result_file, "Sheet","BoxData_t1t2_dif")
writetable(data_frame3, "./Result/t1-t2.csv")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dict_boxplot(All_Hallmarks_t1t2)

xlabel("Hallmarks")
xticks(1:length(Hallmarks))
xticklabels(Hallmarks)
ylim([1,6])

ylabel("Difference of times t1-t2")

% 时间差检验
gene=size(All_Hallmarks_t1t2,1);
p_result = [];
for i=1:gene
    data = All_Hallmarks_t1t2{i};
    [p,~,stats] = signrank(data, 0, 'tail', 'right');
    p_result = [p_result; stats.signedrank, p];

    if p < 0.001
        sig_text = '***';
    elseif p < 0.01
        sig_text = '**';
    elseif p < 0.05
        sig_text = '*';
    else
        sig_text = 'ns'; % 'ns' 代表 'not significant' (不显著)
    end

    y_position = max(data) + 0.2;

    text(i, y_position, sig_text, ...
         'HorizontalAlignment', 'center', ... % 水平居中
         'FontSize', 14, ...                   % 字体大小
         'FontWeight', 'bold', ...             % 粗体
         'Color', 'red');                      % 颜色
end

data_frame4 = table({"sign", "p-value"}',p_result');
writetable(data_frame4, result_file, "Sheet","BoxData_t1t2_test")

print([file_path,'Hallmarks_t2t1_2'],figure_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNB line
figure
hold on 
plot(1:size(DNB_score_dind,2),DNB_score_dind)

for i = 1:size(local_dnb,2)
    plot(local_dnb(:,i))
end

xlabel("Time")
ylabel("DIND")

xt=0:ceil(size(DNB_score_dind,2)/10):size(DNB_score_dind,2);
xt(1)=1;
xlab=1:size(DNB_score_dind,2);
xlab=xlab(xt);

xticks(xt)
xticklabels(xt)
title('DIND trends over time')

dnb_data =[DNB_score_dind;local_dnb'];
label = {"Total", Hallmarks{:}};
dnb_data_frame = table(label',dnb_data);

writetable(dnb_data_frame, result_file, "sheet","dnb_data")
writetable(dnb_data_frame, "./Result/DNB.csv")
print([file_path,'DIND_2'],figure_type)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dnb_plot(DNB_score_dind, local_dnb, file_path, figure_type)








