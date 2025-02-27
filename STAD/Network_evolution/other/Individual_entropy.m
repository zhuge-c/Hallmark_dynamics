close all
clear
clc

file_dir='./rand1/sample/';
% 获取目录下的所有文件信息
files = dir(fullfile(file_dir, '*.txt'));

% 从文件信息中提取文件名
fileNames = {files.name};

sample_num=length(fileNames);

All_sample=[];
%% 读取样本数据
for i=1:sample_num
   % 读取不同的样本数据
   temp_data=importdata([file_dir,'sample',num2str(i),'.txt']);
   All_sample(:,:,i)=temp_data;
end

% All_sample=All_sample(:,22:end,:);
% All_sample d1:基因, d2:时间, d3:样本
psize=size(All_sample);

%% 计算个体熵
A_individual_entropy=[];
for t=2:psize(2)
    t_individual_entropy=[];
    for s=1:psize(3)
        x=All_sample(:,t,s);
        r=All_sample(:,t-1,s);
        y=x-r;

        % 创建直方图并获取分箱的频率
        [counts, bin_centers] = hist(y);
        
        counts(counts==0)=[];
        % 计算频率
        total_count = sum(counts); % 总数据点数
        bin_probabilities = counts / total_count;
        
        % 计算个体熵
        individual_entropy = -sum(bin_probabilities .* log2(bin_probabilities));
        t_individual_entropy=[t_individual_entropy,individual_entropy];
    end
    A_individual_entropy=[A_individual_entropy;t_individual_entropy];
end

boxplot(A_individual_entropy')









