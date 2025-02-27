function A_individual_entropy=DNB_Individual_Entropy(All_sample)
% All_sample 为所有的样本数据
% 第一维： 基因个数/ID
% 第二维： 时间跨度
% 第三维： 样本个数

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


