function [lower_bound,upper_bound]=min_max_interval(data,method,varargin)

if method==1 % "CI"
    if isempty(varargin)
        alpha= 0.05; % 显著性水平（1 - 置信水平）
    else
        alpha=varargin{1};
    end
    n = length(data); % 数据集样本数
    mean_value = mean(data); % 计算样本均值
    std_dev = std(data); % 计算样本标准差
    
    % 计算 t 分布的临界值（两侧）
    t_critical = tinv(1 - alpha / 2, n - 1);
    
    % 计算置信区间的上下限
    margin_of_error = t_critical * (std_dev / sqrt(n));
    lower_bound = mean_value - margin_of_error;
    upper_bound = mean_value + margin_of_error;
end

if method==2 % "minmax"
    lower_bound = min(nor_data);
    upper_bound = max(nor_data);
end

if method==3 %"quantile"
    if isempty(varargin)
        quant=[0.25, 0.75];
    else
        quant=[varargin{:}];
    end
    % 计算四分位数（25% 和 75% 分位数）
    quartiles = quantile(data, quant);
    lower_bound=quartiles(1);
    upper_bound=quartiles(2);
end