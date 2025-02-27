function data_filtered=Del_out_Box_Plot(data)

% 计算箱型图的上下四分位数和四分位距
Q = quantile(data, [0.25, 0.75]);
Q1 = Q(1); % 下四分位数
Q3 = Q(2); % 上四分位数
IQR = Q3 - Q1; % 四分位距

% 计算离群值的边界
lower_bound = Q1 - 1.5 * IQR; % 下界
upper_bound = Q3 + 1.5 * IQR; % 上界

% 删除离群值
data_filtered = data(data >= lower_bound & data <= upper_bound);

end
