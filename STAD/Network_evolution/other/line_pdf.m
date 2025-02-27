function [x,pdf_normal]=line_pdf(data)
% 假设是数据服从正态分布

% 
% % 绘制直方图
% histogram(data, 'Normalization', 'pdf'); % 'pdf' 参数用于归一化直方图

% 计算正态分布的 PDF
mu = mean(data); % 均值
sigma = std(data); % 标准差
x = linspace(min(data), max(data), 100); % 创建 x 值范围
pdf_normal = (1 / (sigma * sqrt(2 * pi))) * exp(-((x - mu).^2) / (2 * sigma^2));
end
