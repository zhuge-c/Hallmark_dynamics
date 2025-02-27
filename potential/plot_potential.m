
%**************************************************************%
%**************************************************************%
function [xGrid, yGrid, potential]=plot_potential(x,y,gridSize)
% 定义核密度估计的带宽
bandwidth = 0.5;

% % 创建一个网格，用于计算每个点的势能
% gridSize = 50;
xGrid = linspace(min(x), max(x), gridSize);
yGrid = linspace(min(y), max(y), gridSize);
[xGridMesh, yGridMesh] = meshgrid(xGrid, yGrid);

% 计算每个网格点的势能
potential = zeros(gridSize, gridSize);
for i = 1:gridSize
    for j = 1:gridSize
        % 计算每个网格点的核密度估计值
        distances = sqrt((x - xGridMesh(i, j)).^2 + (y - yGridMesh(i, j)).^2);
        % 核密度估计值
        potential(i, j) = sum(exp(-(distances.^2) / (2 * bandwidth^2))) / (length(x) * bandwidth * sqrt(2 * pi));
    end
end


end
