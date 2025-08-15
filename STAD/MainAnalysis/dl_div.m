function js_divergence=dl_div(xi1,f1,xi2,f2)
% 使用线性插值，确保插值结果为正值
common_xi = linspace(min([xi1, xi2]), max([xi1, xi2]), 1000);
f1_interp = interp1(xi1, f1, common_xi, 'linear', 'extrap');
f2_interp = interp1(xi2, f2, common_xi, 'linear', 'extrap');

% 避免对数中的零和负值，添加一个小的正数到密度估计中
epsilon = 1e-10;
f1_interp = f1_interp + epsilon;
f2_interp = f2_interp + epsilon;

% 确保插值结果为正值
f1_interp(f1_interp < 0) = epsilon;
f2_interp(f2_interp < 0) = epsilon;

% 归一化PDF使其总和为1
f1_interp = f1_interp / sum(f1_interp * (common_xi(2) - common_xi(1)));
f2_interp = f2_interp / sum(f2_interp * (common_xi(2) - common_xi(1)));

% 计算KL散度
kl_divergence = sum(f1_interp .* log(f1_interp ./ f2_interp)) * (common_xi(2) - common_xi(1));
% disp(['KL Divergence: ', num2str(kl_divergence)]);

% 计算JS散度
M = 0.5 * (f1_interp + f2_interp);
js_divergence = 0.5 * sum(f1_interp .* log(f1_interp ./ M)) * (common_xi(2) - common_xi(1)) + ...
                0.5 * sum(f2_interp .* log(f2_interp ./ M)) * (common_xi(2) - common_xi(1));
% disp(['JS Divergence: ', num2str(js_divergence)]);