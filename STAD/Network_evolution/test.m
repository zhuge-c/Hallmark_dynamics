
rng('default'); % 确保结果可重复

data1=time1_data;
data2=time2_data;

[coeff1, score1] = pca(data1);
[coeff2, score2] = pca(data2);
% score1=sort(score1,2,'descend');
% score2=sort(score2,2,'descend');

figure;
hold on;

% 绘制第一组数据的 PCA 得分
scatter(score1(:,1), score1(:,2), 'filled', 'b');

% 绘制第二组数据的 PCA 得分
scatter(score2(:,1), score2(:,2), 'filled', 'r');

xlabel('First Principal Component');
ylabel('Second Principal Component');
title('PCA Score Plot with Confidence Ellipses');



% 绘制第一组数据的置信椭圆
drawEllipse(score1(:,1), score1(:,2), 0.95, 'blue');

% 绘制第二组数据的置信椭圆
drawEllipse(score2(:,1), score2(:,2), 0.95, 'red');

hold off;


% 添加绘制置信椭圆的功能
function drawEllipse(x, y, confidence, color)
    data = [x y];
    mean_data = mean(data, 1);
    [V, D] = eig(cov(data));
    [D, order] = sort(diag(D), 'descend');
    V = V(:, order);
    t = linspace(0, 2 * pi, 100);
    ell = [cos(t); sin(t)]' * sqrt(chi2inv(confidence, 2)) * diag(sqrt(D));
    ell = ell * V';
    ell = bsxfun(@plus, ell, mean_data);
    plot(ell(:,1), ell(:,2), 'LineWidth', 2, 'Color', color);
end
