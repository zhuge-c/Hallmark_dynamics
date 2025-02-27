%******************************************%
%*** 区分每一个癌症的正常与癌症状态下的PCA效果
%                                 *********%
%******************************************%
clear
clc
file_path='./result/';
figure_type='-djpeg';
rng('default'); % 确保结果可重复

result_path='./';

num_cancer=15;

load ./All_Expression.mat

load([result_path,'cancer_file'],'cancer_folder','cancer_mat_file')

dt=0.1;

mksize=5;

%***************************************************************%

% All_sample_X are all the expression data of the sample
% d1: gege, d2: time; d3: sample

Hallmarks={'Evading Apoptosis'	'Evading Immune Destruction'	'Genome Instability and Mutation'	'Insensitivity to Anti-Growth Signals'	'Limitless Replicative Potential'	'Reprogramming Energy Metabolism'	'Self-Sufficiency in Growth Signals'	'Sustained Angiogenesis'	'Tissue Invasion and Metastasis'	'Tumor-Promoting Inflammation'};

% figure('WindowState', 'maximized')
for cc=1:1%size(cancer_mat_file,1)

    matfile=cancer_mat_file(cc);

    load(matfile{:},"All_sample_X","nor_t","cancer_t","time_point_id","DNB_score_dind")
    
    time_num=size(All_sample_X,2);
    num_sample=size(All_sample_X,3);

    %* 
    % cancer and normal data
    % *%
    cancer_data=reshape(All_sample_X(:,cancer_t,:),size(All_sample_X,1),size(All_sample_X,3));
    nor_data=reshape(All_sample_X(:,nor_t,:),size(All_sample_X,1),size(All_sample_X,3));

    pca_data=[nor_data';cancer_data'];
    [stand_data,C,S]=normalize(pca_data);
    [coeff, score] = pca(stand_data);
    nor_score=score(1:size(nor_data,2),:);
    cancer_score=score(1+size(nor_data,2):end,:);

    %*************************%
    person_sample=zeros(num_sample,1000,2);
    for i=1:num_sample
        sample=reshape(All_sample_X(:,:,i),size(All_sample_X,1),size(All_sample_X,2));
        stand_sample=normalize(sample','center',C,'scale',S);
        temp_pca_sample=stand_sample*coeff;
        person_sample(i,:,:)=temp_pca_sample(:,1:2);
    end
    ppp=reshape(person_sample(1,:,:),size(person_sample,2),size(person_sample,3));
    px=ppp(:,1);
    py=ppp(:,2);

    %*******************************%
    p1=scatter(nor_score(:,1), nor_score(:,2), ...
        mksize,'MarkerEdgeColor','#3BA997','Marker','.','AlphaData',i/num_cancer);

    p2=scatter(cancer_score(:,1), cancer_score(:,2), ...
        mksize,'MarkerEdgeColor','#B291B5','Marker','.','AlphaData',i/num_cancer);

    [xGrid, yGrid, potential]=plot_potential([nor_score(:,1);cancer_score(:,1)],...
        [nor_score(:,2);cancer_score(:,2)],50);

    figure
    surface(xGrid,yGrid,exp(-potential))
    colormap("turbo")
    view(3)
    grid on
    xlabel("PC 1")
    ylabel("PC 2")
    zlabel("potential")
    xlim([-4,4])
    ylim([-2.5,2.5])
    zlim([0.8,1])
    print([file_path,'potential_3d'],figure_type)


    figure
    colormap("turbo")
    [M,c]=contourf(xGrid, yGrid, exp(-potential),10,"LineColor","none");
    xlabel("PC 1")
    ylabel("PC 2")
    xlim([-4,4])
    ylim([-2.5,2.5])
    print([file_path,'potential_face'],figure_type)

%     figure('Position',[10,10,700,600]);   

    figure
    hold on
    colormap("turbo")
    contourf(xGrid, yGrid, exp(-potential),10,"LineColor","none"); % 20表示等高线的数量
    [U,V] = gradient(potential,0.5,0.5);
    quiver(xGrid,yGrid,U,V)
    xlabel("PC 1")
    ylabel("PC 2")
    xlim([-4,4])
    ylim([-2.5,2.5])

    print([file_path,'potential_fv'],figure_type)

end

%*********************************************************%
path_data=xlsread('minpath.xlsx');

figure
[M,c]=contour(xGrid, yGrid, potential,60,'LineWidth',1);
hold on
plot(path_data(:,1),path_data(:,2),'linewidth',3,'Color','y')
print([file_path,'potential_path'],figure_type)
% score=stand_data*coeff
pca_new=zeros(size(path_data,1),size(score,2));
pca_new(:,1:2)=path_data(:,1:2);
H_stand_data=pca_new*inv(coeff);
Hallmark_exp=H_stand_data.*S+C;

figure('WindowState', 'maximized')
hold on
hallmark_mark={'o','p','s','d','*','h','<','>','^','v'};
for i=1:size(Hallmark_exp,2)
%     plot(Hallmark_exp(:,i),'linestyle','--','marker',hallmark_mark{i},'LineWidth',2,'MarkerSize',8)
    subplot(2,5,i)
    plot(Hallmark_exp(:,i))
    title(Hallmarks{i},'FontSize',11)
    ylim([0,1.5])
end
ylim([0.4,1.3])
% legend(Hallmarks,"fontsize",12,'NumColumns',1,'Location','eastoutside')

print([file_path,'H_exp'],figure_type)


%**************************************************************%
%**************************************************************%
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

%**************************************************************%
%**************************************************************%

function result=create_colors
% 定义组合的值
values =linspace(0.1,0.9,5);

% 使用 ndgrid 生成笛卡尔积
[A, B, C] = ndgrid(values, values, values);

% 展开并组合为矩阵
result = [A(:), B(:), C(:)];

end
