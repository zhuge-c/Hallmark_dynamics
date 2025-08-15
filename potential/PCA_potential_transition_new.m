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

fontsize = 10;

%***************************************************************%

% All_sample_X are all the expression data of the sample
% d1: gege, d2: time; d3: sample

Hallmarks={'Evading Apoptosis'	'Evading Immune Destruction'	'Genome Instability and Mutation'	'Insensitivity to Anti-Growth Signals'	'Limitless Replicative Potential'	'Reprogramming Energy Metabolism'	'Self-Sufficiency in Growth Signals'	'Sustained Angiogenesis'	'Tissue Invasion and Metastasis'	'Tumor-Promoting Inflammation'};

temp_Hallmarks ={};
for i=1:length(Hallmarks)
    temp_Hallmarks ={temp_Hallmarks{:},['H',num2str(i)]};
end

Hallmarks = temp_Hallmarks;

% figure('WindowState', 'maximized')

for cc=13:13%size(cancer_mat_file,1)
    
    disp(cancer_mat_file{cc})

    matfile=cancer_mat_file(cc);

    load(matfile{:},"All_sample_X","nor_t","cancer_t","time_point_id","DNB_score_dind")
    
    All_sample = All_sample_X(:,time_point_id,:);

    time_sample = randi([1, 100], 1, 10000);
    sample_sample = randi([1, 10000], 1, 10000);

    pca_data = zeros(size(All_sample,1), 10000);
    for s = 1:10000
        pca_data(:,s) = reshape(All_sample(:,time_sample(s),sample_sample(s)),size(All_sample,1),[]);
    end
    [stand_data,C,S]=normalize(pca_data');
    [coeff, score] = pca(stand_data);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure 
    hold on
    p1=scatter(score(:,1), score(:,2), ...
        mksize,'MarkerEdgeColor','#3BA997','Marker','.','AlphaData',i/num_cancer);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %*******************************%
    % compute potential
    [xGrid, yGrid, potential]=plot_potential(score(:,1),...
        score(:,2),50);
    figure
    surface(xGrid,yGrid,exp(-potential))
    colormap("turbo")
    view(3)
    grid on
    xlabel("PC 1", "FontSize", fontsize, "FontName",'Times New Roman')
    ylabel("PC 2", "FontSize", fontsize, "FontName",'Times New Roman')
    zlabel("Potential", "FontSize", fontsize, "FontName",'Times New Roman')
    xlim([min(xGrid),max(xGrid)])
    ylim([min(yGrid),max(yGrid)])
    zlim([0.8,1])

    ax=gca;
    ax.FontSize = fontsize;
    fig = gcf;
    fig.PaperUnits = 'centimeters';
    width_cm = 10;
    height_cm = 8;
    fig.PaperSize = [width_cm, height_cm]; 
    fig.PaperPosition = [0, 0, width_cm, height_cm]; 

    print([file_path,'potential_3d'],figure_type,'-r1200')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %*******************************%
    % potential transition 

    mean_nor = [-2.73176, 0.0800895];
    mean_cancer = [2.99423, -0.0198261];
    
    other_node =[];
    for i=1:4
        temp = mean_nor + (mean_cancer - mean_nor)*i/4;
        other_node = [other_node; temp];
    end
    
    total_node = [mean_nor;other_node];

    figure
    hold on
    colormap("turbo")
    contourf(xGrid, yGrid, exp(-potential),10,"LineColor","none"); % 20表示等高线的数量
    [U,V] = gradient(potential,0.5,0.5);
    quiver(xGrid,yGrid,U,V,2)

    plot(total_node(:,1), total_node(:,2), '-w.',"MarkerSize",20)

    for n = 1:size(total_node,1)
        text(total_node(n,1), total_node(n,2)+0.3, ['p_',num2str(n)],...
            "FontSize",fontsize,"FontName","times new roman",'Color','w')
    end
    colorbar;
    xlabel("PC 1", "FontSize", fontsize, "FontName",'Times New Roman')
    ylabel("PC 2", "FontSize", fontsize, "FontName",'Times New Roman')
    xlim([min(xGrid),max(xGrid)])
    ylim([min(yGrid),max(yGrid)])
    ax=gca;
    ax.FontSize = fontsize;
    fig = gcf;
    fig.PaperUnits = 'centimeters';
    width_cm = 10;
    height_cm = 8;
    fig.PaperSize = [width_cm, height_cm]; 
    fig.PaperPosition = [0, 0, width_cm, height_cm]; 

    print([file_path,'potential_fv'],figure_type, '-r1200')

    %*********************************************************%
    %*********************************************************%
    
    h25 = [];
    h50 = [];
    h75 = [];
    
    flag = 0;
    for i=1: size(total_node,1)
        temp_node = total_node(i,:);
        if flag == 1
            neighbor_num =2000;
            temp_distance = sum((score(:,1:2)- temp_node).^2,2);
            [val,ind] = sort(temp_distance);
            neighbor = score(ind(1:neighbor_num),1:2);
            neighbor = [neighbor,zeros(size(neighbor,1),size(score,2)-2)];
        else
            temp_ind = abs(score(:,1)-temp_node(1));
            [val, ind] = sort(temp_ind);
            nei_ind = ind(1:length(ind)*0.1);

            temp_ind2 = abs(score(:,2)-temp_node(2));
            [val2, ind2] = sort(temp_ind2);
            nei_ind2 = ind2(1:length(ind2)*0.1);
            
            total_ind = union(nei_ind, nei_ind2);
            neighbor = score(total_ind,1:2);
            neighbor = [neighbor,zeros(size(neighbor,1),size(score,2)-2)];

        end
    
        temp_stand_data=neighbor*inv(coeff);
        temp_Hallmark_exp=temp_stand_data.*S+C;
    
        p_quartiles = [0.25, 0.5, 0.75];
        Q_quartiles = quantile(temp_Hallmark_exp, p_quartiles);
    
        temp_h25 = Q_quartiles(1,:);
        temp_h50 = Q_quartiles(2,:);
        temp_h75 = Q_quartiles(3,:);
        
        h25 = [h25; temp_h25];
        h50 = [h50; temp_h50];
        h75 = [h75; temp_h75];
    
    end
    
    figure
    
    for i=1:size(h25,2)
        subplot(2,5,i)
        hold on
        x=1:size(h25,1);

%         xconf = [x x(end:-1:1)] ;         
%         yconf = [h25(:,i)', h75(end:-1:1,i)']; 
%         p = fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');
%         
%         plot(x, h50(:,i), '-k.','MarkerSize', 20)

        e = errorbar(x, h50(:,i), h50(:,i)-h25(:,i), h75(:,i)-h50(:,i),...
            '-s','MarkerSize', 3, 'MarkerEdgeColor','red','MarkerFaceColor','red');
        e.CapSize = 6; % 设置帽子大小
        e.LineWidth = 1; % 设置误差棒线宽
    
        ylabel(Hallmarks{i},'FontSize',fontsize, "FontName",'Times New Roman')
        xlabel("Time",'FontSize',fontsize, "FontName",'Times New Roman')
        xtickangle(0)
        ax=gca;
        ax.FontSize = fontsize;
        
        tt = 1:size(total_node,1);
        t_label ={};
        for t = 1:size(total_node,1)
            t_label = {t_label{:}, ['p_', num2str(t)]};
        end
        
        ylim([0.5, 1.6])
        set(gca, 'XTick',tt, 'XTickLabel',t_label)
    end

    fig = gcf;
    fig.PaperUnits = 'centimeters';
    width_cm = 20;
    height_cm = 8;
    fig.PaperSize = [width_cm, height_cm]; 
    fig.PaperPosition = [0, 0, width_cm, height_cm]; 

    print([file_path,'H_exp'],figure_type,'-r1200')

    
    


end

%*********************************************************%
% path_data=xlsread('minpath.xlsx');
% 
% figure
% [M,c]=contour(xGrid, yGrid, potential,60,'LineWidth',1);
% hold on
% plot(path_data(:,1),path_data(:,2),'linewidth',3,'Color','y')
% print([file_path,'potential_path'],figure_type)
% % score=stand_data*coeff
% pca_new=zeros(size(path_data,1),size(score,2));
% pca_new(:,1:2)=path_data(:,1:2);
% H_stand_data=pca_new*inv(coeff);
% Hallmark_exp=H_stand_data.*S+C;



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
