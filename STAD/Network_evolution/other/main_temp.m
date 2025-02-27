clear
clc

tic

%%%%%%%%%%%%%%%%%%%%%%%%
% 运行参数设置
rand_num=1;

t_num=100;

sample_num=500;

%%%%%%%%%%%%%%%%%%%%%%%
%% 网络数据读取

load Hallmarks.mat

% 读取的数据为：第i行为第i个Hallmark对其余Hallmark的调控关系
% 为适应程序，我们需要改成第i行为所有Hallmark对第i个基因的调控
noc=1e4;
[V_orginal,VNT]=xlsread("normal_Positive.csv");

V_orginal=V_orginal'/noc;

[V_cancer,VCT]=xlsread("cancer_Positive.csv");

V_cancer=V_cancer'/noc;

network=(V_cancer+V_orginal)/2;

%%%%%%%%%%%%%%%%%%%%%%
%% 参数设置
gene=size(V_cancer,1); % 总的基因个数

Parameters.alpha=ones(gene,gene);

Parameters.lambda=3.8*ones(1,gene); % 基因的最大表达率

Parameters.n=2; % Hill 系数

Parameters.rho1=0.1; % 最小调控促进能力

Parameters.sigma=0.05; % 噪声扰动强度

Parameters.tau=1; % 随机扰动的相关事件

Parameters.theta=sum(V_orginal)/2;

Parameters.dt=0.1; % 单位时间长度

Parameters.total_t=100; %总的时间跨度

time_seq=Parameters.total_t/Parameters.dt;

dt=time_seq/t_num;

time_id=1:time_seq;

time_point=0:Parameters.dt:Parameters.total_t;
time_point=time_point(1:end-1);

time_point_id=time_id(1:dt:end);

nor_t=floor(time_seq*1/3);
nor_can_t=floor(time_seq*2/3);

%% 初值赋予
% X0_distribution=reshape(All_sample_X(:,end,:),size(All_sample_X,1),size(All_sample_X,3));
% save("X0_distribution","X0_distribution")

load X0_distribution.mat

%% 模型求解
load ./orginal_stable/orginal_distribution.mat

Rand_All_sample=zeros(rand_num,gene,t_num,sample_num);

for ran=1:rand_num
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % All_sample_X 为所有输入样本表达数据,  d1:基因,  d2:时间, d3:样本
    % All_sample 为选取时间点数据， 时间点为time_point_id
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    All_sample=zeros(gene,t_num,sample_num);
    All_sample_X=zeros(gene,time_seq,sample_num);
    
    % stationary_distribution(:,:) = All_sample_X(:,end,:);
    
    %%%%%%%%%%%%%%% 文件夹管理
    rand_dir=['./rand',num2str(ran),'/'];
    if exist(rand_dir, 'dir') ~= 7
            mkdir(rand_dir)
    end

    %% 随机生成 sample_num 个样本
    for r=1:sample_num
        clc
        r
        % 初始 eta 与 X0 设置
        eta0=0*rand(1,gene);

        X0 = orginal_distribution(:,randsample(size(orginal_distribution,2),1,true))';
    
        % 网络变化
        V1=V_orginal;
        V2=V_cancer;
    
        %% 模型计算
        [All_X,All_eta]=Node_SDE_modelV2(V1,V2,X0,eta0,Parameters);
        
        mean_all_X=All_X(time_point_id,:);
        for i=1:length(time_point_id)
            mean_all_X(i,:)=mean(All_X(time_point_id(i):time_point_id(i)+9,:));
        end

        All_sample(:,:,r)=mean_all_X';
%         All_sample(:,:,r)=All_X(time_point_id,:)';
        All_sample_X(:,:,r)=All_X';

        %% 保存样本数据

        %%%%%%%%%%%%%%% 文件夹管理
        sample_dir=[rand_dir,'sample/'];
        if exist(sample_dir, 'dir') ~= 7
            mkdir(sample_dir)
        end

        writematrix(All_X',[sample_dir,'sample',num2str(r),'.txt'],"Delimiter",'\t');
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 保存每次随机实验结果
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rand_All_sample(ran,:,:,:)=All_sample;

    %%%%%%%%%%%%%%% 文件夹管理
    fig_dir=[rand_dir,'fig/'];
    if exist(fig_dir,'dir')~=7
        mkdir(fig_dir)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 简单的网络随时间的变化
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    for g=1:gene
        temp_d=reshape(All_sample_X(g,:,:),size(All_sample_X,2),size(All_sample_X,3));
        subplot(2,5,g)
        hold on
        for i=1:sample_num
            plot(Parameters.dt:Parameters.dt:Parameters.total_t,temp_d(:,i))
        end
        title(Hallmarks{g})
        xlabel("Time")
        ylabel("Expression of Hallmark")
        ylim([0,max(max(temp_d))*1.1])
    end
    sgtitle("初值随机扰动, 考虑随机扰动和网络变化")
    print([fig_dir,'network_evolution'],'-depsc')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% I-score
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshold=0.15;
    xt_iscore=11:size(All_sample,2);
    DNB_Iscore=DNB_I_score(All_sample(:,xt_iscore,:),network,threshold);


    % 输出 DNB 结果
    figure
    plot(DNB_Iscore,'.-')
    xlabel("Time")
    ylabel("DNB value")
    title(['Threshold=',num2str(threshold)])

    xt=0:ceil(size(xt_iscore,2)/9):size(xt_iscore,2);
    xt(1)=1;
    xticks(xt)
    xticklabels(xt_iscore(xt))
    
    print([fig_dir,'DNB',num2str(ran)],'-depsc')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 个体熵随时间变化趋势
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A_Individual_Entropy=DNB_Individual_Entropy(All_sample);
    
    %%% 箱图
    figure
    hold on
    boxplot(A_Individual_Entropy')
    xlabel("Time")
    ylabel("Individual Entropy")
    title("个体熵随时间变化的趋势")

    xt=[0:ceil(size(A_Individual_Entropy,1)/4):size(A_Individual_Entropy,1),size(A_Individual_Entropy,1)];
    xt(1)=2;
    xlab=1:size(All_sample,2);
    xlab=xlab(xt);
    xticks(xt)
    xticklabels(xlab)
    
    print([fig_dir,'IE',num2str(ran)],'-depsc')
    
end

save Result0.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 计算DIND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DNB_score_dind=DNB_DIND(Rand_All_sample,8);

figure
plot(1:size(DNB_score_dind,2),DNB_score_dind)
xlabel("Time")
ylabel("DIND")

xt=0:ceil(size(DNB_score_dind,2)/10):size(DNB_score_dind,2);
xt(1)=1;
xlab=1:size(DNB_score_dind,2);
xlab=xlab(xt);

xticks(xt)
xticklabels(xt)
title("DIND随时间变化的趋势")
print(['fig/','DIND'],'-depsc')


toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1

figure
color_order = [rand(sample_num,1),rand(sample_num,1),rand(sample_num,1)];
for g=1:gene
    temp_d=reshape(All_sample_X(g,:,:),size(All_sample_X,2),size(All_sample_X,3));
    subplot(2,5,g)
    hold on
    yyaxis left
    set(gca, 'ColorOrder', color_order);
    ylabel("Expression of Hallmark")
    for i=1:sample_num
        plot(Parameters.dt:Parameters.dt:Parameters.total_t,temp_d(:,i))
    end
    ylim([0,max(max(temp_d))*1.1])

    yyaxis right
    plot(1:size(DNB_score_dind,2),DNB_score_dind)
    ylabel("DNB-DIND")
    title(Hallmarks{g})
    xlabel("Time")
    
end
print(['fig/','DIND_Expression'],'-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2
figure
Box_Data=[];
for i=1:10
    cancer=reshape(Rand_All_sample(1,i,end,:),1,[]);
    nor=reshape(Rand_All_sample(1,i,10,:),1,[]);
    Box_Data=[Box_Data;nor;cancer];
end

h=boxplot(Box_Data','BoxStyle','filled');
colors=['r','b'];
for hi=1:size(h,2)
    set(h(:,hi), 'Color', colors(mod(hi,2)+1)) % 设置每个箱体的颜色
end
xlabel("Normal/Cancer of Hallmark")
ylabel("Expression of Hallmarks")
print(['fig/','Hallmarks_Expression'],'-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 3

[dnb_max_v,dnb_max_i]=max(DNB_score_dind);

All_gene_t2=containers.Map;
for g=1:gene
    nor_data=reshape(All_sample_X(g,floor(nor_t/3):nor_t,:),1,[]);
    nor_max=max(nor_data);
    nor_min=min(nor_data);
    x_sample_t2_id=[];
    for j=1:sample_num
        for t=nor_t+1:time_seq
            if All_sample_X(g,t,j)<nor_min || All_sample_X(g,t,j)>nor_max
                x_sample_t2_id=[x_sample_t2_id,t];
                break;
            end
        end
    end
    x_sample_t2=time_point(x_sample_t2_id)-dnb_max_i;
    All_gene_t2(num2str(g))=x_sample_t2;
end

Dict_boxplot(All_gene_t2)
xlabel("Hallmarks")
ylabel("Difference of times t2-t1")

print(['fig/','Hallmarks_t2t1'],'-depsc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 4











