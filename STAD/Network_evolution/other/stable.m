clear
clc
tic
%%%%%%%%%%%%%%%%%%%%%%%%
% 运行参数设置
rand_num=1;

t_num=100;

sample_num=1000;

%%%%%%%%%%%%%%%%%%%%%%%
%% 网络数据读取

load ./orginal_data/Hallmarks.mat

% 读取的数据为：第i行为第i个Hallmark对其余Hallmark的调控关系
% 为适应程序，我们需要改成第i行为所有Hallmark对第i个基因的调控
noc=1e4;
[V_orginal,VNT]=xlsread("./orginal_data/normal_Positive.csv");

V_orginal=V_orginal'/noc;

[V_cancer,VCT]=xlsread("./orginal_data/cancer_Positive.csv");

V_cancer=V_cancer'/noc;

network=(V_cancer+V_orginal)/2;

primary_data=xlsread("./orginal_data/HM_normal_expression.csv");

primary_data=primary_data/1e3;

%%%%%%%%%%%%%%%%%%%%%%
%% 参数设置
gene=size(V_cancer,1); % 总的基因个数

Parameters.alpha=ones(gene,gene);

Parameters.lambda=3.8*ones(1,gene); % 基因的最大表达率

Parameters.n=2; % Hill 系数

Parameters.rho1=0.1; % 最小调控促进能力

Parameters.sigma=0.05; % 噪声扰动强度

Parameters.tau=1; % 随机扰动的相关事件

Parameters.theta=mean(sum(V_orginal));

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

cancer_t=floor(time_seq*2/3)+1;

%% 初值赋予
% X0_distribution=reshape(All_sample_X(:,end,:),size(All_sample_X,1),size(All_sample_X,3));
% save("X0_distribution","X0_distribution")

% load X0_distribution.mat

%% 模型求解
load ./orginal_data/stationary_distribution.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All_sample_X 为所有输入样本表达数据,  d1:基因,  d2:时间, d3:样本
% All_sample 为选取时间点数据， 时间点为time_point_id
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All_sample=zeros(gene,t_num,sample_num);
All_sample_X=zeros(gene,time_seq,sample_num);

% stationary_distribution(:,:) = All_sample_X(:,end,:);

%%%%%%%%%%%%%%% 文件夹管理
rand_dir=['./orginal_stable/'];
if exist(rand_dir, 'dir') ~= 7
        mkdir(rand_dir)
end

%% 随机生成 sample_num 个样本
for r=1:sample_num
    clc
    r
    % 初始 eta 与 X0 设置
    eta0=0*rand(1,gene);

    X0 = primary_data;

    % 网络变化
    V1=V_orginal;
    V2=V_orginal;

    %% 模型计算
    [All_X,All_eta]=Node_SDE_modelV2(V1,V2,X0,eta0,Parameters);
    
    mean_all_X=All_X(time_point_id,:);
    for i=1:length(time_point_id)
        mean_all_X(i,:)=mean(All_X(time_point_id(i):time_point_id(i)+9,:));
    end

    All_sample_X(:,:,r)=All_X';


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
print([rand_dir,'network_evolution'],'-depsc')

orginal_distribution=reshape(All_sample_X(:,end,:),size(All_sample_X,1),[]);

save("./orginal_stable/orginal_distribution","orginal_distribution")

