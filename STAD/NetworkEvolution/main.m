clear
clc

tic

%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters setting
t_num=100;
sample_num=10000;
figure_type='-djpeg';
%% group setting
group_num=50;
cell_num=1000;
%%%%%%%%%%%%%%%%%%%%%%%
%% Reading network data
load ./orginal_data/Hallmarks.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the load data format:
% i-th row: The regulation relationship of the i-th Hallmark to the other Hallmarks.
%   to adjust the algorithm, we set the i-th row as the regulation of all
%   Hallmarks to the i-th Hallmark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noc=1;

[V_orginal,VNT]=xlsread("./orginal_data/normal_Positive.csv");
V_orginal=V_orginal'/noc;
[V_cancer,VCT]=xlsread("./orginal_data/cancer_Positive.csv");
V_cancer=V_cancer'/noc;
network=(V_cancer+V_orginal)/2;
%%%%%%%%%%%%%%%%%%%%%%
%% Parameters setting
gene=size(V_cancer,1); % the number of all genes
Parameters.alpha=ones(gene,gene)*1e4;
Parameters.lambda=3.8*ones(1,gene); % The maximum rate of gene expression
Parameters.n=2; % Hill coefficient
Parameters.rho1=0.1; % Minimum regulation and promotion capacity
Parameters.sigma=0.05; % Noise disturbance intensity
Parameters.tau=1; % Random disturbance factor
Parameters.theta=mean(sum(V_orginal));
Parameters.dt=0.1; % unit time
Parameters.total_t=100; % the period of time

time_seq=Parameters.total_t/Parameters.dt;
dt=time_seq/t_num;
time_id=1:time_seq;
time_point=0:Parameters.dt:Parameters.total_t;
time_point=time_point(1:end-1);

time_point_id=time_id(1:dt:end);
nor_t=floor(time_seq*1/3);
nor_can_t=floor(time_seq*2/3);
cancer_t=floor(time_seq*2/3)+1;

%% model orginal solution
primary_data=xlsread('./orginal_data/HM_normal_expression.csv');
primary_data=primary_data/1e3;
primary_sample_num=1000;
orginal_distribution=Stable_State_Compute(gene, Hallmarks, V_orginal, time_seq, primary_sample_num, time_point_id, primary_data, Parameters);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All_sample_X are all the expression data of the sample
% d1: gege, d2: time; d3: sample
% All_sample are the selected time slot data, the time points included in 'time_point_id'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All_sample=zeros(gene,t_num,sample_num);
All_sample_X=zeros(gene,time_seq,sample_num);
All_sample_V=zeros(sample_num,t_num,gene,gene);
%% Number of randomly generated samples
for r=1:sample_num
    clc
    r
    %% initial 'eta' and 'X0'
    eta0=0*rand(1,gene);
    X0 = orginal_distribution(:,randsample(size(orginal_distribution,2),1,true))';
    %% normal and cancer network
    V1=V_orginal;
    V2=V_cancer;
    %% model soolution
    [All_X,All_eta,All_v]=Node_SDE_modelV2(V1,V2,X0,eta0,Parameters);
    mean_all_X=All_X(time_point_id,:);
    mean_all_v=All_v(time_point_id,:,:);
    for i=1:length(time_point_id)
        mean_all_X(i,:)=mean(All_X(time_point_id(i):time_point_id(i)+9,:));
        mean_all_v(i,:,:)=mean(All_v(time_point_id(i):time_point_id(i)+9,:,:));
    end
    All_sample(:,:,r)=mean_all_X';
    All_sample_X(:,:,r)=All_X';
    All_sample_V(r,:,:,:)=mean_all_v;

end

%% Group data
group_id=[];
for i=1:group_num
    rs=randperm(sample_num,cell_num);
    group_id=[group_id;rs];
end

%% Result analysis
time_str=datestr(datetime('now'),'yyyy-mm-dd-HH-MM-SS');

file_path=['./Result/fig',time_str,'/'];

if ~exist(file_path, 'dir')
    mkdir(file_path);
    disp(['Folder ' file_path ' is created']);
else
    disp(['Folder ' file_path ' is existing']);
end
figure_type='-djpeg';

toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute DIND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_path = './Result/result/';

Rand_All_sample=zeros(group_num,gene,t_num,cell_num);
for i=1:group_num
    rs=group_id(i,:);
%     rs = randperm(sample_num,cell_num);
    Rand_All_sample(i,:,:,:)=All_sample(:,:,rs);
end

reference_num=6;
edge_weights_threshold=0.25;
addpath('./DIND/');
[DNB_score_dind,local_dnb]=Compute_DNB(Rand_All_sample,reference_num,edge_weights_threshold,figure_type,file_path);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save([file_path,'Result',time_str,'.mat'])

addpath('./MainAnalysis/');

figure_type='-dpng';
% figure 1
Plot_Hallmarks_Expression(Hallmarks,All_sample_X,DNB_score_dind,...
    Parameters.dt,Parameters.total_t,figure_type,file_path)










