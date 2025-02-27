%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% add the script file path of the functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

addpath('./MainAnalysis/');

figure_type='-djpeg';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hallmark expression and DNB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rand_All_sample_X=zeros(group_num,gene,time_seq,cell_num);

for i=1:group_num
    rs=group_id(i,:);
    Rand_All_sample_X(i,:,:,:)=All_sample_X(:,:,rs);
end

Plot_Expression_DNB(Hallmarks,All_sample_X,DNB_score_dind,Parameters.dt,Parameters.total_t,figure_type,file_path)

Plot_Expression_Local_DNB(Hallmarks,All_sample_X,local_dnb,Parameters.dt,Parameters.total_t,figure_type,file_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the specificed point of t1, t2, the expression of Hallmarks and
% DINDs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
All_Hallmarks_t1t2=Plot_DNB_Expression_Time_Difference(DNB_score_dind,Rand_All_sample_X,nor_t,cancer_t,...
    time_point,Hallmarks,figure_type,file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Last Result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare the distribution of Hallmark with normal and cancer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Plot_Norm_Cancer_Difference_bar_and_box(All_sample_X,figure_type,Hallmarks,file_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Plot distribution
% 2. Compute Area
% 3. Compute peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cross_areas, cancer_peak, nor_peak]=Plot_Norm_Cancer_Difference_DF(Hallmarks, All_sample_X,nor_t,cancer_t,figure_type,file_path);

[cross_areas, cancer_peak, nor_peak]=Plot_Norm_Cancer_Difference_fill_DF(Hallmarks, All_sample_X,nor_t,cancer_t,figure_type,file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Area bar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('WindowState', 'maximized')

color=[0 0.4470 0.7410;0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
hold on
b2=bar(1:gene,cross_areas,0.2);
b2.FaceColor=color(3,:);
ylabel("cross areas")
xlabel("Hallmarks")
xticks(1:length(Hallmarks))
xticklabels(Hallmarks)

print([file_path,'areas-peak'],figure_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot two Hallmarks expression
% Plot two Hallmarks distribution overleap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[maxval,maxindex]=max(cross_areas);
[minval,minindex]=min(cross_areas);
analysis_gene=[minindex,maxindex];


Plot_Exp_DNB(Hallmarks,All_sample_X,DNB_score_dind,analysis_gene, Parameters.dt,Parameters.total_t,figure_type,file_path)

Plot_Exp_DNB3(Hallmarks,All_sample_X,DNB_score_dind,analysis_gene, Parameters.dt,Parameters.total_t,figure_type,file_path)

Plot_Interval_Exp_DNB(Hallmarks,All_sample_X(:,:,group_id(1,1:100)),DNB_score_dind,analysis_gene, nor_t,cancer_t, time_point, Parameters.dt,Parameters.total_t,figure_type,file_path)

Plot_Norm_Cancer_Diff_DF(Hallmarks, All_sample_X,nor_t,cancer_t, analysis_gene, figure_type,file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'Resultdata.xlsx';
A = {'Cross area','normal peak', 'cancer peak'};
B=num2cell([cross_areas, nor_peak, cancer_peak]);
A=[A;B];

xlswrite(filename,A,'cross area')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time1_data=reshape(All_sample_X(:,nor_t,:),gene,[])';
time2_data=reshape(All_sample_X(:,cancer_t,:),gene,[])';

save('Ana_Data.mat','time2_data','time1_data')


SignificantAnalysis(All_sample_X, nor_t, cancer_t, figure_type,Hallmarks,file_path)



