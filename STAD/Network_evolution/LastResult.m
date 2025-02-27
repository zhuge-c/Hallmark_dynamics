%% add the script file path of the functions
clc

addpath('./MainAnalysis/');

figure_type='-djpeg';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Plot distribution
% 2. Compute Area
% 3. Compute peak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[cross_areas, cancer_peak, nor_peak]=Plot_Norm_Cancer_Difference_DF(Hallmarks, All_sample_X,nor_t,cancer_t,figure_type,file_path);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot peak bar
% Plot Area bar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('WindowState', 'maximized')

color=[0 0.4470 0.7410;0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
hold on
yyaxis left
ylabel("Probability")
b1=bar((1:gene)-0.1,[cancer_peak,nor_peak],0.8); 
b1(1).FaceColor=color(1,:);
b1(2).FaceColor=color(2,:);

yyaxis right
b2=bar((1:gene)+0.3,cross_areas,'BarWidth', 0.2);
b2.FaceColor=color(3,:);
ylabel("cross areas")
xlabel("Hallmarks")

legend('Cancer','Normal', 'Cross Area')
print([file_path,'areas-peak'],figure_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot two Hallmarks expression
% Plot two Hallmarks distribution overleap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_gene=[3,4];

Plot_Exp_DNB(Hallmarks,All_sample_X,DNB_score_dind,analysis_gene, Parameters.dt,Parameters.total_t,figure_type,file_path)

Plot_Norm_Cancer_Diff_DF(Hallmarks, All_sample_X,nor_t,cancer_t, analysis_gene, figure_type,file_path);
