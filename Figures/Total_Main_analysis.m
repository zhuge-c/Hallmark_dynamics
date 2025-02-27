%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finaly results analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
addpath('./MainAnalysis/');

%% get all the folders
directory = './DATA/';

folders = dir(directory);

folderNames = {folders([folders.isdir]).name};

cancer_folder=folderNames(3:end);

%% get .mat file path
cancer_mat_file={};
for cancer_ana=cancer_folder
    temp_path=[directory,cancer_ana{:},'/Network_evolution/Result/'];
    
    % get the .mat files
    mat_files = dir(fullfile(temp_path, '*/*.mat'));
    
    % get total name
    mat_fileNames = {mat_files.name};
    filePaths = fullfile({mat_files.folder}, mat_fileNames);
    cancer_mat_file=[cancer_mat_file;filePaths];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get the 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cc=13 % 1:size(cancer_mat_file,1)
    mat_data=cancer_mat_file{cc};
    load(mat_data,'nor_t','cancer_t','All_sample_X','Hallmarks','DNB_score_dind','group_id',...
        'group_num','gene','time_seq','cell_num','time_point','Parameters');
    
    %% compute Rand_All_sample
    Rand_All_sample_X=zeros(group_num,gene,time_seq,cell_num);
    for i=1:group_num
        rs=group_id(i,:);
        Rand_All_sample_X(i,:,:,:)=All_sample_X(:,:,rs);
    end
    result_path=['./Result/',cancer_folder{cc},'/'];

    if ~exist(result_path, 'dir')
        mkdir(result_path);
        disp(['Folder ', result_path, ' is created']);
    else
        disp(['Folder ', result_path, ' is existing']);
    end
    figure_type='-djpeg';

    [cross_areas, cancer_peak, nor_peak]=Plot_Norm_Cancer_Difference_fill_DF(Hallmarks, ...
        All_sample_X,nor_t,cancer_t,figure_type,result_path, [10,10,1000,800]);
    
    [maxval,maxindex]=max(cross_areas);
    [minval,minindex]=min(cross_areas);
    analysis_gene=[minindex,maxindex];
    %% figure 2

    Plot_Exp_DNB4(Hallmarks,All_sample_X,DNB_score_dind,analysis_gene, ...
        Parameters.dt,Parameters.total_t,figure_type,result_path)

    Plot_Norm_Cancer_Diff_DF(Hallmarks, All_sample_X,nor_t,cancer_t, ...
        analysis_gene, figure_type,result_path);

    Plot_Norm_Cancer_Difference_bar_and_box(All_sample_X,figure_type,Hallmarks,result_path)

    %% figure 4
    Plot_Interval_Exp_DNB(Hallmarks,All_sample_X(:,:,group_id(1,1:100)),...
        DNB_score_dind,analysis_gene, nor_t,cancer_t, time_point, ...
        Parameters.dt,Parameters.total_t,figure_type,result_path)

    All_Hallmarks_t1t2=Plot_DNB_Expression_Time_Difference(DNB_score_dind,...
        Rand_All_sample_X,nor_t,cancer_t, time_point,Hallmarks,figure_type,result_path);


    clear nor_t cancer_t All_sample_X Hallmarks DNB_score_dind group_id Rand_All_sample_X
    clear group_num gene time_seq cell_num time_point
end

close all

