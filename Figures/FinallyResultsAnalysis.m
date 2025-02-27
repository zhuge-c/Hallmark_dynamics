%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finaly results analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
addpath('./MainAnalysis/');

%% get all the folders
directory = './DATA/';

%  dir: get all the files under directory
folders = dir(directory);

% screen the folder 
folderNames = {folders([folders.isdir]).name};

% cancer folder
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
Cancer_Hallmarks_Expression=[];
All_Hallmarks_Expression=[];
Cancer_overleap=[];
Cancer_Time_diference=[];
for mat_data=cancer_mat_file'
    load(mat_data{:},'nor_t','cancer_t','All_sample_X','Hallmarks','DNB_score_dind','group_id',...
        'group_num','gene','time_seq','cell_num','time_point');
    % get hallmarks expression
    mean_all_samp_x=mean(All_sample_X,3);
    temp_hallmark_expression=[mean_all_samp_x(:,nor_t);mean_all_samp_x(:,cancer_t)]';
    Cancer_Hallmarks_Expression=[Cancer_Hallmarks_Expression;temp_hallmark_expression];

    temp_nor=reshape(All_sample_X(:,nor_t,:),size(All_sample_X,1),[]);
    temp_cancer=reshape(All_sample_X(:,cancer_t,:),size(All_sample_X,1),[]);
    temp_all_harmarks=[temp_nor;temp_cancer]';
    All_Hallmarks_Expression=[All_Hallmarks_Expression;temp_all_harmarks];

    % get  overleap 
    [cross_areas, cancer_peak, nor_peak]=Norm_Cancer_Difference_DF(Hallmarks,All_sample_X,nor_t,cancer_t);
    Cancer_overleap=[Cancer_overleap;cross_areas'];

    % get time diference
    Rand_All_sample_X=zeros(group_num,gene,time_seq,cell_num);
    
    for i=1:group_num
        rs=group_id(i,:);
        Rand_All_sample_X(i,:,:,:)=All_sample_X(:,:,rs);
    end
    All_Hallmarks_t1t2=DNB_Expression_Time_Difference(DNB_score_dind,Rand_All_sample_X,nor_t,cancer_t,...
        time_point,Hallmarks);
    All_Hallmarks_t1t2=cellfun(@(x) mean(x),All_Hallmarks_t1t2);
    Cancer_Time_diference=[Cancer_Time_diference;All_Hallmarks_t1t2'];

    clear nor_t cancer_t All_sample_X Hallmarks DNB_score_dind group_id Rand_All_sample_X
    clear group_num gene time_seq cell_num time_point
end

close all

result_path='./Result/';
save([result_path,'All_Expression.mat'],"All_Hallmarks_Expression")

save([result_path,'cancer_file'],'cancer_folder','cancer_mat_file')

cancer_expression_data=table(cancer_folder',Cancer_Hallmarks_Expression);
writetable(cancer_expression_data,[result_path,'cancer_expression.csv'],'Encoding', 'GBK')

cancer_overleap_data=table(cancer_folder',Cancer_overleap);
writetable(cancer_overleap_data,[result_path,'cancer_overleap.csv'],'Encoding', 'GBK')

cancer_time_diference_data=table(cancer_folder',Cancer_Time_diference);
writetable(cancer_time_diference_data,[result_path,'cancer_time_t1t2.csv'],'Encoding', 'GBK')

