function All_gene_t2=Plot_DNB_Expression_Time_Difference(DNB_score_dind,Rand_All_sample_X,nor_t,cancer_t,...
    time_point,Hallmarks,figure_type,file_path)

group_num=size(Rand_All_sample_X,1);
gene=size(Rand_All_sample_X,2);
cell_num=size(Rand_All_sample_X,4);
time_seq=size(Rand_All_sample_X,3);

%%%%%%%%%%%%%%%%%%%%%%%%
% find the t1: DNB peak time
%%%%%%%%%%%%%%%%%%%%%%%%
[val,index]=findpeaks(DNB_score_dind,'SortStr','descend');
index2=index(1:2);
val2=val(1:2);
[v,ind]=min(index2);
dnb_max_v=val2(ind);
dnb_max_i=index2(ind);

%%%%%%%%%%%%%%%%%%%%%%%%
% find the t2: Hallmarks expression peak time
%%%%%%%%%%%%%%%%%%%%%%%%

time_windows=5;

All_gene_t2=cell(gene,1);
for g=1:gene
    nor_data=reshape(Rand_All_sample_X(:,g,floor(nor_t/3):nor_t,:),1,[]);
    cancer_data=reshape(Rand_All_sample_X(:,g,floor(cancer_t+nor_t/3):end,:),1,[]);
    nor_max=mean(nor_data)+0.15*(mean(cancer_data)-mean(nor_data));

    x_sample_t2=[];
    for i=1:group_num
        for t=nor_t+1:time_seq
            temp_data=reshape(Rand_All_sample_X(i,g,(t-time_windows):t,:),1,[]);
            if mean(temp_data)>nor_max
                x_sample_t2_id=t;
                break;
            end
        end
        temp_t2=time_point(x_sample_t2_id)-dnb_max_i;
        x_sample_t2=[x_sample_t2,temp_t2];
    end
        All_gene_t2{g}=x_sample_t2;
end

Dict_boxplot(All_gene_t2)

xlabel("Hallmarks")
xticks(1:length(Hallmarks))
xticklabels(Hallmarks)
ylabel("Difference of times t1-t2")

print([file_path,'Hallmarks_t2t1'],figure_type)
