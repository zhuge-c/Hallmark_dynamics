function orginal_distribution=Stable_State_Compute(gene, Hallmarks, V_orginal, time_seq, sample_num, time_point_id, primary_data, Parameters)

All_sample_X=zeros(gene,time_seq,sample_num);

rand_dir='./orginal_stable/';
if exist(rand_dir, 'dir') ~= 7
        mkdir(rand_dir)
end

for r=1:sample_num
    %% initial 'eta0' and 'X0'
    eta0=0*rand(1,gene);
    X0 = primary_data;

    %% network
    V1=V_orginal;
    V2=V_orginal;

    %% Ä£ÐÍ¼ÆËã
    [All_X,All_eta]=Node_SDE_modelV2(V1,V2,X0,eta0,Parameters);
    
    mean_all_X=All_X(time_point_id,:);
    for i=1:length(time_point_id)
        mean_all_X(i,:)=mean(All_X(time_point_id(i):time_point_id(i)+9,:));
    end
    All_sample_X(:,:,r)=All_X';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the variation with time of the Hallmarks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('WindowState', 'maximized')
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

print([rand_dir,'network_evolution'],'-djpeg')
close('all')
orginal_distribution=reshape(All_sample_X(:,end,:),size(All_sample_X,1),[]);


