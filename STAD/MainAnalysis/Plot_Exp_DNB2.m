function Plot_Exp_DNB2(Hallmarks,All_sample_X,DNB_score_dind, gene, dt,total_t,figure_type,file_path)


sample_num=size(All_sample_X,3);
time_num=size(All_sample_X,2);

pos=[100,0,400,800 ];
% figure('WindowState', 'maximized')
figure('Position',pos)

color_order = [rand(sample_num,1),rand(sample_num,1),rand(sample_num,1)];

%% Compute DNB time
[val,index]=findpeaks(DNB_score_dind,'SortStr','descend');
index2=index(1:2);
val2=val(1:2);
[v,ind]=sort(index2);
dnb_1_v=val2(ind(1));
dnb_1_i=index2(ind(1));
dnb_2_v=val2(ind(2));
dnb_2_i=index2(ind(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title set
t = tiledlayout(3,1,'TileSpacing','Compact');
maxy=max(reshape(All_sample_X(gene,:,:),1,[]));

for g=gene
    temp_d=reshape(All_sample_X(g,:,:),time_num,sample_num);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    hold on
    set(gca, 'ColorOrder', color_order);
    for i=1:sample_num
        plot(dt:dt:total_t,temp_d(:,i))
    end
    maxy=max(reshape(temp_d,1,[]))*1.1;
    plot(dnb_1_i*ones(1,2),[0,maxy],'k-','LineWidth',1)
    plot(dnb_2_i*ones(1,2),[0,maxy],'k-','LineWidth',1)

    title(Hallmarks{g})
    hold off
    xticklabels({})
    ylim([0,maxy])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
hold on
plot(1:size(DNB_score_dind,2),DNB_score_dind,'Color','#EDB120','linewidth',2)
plot(dnb_1_i*ones(1,2),[0,dnb_1_v],'k-','LineWidth',1)
plot(dnb_2_i*ones(1,2),[0,dnb_2_v],'k-','LineWidth',1)
hold off
xlabel(t,'Time')
ylabel(t,'Expression of Hallmarks')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print result
print([file_path,'DNB_Exp2'],figure_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',pos)
% title set
t2 = tiledlayout(3,1,'TileSpacing','Compact');
max2y=max(reshape(All_sample_X(gene,:,:),1,[]));

for g=gene
    temp_d2=reshape(All_sample_X(g,:,:),time_num,sample_num);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nexttile
    hold on
    set(gca, 'ColorOrder', color_order);

    plot(dt:dt:total_t,mean(temp_d2,2),'LineWidth',2)

    max2y=max(reshape(temp_d2,1,[]))*1.1;
    plot(dnb_1_i*ones(1,2),[0,max2y],'k-','LineWidth',1)
    plot(dnb_2_i*ones(1,2),[0,max2y],'k-','LineWidth',1)

    title(Hallmarks{g})
    hold off
    xticklabels({})
    ylim([0,max2y])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
hold on
plot(1:size(DNB_score_dind,2),DNB_score_dind,'Color','#EDB120','linewidth',2)
plot(dnb_1_i*ones(1,2),[0,dnb_1_v],'k-','LineWidth',1)
plot(dnb_2_i*ones(1,2),[0,dnb_2_v],'k-','LineWidth',1)
hold off
xlabel(t2,'Time')
ylabel(t2,'Expression of Hallmarks')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print result
print([file_path,'DNB_Exp3'],figure_type)
