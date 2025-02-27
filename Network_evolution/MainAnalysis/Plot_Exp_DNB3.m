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
[maxval,maxindex]=max(val);
dnb_1_v=val2(ind(1));
dnb_1_i=index2(ind(1));
dnb_2_v=val2(ind(2));
dnb_2_i=index2(ind(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title set
t = tiledlayout(2,1,'TileSpacing','Compact');
maxy=max(reshape(All_sample_X(gene,:,:),1,[]))*1.1;

meany=mean(reshape(All_sample_X(gene,:,:),1,[]));

ax1=nexttile;
hold on
set(gca, 'ColorOrder', color_order);
for g=gene
    temp_d=reshape(All_sample_X(g,:,:),time_num,sample_num);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:sample_num
        plot(ax1,dt:dt:total_t,temp_d(:,i),'LineStyle','-','linewidth',0.6)
    end
    % 添加标记文本
    if mean(reshape(temp_d,1,[]))> meany
        text(5, maxy*0.95,Hallmarks{g})
        annotation('textarrow', [0.5 0.5], [0.9 0.85]);
    else
        text(10,maxy*0.05,Hallmarks{g})
        annotation('textarrow', [0.5 0.5], [0.6 0.65]);
    end

end
plot(ax1,dnb_1_i*ones(1,2),[0,maxy],'color','#D9D9D9','LineStyle','--','LineWidth',1)
plot(ax1,dnb_2_i*ones(1,2),[0,maxy],'color','#D9D9D9','LineStyle','--','LineWidth',1)

hold off
xticklabels({})
ylabel(ax1,'Expression of Hallmarks')
ylim([0,maxy])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2=nexttile;
hold on
plot(ax2,1:size(DNB_score_dind,2),DNB_score_dind,'Color','#407BD0','linestyle','-','linewidth',2)
plot(ax2,dnb_1_i*ones(1,2),[0,maxval*1.1],'color','#D9D9D9','LineStyle','--','LineWidth',1)
plot(ax2,dnb_2_i*ones(1,2),[0,maxval*1.1],'color','#D9D9D9','LineStyle','--','LineWidth',1)
ylim([0,maxval*1.1])
hold off
xlabel(t,'Time')
ylabel(ax2,'DNB')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print result
print([file_path,'DNB_Exp2'],figure_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Position',pos)
% title set
t2 = tiledlayout(2,1,'TileSpacing','Compact');
max2y=max(reshape(All_sample_X(gene,:,:),1,[]));

nexttile
hold on

colors2={'#A24D56','#407BD0'};
for j=1:length(gene)
    g=gene(j);
    temp_d2=reshape(All_sample_X(g,:,:),time_num,sample_num);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(dt:dt:total_t,mean(temp_d2,2),'color',colors2{j},'LineWidth',2)

end

plot(dnb_1_i*ones(1,2),[0,max2y],'color','#D9D9D9','LineStyle','--','LineWidth',1)
plot(dnb_2_i*ones(1,2),[0,max2y],'color','#D9D9D9','LineStyle','--','LineWidth',1)
ylabel("Expression of Hallmarks")
hold off
xticklabels({})
ylim([0,max2y])
legend(Hallmarks{gene},'','Location','southeast')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
hold on
plot(1:size(DNB_score_dind,2),DNB_score_dind,'Color','#407BD0','linewidth',2)
plot(dnb_1_i*ones(1,2),[0,maxval*1.1],'color','#D9D9D9','LineStyle','--','LineWidth',1)
plot(dnb_2_i*ones(1,2),[0,maxval*1.1],'color','#D9D9D9','LineStyle','--','LineWidth',1)
hold off
xlabel(t2,'Time')
ylabel('DNB')
ylim([0,maxval*1.1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print result
print([file_path,'DNB_Exp3'],figure_type)
