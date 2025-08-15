function Plot_Exp_DNB4(Hallmarks,All_sample_X,DNB_score_dind, gene, dt,total_t,figure_type,file_path)


sample_num=size(All_sample_X,3);
time_num=size(All_sample_X,2);

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

maxy=max(reshape(All_sample_X(gene,:,:),1,[]))*1.1;

meany=mean(reshape(All_sample_X(gene,:,:),1,[]));

pos=[100,0,500,400 ];
% figure('WindowState', 'maximized')

figure('Position',pos)
fontSize=16;

hold on
set(gca, 'ColorOrder', color_order);
for g=gene
    temp_d=reshape(All_sample_X(g,:,:),time_num,sample_num);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:sample_num
        plot(dt:dt:total_t,temp_d(:,i),'LineStyle','-','linewidth',0.6)
    end
    % 添加标记文本
    if mean(reshape(temp_d,1,[]))> meany
        text(5, maxy*0.95,Hallmarks{g},'fontsize',fontSize,'fontname','times new roman')
        annotation('textarrow', [0.4 0.5], [0.8 0.75]);
    else
        text(10,maxy*0.05,Hallmarks{g},'fontsize',fontSize,'fontname','times new roman')
        annotation('textarrow', [0.5 0.5], [0.2 0.3]);
    end

end
plot(dnb_1_i*ones(1,2),[0,maxy],'color','#D9D9D9','LineStyle','--','LineWidth',1)
plot(dnb_2_i*ones(1,2),[0,maxy],'color','#D9D9D9','LineStyle','--','LineWidth',1)

hold off
ylabel('Expression of Hallmarks','fontsize',fontSize,'fontname','times new roman')
ylim([0,maxy])
xlabel('Time','fontsize',fontSize,'fontname','times new roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print result
print([file_path,'DNB_Exp4'],figure_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% title set

max2y=max(reshape(All_sample_X(gene,:,:),1,[]));

figure('Position',[100,0,500,700])
hold on

colors2={'#A24D56','#407BD0'};
yyaxis left
for j=1:length(gene)
    g=gene(j);
    temp_d2=reshape(All_sample_X(g,:,:),time_num,sample_num);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot(dt:dt:total_t,mean(temp_d2,2),'color',colors2{j},'LineStyle','-','LineWidth',2)

end

ylabel("Expression of Hallmarks")
ylim([0,max2y])

yyaxis right

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DNB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(1:size(DNB_score_dind,2),DNB_score_dind,'Color','#ED7C72','linewidth',2)

plot(dnb_1_i*ones(1,2),[0,maxval*5],'color','#D9D9D9','LineStyle','--','LineWidth',1)
plot(dnb_2_i*ones(1,2),[0,maxval*5],'color','#D9D9D9','LineStyle','--','LineWidth',1)
hold off
xlabel('Time','fontsize',fontSize,'fontname','times new roman')
ylabel('DNB','fontsize',fontSize,'fontname','times new roman')
ylim([0,maxval*2.5])
legend({Hallmarks{gene},'DNB','',''},'Location','southoutside','fontsize',fontSize,'fontname','times new roman')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print result
print([file_path,'DNB_Exp4'],figure_type)
