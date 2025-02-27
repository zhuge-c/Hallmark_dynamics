function Plot_Exp_DNB(Hallmarks,All_sample_X,DNB_score_dind, gene, dt,total_t,figure_type,file_path)


sample_num=size(All_sample_X,3);
time_num=size(All_sample_X,2);

figure('WindowState', 'maximized')
color_order = [rand(sample_num,1),rand(sample_num,1),rand(sample_num,1)];

subplot_num=0;
for g=gene
    temp_d=reshape(All_sample_X(g,:,:),time_num,sample_num);
    subplot_num=subplot_num+1;
    subplot(1,2,subplot_num)
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yyaxis left
    set(gca, 'ColorOrder', color_order);
    ylabel("Expression of Hallmark")
    for i=1:sample_num
        plot(dt:dt:total_t,temp_d(:,i))
    end
    ylim([0,max(max(temp_d))*1.1])
    ax = gca; % get  the current fig handle
    ax.YColor = 'b';  % set the y axis color

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yyaxis right
    plot(1:size(DNB_score_dind,2),DNB_score_dind,'Color','#EDB120','linewidth',2)

    ylabel("DNB-DIND")
    ax = gca; % get  the current fig handle
    ax.YColor = '#EDB120';  % set y axis color

    ax.XColor = 'k';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    title(Hallmarks{g})
    xlabel("Time")
end
%% print result
print([file_path,'DNB_Exp'],figure_type)
