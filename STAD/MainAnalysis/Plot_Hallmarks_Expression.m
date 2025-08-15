function Plot_Hallmarks_Expression(Hallmarks,All_sample_X,DNB_score_dind,dt,total_t,figure_type,file_path)


sample_num=size(All_sample_X,3);
time_num=size(All_sample_X,2);
gene=size(All_sample_X,1);

% figure('WindowState', 'maximized')
% figure("Position",[10,10,1200,600])
figure

color_order = [rand(sample_num,1),rand(sample_num,1),rand(sample_num,1)];
fontSize = 10;

for g=1:gene
    temp_d=reshape(All_sample_X(g,:,:),time_num,sample_num);
    subplot(2,5,g)
    hold on
    set(gca, 'ColorOrder', color_order);
    ylabel("Hallmark level",'fontsize',fontSize,'fontname','times new roman')
    for i=1:sample_num
        plot(dt:dt:total_t,temp_d(:,i))
    end
%     ylim([0,max(max(temp_d))*1.1])
    ylim([0,1.8])
    title(Hallmarks{g})
    xlabel("Time",'fontsize',fontSize,'fontname','times new roman')

    ax=gca;
    ax.FontSize = fontSize;

    pos1 = get(ax, 'Position'); % 或者 pos1 = ax1.Position;
    if i>5
        new_pos1 = [pos1(1), pos1(2)+pos1(4)*0.2, pos1(3)*0.9, 0.8*pos1(4)];
    else
        new_pos1 = [pos1(1), pos1(2)+pos1(4)*0.1, pos1(3)*0.9, 0.8*pos1(4)];
    end
    set(ax, 'Position', new_pos1); 

end
%% print result

fig = gcf;

fig.PaperUnits = 'centimeters';
width_cm = 21;
height_cm = 10;
fig.PaperSize = [width_cm, height_cm]; 
fig.PaperPosition = [0, 0, width_cm, height_cm]; 

print([file_path,'Hallmark_Expression'],figure_type)
