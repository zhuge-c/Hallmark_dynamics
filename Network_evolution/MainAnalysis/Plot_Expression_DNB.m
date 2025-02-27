function Plot_Expression_DNB(Hallmarks,All_sample_X,DNB_score_dind,dt,total_t,figure_type,file_path)


sample_num=size(All_sample_X,3);
time_num=size(All_sample_X,2);
gene=size(All_sample_X,1);

figure('WindowState', 'maximized')
color_order = [rand(sample_num,1),rand(sample_num,1),rand(sample_num,1)];

for g=1:gene
    temp_d=reshape(All_sample_X(g,:,:),time_num,sample_num);
    subplot(2,5,g)
    hold on
    yyaxis left
    set(gca, 'ColorOrder', color_order);
    ylabel("Expression of Hallmark")
    for i=1:sample_num
        plot(dt:dt:total_t,temp_d(:,i))
    end
    ylim([0,max(max(temp_d))*1.1])

    yyaxis right
    plot(1:size(DNB_score_dind,2),DNB_score_dind)
    ylabel("DNB-DIND")
    title(Hallmarks{g})
    xlabel("Time")
end
%% print result
print([file_path,'DIND_Expression'],figure_type)
