function cross_areas=Plot_Norm_Cancer_Diff_DF(Hallmarks, All_sample_X, nor_t, cancer_t, gene, figure_type,file_path)

figure('WindowState', 'maximized')

cross_areas = zeros(length(gene), 1); % initial the area

subplot_num=0;

for i=gene
    subplot_num=subplot_num+1;
    subplot(1,2,subplot_num)
    hold on
    nor_data=reshape(All_sample_X(i,floor(nor_t/3):nor_t,1:100),1,[]);
    cancer_data=reshape(All_sample_X(i,floor(nor_t/3+cancer_t):end,1:100),1,[]);
    
    [nor_pdf,nor_xi]=ksdensity(nor_data);
    nor_area=trapz(nor_xi,nor_pdf);

    [cancer_pdf,cancer_xi]=ksdensity(cancer_data);
    cancer_area=trapz(cancer_xi,cancer_pdf);

%     plot(nor_xi, nor_pdf,'k',cancer_xi,cancer_pdf,'r');

    fill([nor_xi fliplr(nor_xi)], [nor_pdf zeros(1, length(nor_pdf))], 'b', 'LineStyle', 'none','FaceAlpha', 0.5);
    fill([cancer_xi fliplr(cancer_xi)], [cancer_pdf zeros(1, length(cancer_pdf))], 'r', 'LineStyle', 'none','FaceAlpha', 0.5);
    
    xlabel("Values")
    ylabel("Frequency")
    title(Hallmarks{i})
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Compute cross area
    %%%%%%%%%%%%%%%%%%%%%%%
    cross_areas(i)=distribution_cross_area(nor_xi,nor_pdf,cancer_xi,cancer_pdf);

%     %%%%%%%%%%%%%%%%%%%%%%
%     % test
%     %%%%%%%%%%%%%%%%%%%%%%
%     [p,h] = signrank(nor_data, cancer_data);

end

legend("nor","cancer")

sgtitle("The distribution of Hallmarks")

print([file_path,'DF'],figure_type)


