function [cross_areas, cancer_peak, nor_peak]=Plot_Norm_Cancer_Difference_fill_DF(Hallmarks, All_sample_X,nor_t,cancer_t,figure_type,file_path)
gene=size(All_sample_X,1);

figure('WindowState', 'maximized')
bar_num=100;

cross_areas = zeros(gene, 1); % initial the area
cancer_peak=zeros(gene,1);
nor_peak=zeros(gene,1);

for i=1:gene
    subplot(2,5,i)
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
    ylabel("Probability")
    title(Hallmarks{i})

    cross_areas(i)=distribution_cross_area(nor_xi,nor_pdf,cancer_xi,cancer_pdf);
    cancer_peak(i)=max(cancer_pdf);
    nor_peak(i)=max(nor_pdf);

end
legendHandle=legend("Normal","Cancer",'Location', 'northeastoutside','NumColumns',2);

sgtitle("The distribution of Hallmarks")

set(legendHandle, 'Position', [0.25, 0.01, 0.5, 0.05], 'Units', 'normalized');

print([file_path,'DF_fill'],figure_type)


