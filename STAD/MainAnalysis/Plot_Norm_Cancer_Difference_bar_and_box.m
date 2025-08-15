function [hist_data, Box_Data]=Plot_Norm_Cancer_Difference_bar_and_box(All_sample_X,figure_type,Hallmarks,file_path)

gene=size(All_sample_X,1);


Box_Data=zeros(1,size(All_sample_X,2));
sn=size(Box_Data,2);
for i=1:gene
    cancer=reshape(All_sample_X(i,end,:),1,[]);
    cancer=Del_out_Box_Plot(cancer);
    nor=reshape(All_sample_X(i,10,:),1,[]);
    nor=Del_out_Box_Plot(nor);
    sn=min([size(cancer,2),size(nor,2),sn]);
    Box_Data=[Box_Data(:,1:sn);nor(1:sn);cancer(1:sn)];
end

Box_Data(1,:)=[];

mean_Data=mean(Box_Data,2);
hist_data=[mean_Data(1:2:end),mean_Data(2:2:end)];
%% bar plot
figure('WindowState', 'maximized')
bar(hist_data)
legend("Normal","Cancer")
xlabel("Hallmark")
xticks(1:length(Hallmarks))
xticklabels(Hallmarks)
ylabel("Expression of Hallmarks")
print([file_path,'Hallmarks_Expression_bar'],figure_type)

%% box_plot
width=0.2;
rp=2/3;
position=zeros(1,size(Box_Data,1));
position(1:2:end)=(1:10)-width*rp;
position(2:2:end)=(1:10)+width*rp;

figure('WindowState', 'maximized')
h=boxplot(Box_Data','Symbol','o','OutlierSize',3,'Colors',[0,0,0],...
    'positions',position,'Widths',width);

% h=boxplot(Box_Data','MedianStyle', 'target','PlotStyle','compact','positions',[1 1.1 2 2.1 3 3.1 4 4.1 5 5.1 6 6.1 7 7.1 8 8.1 9 9.1 10 10.1 ]);

colors=['r','b'];
% colors={'#0072BD', '#A2142F'};

% set the line width
lineObj=findobj(gca,'Type','Line');
for i=1:length(lineObj)
    lineObj(i).LineWidth=1;
    lineObj(i).MarkerFaceColor=[1,1,1].*.3;
    lineObj(i).MarkerEdgeColor=[1,1,1].*.3;
end

% set the color of the box 
boxObj=findobj(gca,'Tag','Box');
for i=1:length(boxObj)
    patch(boxObj(i).XData,boxObj(i).YData,colors(mod(i,2)+1),'FaceAlpha',0.5,...
        'LineWidth',1.1);
end
% for hi=1:size(h,2)
%     set(h(:,hi), 'Color', colors{mod(hi,2)+1}) % set the color
% end


xlabel("Normal/Cancer of Hallmark")
xticks(1:length(Hallmarks))
xticklabels(Hallmarks)
ylabel("Expression of Hallmarks")
legend("Cancer","Normal",'fontsize',16,'fontname','Times New Roman','location','northeast')

print([file_path,'Hallmarks_Expression_box'],figure_type)




