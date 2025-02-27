function SignificantAnalysis(All_sample_X, nor_t, cancer_t, figure_type,Hallmarks,file_path)

gene=size(All_sample_X,1);

Box_Data=zeros(1,size(All_sample_X,2));
sn=size(Box_Data,2);

for i=1:gene
    cancer=reshape(All_sample_X(i,cancer_t,:),1,[]);
    cancer=Del_out_Box_Plot(cancer);
    nor=reshape(All_sample_X(i,nor_t,:),1,[]);
    nor=Del_out_Box_Plot(nor);
    sn=min([size(cancer,2),size(nor,2),sn]);
    Box_Data=[Box_Data(:,1:sn);nor(1:sn);cancer(1:sn)];
end

Box_Data(1,:)=[];

p_value=zeros(10,1);
for i=1:10
    nor=Box_Data(2*i-1,:);
    cancer=Box_Data(2*i,:);
    [h,p]=ttest(nor,cancer);
    p_value(i)=p;
end

%% box_plot
width=0.2;
rp=2/3;
position=zeros(1,size(Box_Data,1));
position(1:2:end)=(1:10)-width*rp;
position(2:2:end)=(1:10)+width*rp;

figure('WindowState', 'maximized')
h=boxplot(Box_Data','Symbol','o','OutlierSize',3,'Colors',[0,0,0],...
    'positions',position,'Widths',width);
colors=['r','b'];

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


xlabel("Normal/Cancer of Hallmark")
xticks(1:length(Hallmarks))
xticklabels(Hallmarks)
ylabel("Expression of Hallmarks")
legend("Cancer","Normal",'fontsize',16,'fontname','Times New Roman','location','northeast')
