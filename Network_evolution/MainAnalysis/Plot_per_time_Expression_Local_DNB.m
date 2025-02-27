function Plot_per_time_Expression_Local_DNB(Hallmarks,All_sample_X,DNB_score_dind,dt,total_t,figure_type,...
    Rand_All_sample_X,group_num,cell_num,nor_t,cancer_t,time_point,file_path)

x_position=-0.08;
sample_num=size(All_sample_X,3);
time_num=size(All_sample_X,2);
gene=size(All_sample_X,1);

for g=1:gene
    figure('WindowState', 'maximized')
    % cell
    
    per_cell=reshape(Rand_All_sample_X(1,g,:,1),time_num,[]);
    % find the mutation point
    nor_data=reshape(Rand_All_sample_X(1,g,floor(nor_t/3):nor_t,1),1,[]);
    cancer_data=reshape(Rand_All_sample_X(1,g,floor(cancer_t+nor_t/3):end,1),1,[]);
    
    nor_max=mean(nor_data)+0.15*(mean(cancer_data)-mean(nor_data));

    for t=nor_t+1:time_num
        if per_cell(t)>nor_max
            id=time_point(t);
            value=per_cell(t);
            break;
        end
    end
    subplot(3,1,1)
    hold on
    plot(dt:dt:total_t,per_cell)
    plot(id*ones(1,2),[0,value],'k--')
    yyaxis left
    ylim([min(per_cell)*0.9,max(per_cell)*1.1])
    ylabel("Expression of Hallmark")

    yyaxis right
    plot(1:size(DNB_score_dind(:,g),1),DNB_score_dind(:,g))
    [t2_1,t2_1_index]=max(DNB_score_dind(:,g));
    plot(t2_1_index*ones(1,2),[0,t2_1],'k--')
    ylabel("DNB-DIND")
    
    if id==t2_1_index
        m_index=[0,10,20,50:10:100];
        log_ind=m_index==id;
        m_index(log_ind)=[];
        xl=arrayfun(@num2str, m_index, 'UniformOutput', false);
        xl=[xl,{'t1/t2'}];
        [xt,xt_id]=sort([m_index,id]);
        xl=xl(xt_id);
    else
        m_index=[0,10,20,50:10:100];
        log_ind=m_index==id;
        m_index(log_ind)=[];
        log_ind=m_index==t2_1_index;
        m_index(log_ind)=[];
        xl=arrayfun(@num2str, m_index, 'UniformOutput', false);
        xl=[xl,{'t1'},{'t2'}];
        [xt,xt_id]=sort([m_index,id,t2_1_index]);
        xl=xl(xt_id);
    end

    xticks(xt)
    xticklabels(xl)
    xlabel("Time")
    text('String', 'Cell', 'Units', 'normalized', 'Position', [x_position, 0.5], 'Rotation', 90, 'HorizontalAlignment', 'center');


    % human
    per_human=reshape(Rand_All_sample_X(1,g,:,:),time_num,cell_num);
    per_human=mean(per_human,2);

    % find the mutation point
    nor_data=reshape(Rand_All_sample_X(1,g,floor(nor_t/3):nor_t,:),1,[]);
    cancer_data=reshape(Rand_All_sample_X(1,g,floor(cancer_t+nor_t/3):end,:),1,[]);
    
    nor_max=mean(nor_data)+0.15*(mean(cancer_data)-mean(nor_data));

    for t=nor_t+1:time_num
        if per_human(t)>nor_max
            id=time_point(t);
            value=per_human(t);
            break;
        end
    end
    subplot(3,1,2)
    hold on
    plot(dt:dt:total_t,per_human)
    plot(id*ones(1,2),[0,value],'k--')
    yyaxis left
    ylim([min(per_human)*0.9,max(per_human)*1.1])
    ylabel("Expression of Hallmark")

    yyaxis right
    plot(1:size(DNB_score_dind(:,g),1),DNB_score_dind(:,g))
    [t2_1,t2_1_index]=max(DNB_score_dind(:,g));
    plot(t2_1_index*ones(1,2),[0,t2_1],'k--')
    ylabel("DNB-DIND")


    
    if id==t2_1_index
        m_index=[0,10,20,50:10:100];
        log_ind=m_index==id;
        m_index(log_ind)=[];
        xl=arrayfun(@num2str, m_index, 'UniformOutput', false);
        xl=[xl,{'t1/t2'}];
        [xt,xt_id]=sort([m_index,id]);
        xl=xl(xt_id);
    else
        m_index=[0,10,20,50:10:100];
        log_ind=m_index==id;
        m_index(log_ind)=[];
        log_ind=m_index==t2_1_index;
        m_index(log_ind)=[];
        xl=arrayfun(@num2str, m_index, 'UniformOutput', false);
        xl=[xl,{'t1'},{'t2'}];
        [xt,xt_id]=sort([m_index,id,t2_1_index]);
        xl=xl(xt_id);
    end

    xticks(xt)
    xticklabels(xl)
    xlabel("Time")
    text('String', 'Human', 'Units', 'normalized', 'Position', [x_position, 0.5], 'Rotation', 90, 'HorizontalAlignment', 'center');


    % group
    per_group=reshape(Rand_All_sample_X(:,g,:,:),group_num,time_num,cell_num);
    per_group=mean(per_group,[1,3]);

     % find the mutation point
    nor_data=reshape(Rand_All_sample_X(:,g,floor(nor_t/3):nor_t,:),1,[]);
    cancer_data=reshape(Rand_All_sample_X(:,g,floor(cancer_t+nor_t/3):end,:),1,[]);
    
    nor_max=mean(nor_data)+0.15*(mean(cancer_data)-mean(nor_data));

    for t=nor_t+1:time_num
        if per_group(t)>nor_max
            id=time_point(t);
            value=per_group(t);
            break;
        end
    end

    subplot(3,1,3)
    hold on
    plot(dt:dt:total_t,per_group)
    plot(id*ones(1,2),[0,value],'k--')
    yyaxis left
    ylim([min(per_group)*0.9,max(per_group)*1.1])
    ylabel("Expression of Hallmark")

    yyaxis right
    plot(1:size(DNB_score_dind(:,g),1),DNB_score_dind(:,g))
    [t2_1,t2_1_index]=max(DNB_score_dind(:,g));
    plot(t2_1_index*ones(1,2),[0,t2_1],'k--')
    ylabel("DNB-DIND")
    

    if id==t2_1_index
        m_index=[0,10,20,50:10:100];
        log_ind=m_index==id;
        m_index(log_ind)=[];
        xl=arrayfun(@num2str, m_index, 'UniformOutput', false);
        xl=[xl,{'t1/t2'}];
        [xt,xt_id]=sort([m_index,id]);
        xl=xl(xt_id);
    else
        m_index=[0,10,20,50:10:100];
        log_ind=m_index==id;
        m_index(log_ind)=[];
        log_ind=m_index==t2_1_index;
        m_index(log_ind)=[];
        xl=arrayfun(@num2str, m_index, 'UniformOutput', false);
        xl=[xl,{'t1'},{'t2'}];
        [xt,xt_id]=sort([m_index,id,t2_1_index]);
        xl=xl(xt_id);
    end
    
    xticks(xt)
    xticklabels(xl)
    xlabel("Time")
    text('String', 'Group', 'Units', 'normalized', 'Position', [x_position, 0.5], 'Rotation', 90, 'HorizontalAlignment', 'center');

    sgtitle(Hallmarks{g})

print([file_path,'Per_local_DIND_Expression',Hallmarks{g}],figure_type)
end

