function time_val = Plot_Interval_Exp_Local_DNB2(Hallmarks,All_sample_X,Local_DNB, ...
    nor_t,cancer_t, time_point, dt,total_t,figure_type,file_path)

fontSize=16;
sample_num=size(All_sample_X,3);
time_num=size(All_sample_X,2);

% figure('WindowState', 'maximized')
figure('Position',[10,10,1000,550])
color_order = [rand(sample_num,1),rand(sample_num,1),rand(sample_num,1)];



%% 
subplot_num=0;

time_val = [];
for g=1:10
    %% compute DNB t1
    DNB_score_dind=Local_DNB(:,g);
    [val,index]=findpeaks(DNB_score_dind,'SortStr','descend');
    index2=index(1:2);
    val2=val(1:2);
    [v,ind]=min(index2);
    dnb_max_v=val2(ind);
    dnb_max_i=index2(ind);
    dnbmax=max(DNB_score_dind);
    time_windows=1;

    %%
    temp_d=reshape(All_sample_X(g,:,:),time_num,sample_num);
    nor_data=reshape(temp_d(floor(nor_t/3):nor_t,:),1,[]);
    cancer_data=reshape(temp_d(floor(cancer_t+nor_t/3):end,:),1,[]);
    nor_max=median(nor_data)+0.2*(median(cancer_data)-median(nor_data));

    time_windows_step=ceil(time_windows/dt);

    for t=nor_t+1:time_num
        temp_data=reshape(temp_d((t-time_windows_step):t,:),1,[]);
        if mean(temp_data)>nor_max
            x_sample_t2_id=t;
            break;
        end
    end


    subplot_num=subplot_num+1;
    subplot(2,5,subplot_num)
    hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % left
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CI=interval_y(temp_d);    
    x=dt:dt:total_t;
    xconf = [x x(end:-1:1)] ;         
    yconf = [CI(:,1)', CI(end:-1:1,2)']; 
    
    yyaxis left

    p = fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');

    plot(x,mean(temp_d,2),'color','#276C9E','LineStyle','-','linewidth',2)
    
    plot(x,ones(1,length(x))*nor_max,'color','#8D4146','linestyle','--','linewidth',1) % plot threshold
    
    meanY=mean(temp_d,2);
    plot(time_point(x_sample_t2_id)*ones(1,2),[0,max(max(temp_d))*1.1],...
        'color','#456990','LineStyle','--','linewidth',1) % plot max time

    ylabel("Expression of Hallmark",'fontsize',fontSize,'fontname','times new roman')
    ylim([0,max(max(temp_d))*1.1])
    ax = gca; % get  the current fig handle
    ax.YColor = '#276C9E';  % set the y axis color

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % right
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    yyaxis right
    plot(1:size(DNB_score_dind,1),DNB_score_dind,'Color','#A24D56','linestyle','-','linewidth',2)
    plot(dnb_max_i*ones(1,2),[0,dnbmax*2.5],'color','#EF767A','LineStyle','--','linewidth',1)

    ylabel("DNB-DIND",'fontsize',fontSize,'fontname','times new roman')
    ylim([0,dnbmax*2.5])
    ax = gca; % get  the current fig handle
    ax.YColor = '#A24D56';  % set y axis color

    ax.XColor = 'k';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    title(Hallmarks{g},'fontsize',fontSize,'fontname','times new roman')
    id=time_point(x_sample_t2_id);
    t2_1_index=dnb_max_i;
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
        xl=[xl,{'t2'},{'t1'}];
        [xt,xt_id]=sort([m_index,id,t2_1_index]);
        xl=xl(xt_id);
    end

    time_val = [time_val; t2_1_index, id];

    xticks(xt)
    xticklabels(xl)
    xlabel("Time",'fontsize',fontSize,'fontname','times new roman')
end
%% print result
print([file_path,'DNB_Exp_interval'],figure_type)

end

function CI=interval_y(y)

flag=3;

if flag==1
    CI=[max(y,[],2),min(y,[],2)];
end

if flag==2
    % compute std and mean
    meanY = mean(y,2);
    stdY = std(y,0,2);
    n = size(y,2);
    
    % compute 95% CI
    se = stdY / sqrt(n);
    CI = meanY +  se* tinv([0.025  0.975], n-1);  % tinv compute t distribution quantle score
end

if flag==3
    % 排序数据
    sorted_data = sort(y,2);
    
    % 计算2.5%和97.5%的分位数位置
    lower_index = ceil(0.025 * size(sorted_data,2));
    upper_index = floor(0.975 * size(sorted_data,2));
    
    % 提取分位数值
    lower_bound = sorted_data(:,lower_index);
    upper_bound = sorted_data(:,upper_index);
    CI=[lower_bound, upper_bound];
end

end