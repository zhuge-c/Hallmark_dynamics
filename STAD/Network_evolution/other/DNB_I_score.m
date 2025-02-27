function Inconsis_I=DNB_I_score(All_sample,network,threshold)

% All_sample 为所有输入样本表达数据
% d1:基因
% d2:时间
% d3:样本

% network 为所研究的基因调控网络的拓扑结构

% threshold 为度量边变化差异的阈值

psize=size(All_sample);

%% 读取网络数据
node_num=size(network,1);

adjcent_network={};
for n=1:node_num
    temp={n};
    for a=1:node_num
        if a==n
            continue
        end
        if network(n,a)>0
            temp=[temp,a];
        end
    end
    adjcent_network{n}=temp;
end

%% 计算I-score
time_num=size(All_sample,2);

Inconsis_I=zeros(time_num-2,1);
aver_change=zeros(time_num-1,1);
test_pcc=zeros(time_num,1);

%% 测试变量
A_beta=[];
for t=1:time_num
    for na=1:node_num
        %% 读取节点网络信息
        % 获得所有二阶邻点
        edge_list=[];
        center=adjcent_network{na}{1};
        en=0;
        for n1=2:length(adjcent_network{na})
            nei1=adjcent_network{na}{n1};
            en=en+1;
            edge_list(en,:)=[center,nei1];
            for n2=n1+1:length(adjcent_network{na})
                nei2=adjcent_network{na}{n2};
                if isempty(find(cell2mat(adjcent_network{n1})==nei2))==0
                    en=en+1;
                    edge_list(en,:)=[nei1,nei2];
                end
            end
        end
        %% 计算每条边相关系数
        if t>=2
            diff_edges_num=0;
            for i=1:en
                pre_pcc=abs(corr(reshape(All_sample(edge_list(i,1),t-1,:),psize(3),1),...
                    reshape(All_sample(edge_list(i,2),t-1,:),psize(3),1))); % 前一时刻相关性

                post_pcc=abs(corr(reshape(All_sample(edge_list(i,1),t,:),psize(3),1),...
                    reshape(All_sample(edge_list(i,2),t,:),psize(3),1))); % 当前时刻相关性

                test_edge_pcc(t,i)=post_pcc;

                if abs(pre_pcc-post_pcc)>threshold %% 差异变化较大的边数

                    diff_edges_num=diff_edges_num+1;
                end
                
            end
            test_pcc(t)=test_pcc(t)+mean(test_edge_pcc(t,1:en));

            obs_seq(t-1)=diff_edges_num;

            new_seq(na,1:t-1)=myquantile_6p(obs_seq(1:t-1),en);

            aver_change(t-1)=aver_change(t-1)+new_seq(na,t-1);

        end
    end
    test_pcc(t)=test_pcc(t)/na;
    
    %% 训练隐马尔可夫
    if t>=3
        A_beta=[A_beta,new_seq];
        cell_new_seq=mat2cell(new_seq(:,1:t-2),ones(1,node_num))'; %% 取出多组数据
        init_trans=[0.5,0.5;0.5,0.5];
        for i=1:psize(3)
            init_emiss(1,i)=1/psize(3);
            init_emiss(2,i)=1/psize(3);
        end
        [state_transi,emission]=hmmtrain(cell_new_seq,init_trans,...
            init_emiss,'ALGORITHM','Viterbi');
        %% 计算I-score
        pi=[0.5,0.5];
        aver_pt=0;
        for i=1:node_num
            [beta,pt]=pr_hmm2(new_seq(i,t-2:t-1),state_transi,emission,pi);
            p_at_t=beta(1,2);
            aver_pt=aver_pt+p_at_t/node_num;
        end
        Inconsis_I(t-2)=Inconsis_I(t-2)+(1-aver_pt);
%         t, aver_pt
    end
end
