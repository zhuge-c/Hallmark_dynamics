function [DNB,Local_KL]=DNB_DIND(data,start_point,reference_num,edge_weights_threshold)

% add the script file path of the functions
addpath('./DIND/');

round_num=size(data,1);
gene_num=size(data,2);
time_num=size(data,3);
sample_num=size(data,4);

n_data=zscore(reshape(data,round_num,gene_num,[]),0,3);
n_data=reshape(n_data,round_num,gene_num,time_num,[]);

round_edge_weights_ts = cell(round_num, 1);

%% iteration of roundnum
for r = 1:round_num
    % Compute edge weights for each time stage and save to CSV
    edge_weights_ts = cell(time_num, 1);

    for t_i = 1:time_num

        % Initialize Pearson correlation coefficient threshold
        initial_pcc_thres = 0;

        % 计算初始相关性网络
        t_data=reshape(n_data(r, :, :, :),size(n_data,2),size(n_data,3),size(n_data,4));

        [stage_data, stage_gene_pccs, ~, stage_node_neighbour_counts] = ...
            compute_stage_initial_graph(t_data, t_i, initial_pcc_thres);
%             stage_data % t_i time slot gene expression
%             stage_gene_pccs % t_i time slot coefficient of correlation between genes
%             stage_node_neighbour_counts % neighbor of each node

        % compute the edge weight ,'rho(x,y|z)'

        edge_weights = compute_edge_weights(stage_data, stage_gene_pccs, ...
            stage_node_neighbour_counts, 0, reference_num);

        edge_weights_ts{t_i} = edge_weights;
    end
    round_edge_weights_ts{r} = edge_weights_ts;
end

% Edge weights threshold
% edge_weights_threshold = 0.25;
round_adjacent_matrixs = cell(round_num, 1);

%% Iterate over rounds
for r = 1:round_num
    adjacent_matrixs = cell(time_num, 1);
    
    % Iterate over time stages
    for t_i = 1:time_num
        edge_weights = round_edge_weights_ts{r}{t_i};
        adjacent_matrix = abs(edge_weights) > edge_weights_threshold;
        diag_index=sub2ind(size(adjacent_matrix),1:size(adjacent_matrix,1),1:size(adjacent_matrix,1));
        adjacent_matrix(diag_index) = true;

        for i = 1:size(adjacent_matrix, 1) - 1
            adjacent_matrix(i+1:end, i) = adjacent_matrix(i, i+1:end);
        end

        adjacent_matrixs{t_i} = adjacent_matrix;

    end

    round_adjacent_matrixs{r} = adjacent_matrixs;
end

% node_expression = reshape(mean(data,4),size(data,1),size(data,2),[]);
% save("DNB_network.mat", "round_adjacent_matrixs","round_edge_weights_ts", "node_expression")

round_max_mean_KLs = zeros(round_num, time_num);
round_states_KL = cell(round_num, 1);
local_KL=zeros(round_num,time_num,gene_num);

%% Iterate over rounds
for r = 1:round_num
    states_KL = cell(time_num, 1);
    max_mean_KLs = zeros(time_num, 1);

    % Iterate over time stages
    for t_idx = start_point:time_num
        normal_state_data = reshape(n_data(r, :, t_idx-1, :),size(n_data,2),size(n_data,4));

        test_state_data = squeeze(n_data(r, :, t_idx, :));

        state_KL = compute_genes_KL(normal_state_data, test_state_data, round_adjacent_matrixs{r}{t_idx});

        state_KL(state_KL < 0) = 0;

        [~, state_max_idxs] = sort(state_KL);

        max_mean_KLs(t_idx) = mean(state_KL(state_max_idxs));

        states_KL{t_idx} = state_KL;
        
        local_KL(r,t_idx,:)=state_KL;
    end

    round_max_mean_KLs(r, :) = max_mean_KLs;

    
    round_states_KL{r} = states_KL;
end
%% Result
Local_KL=reshape(mean(local_KL,1),time_num,gene_num);

DNB=mean(round_max_mean_KLs, 1);


