% Function to compute stage initial graph
function [stage_data, stage_gene_pccs, stage_initial_graph, stage_node_neighbour_counts] = compute_stage_initial_graph(data, time_stage, initial_pcc_thres)
    
stage_data = reshape(data(:, time_stage, :),size(data,1),[]);

if ~isempty(find(std(stage_data)==0))

    disp("std=0")
    
end

stage_gene_pccs = corr(stage_data');

stage_gene_pccs(isnan(stage_gene_pccs))=0;

stage_initial_graph = zeros(size(stage_gene_pccs));

stage_initial_graph(abs(stage_gene_pccs) >= initial_pcc_thres) = 1.0;

stage_node_neighbour_counts = sum(stage_initial_graph, 2) - 1;

end