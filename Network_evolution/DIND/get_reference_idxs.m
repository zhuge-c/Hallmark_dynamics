function reference_idxs = get_reference_idxs(a_idx, stage_gene_pccs, invalid_node_bool_idxs, reference_count)
    node_count = size(stage_gene_pccs, 1);
    assert(a_idx >= 0 && a_idx < node_count );

    % compute the pcc
    a_pcc = stage_gene_pccs(a_idx, :);
    b_pccs = stage_gene_pccs(a_idx + 1:end, :);

    mean_pcc = (a_pcc + b_pccs) / 2.0;

    mean_pcc(:, invalid_node_bool_idxs) = -1.0;
    mean_pcc(:, a_idx ) = -1.0;

    % self ind
    self_ind=sub2ind(size(mean_pcc),1:node_count - a_idx, a_idx + 1:node_count);
    
    mean_pcc(self_ind) = -1.0;

    [~, reference_idxs] = sort(mean_pcc, 2, 'descend');
    
    reference_idxs = reference_idxs(:, 1:reference_count);

end
