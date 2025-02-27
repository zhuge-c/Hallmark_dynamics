% Function to compute KL divergence for genes
function genes_KL = compute_genes_KL(stage_data_1, stage_data_2, adjacent_matrix_2)
    genes_KL = zeros(size(stage_data_1, 1), 1);

    for gene_idx = 1:size(stage_data_1, 1)
        gene_neighbours = find(adjacent_matrix_2(gene_idx, :));

        if numel(gene_neighbours) >= 2
            try
                local_genes_1 = stage_data_1(gene_neighbours, :);
                local_genes_2 = stage_data_2(gene_neighbours, :);

                local_genes_mean_1 = mean(local_genes_1, 2);
                local_genes_mean_2 = mean(local_genes_2, 2);

                if ~isempty(find(std(local_genes_1')==0)) || ~isempty(find(std(local_genes_2')==0))
                    disp("Gene expression is constant")
                end

                local_genes_cov_1 = cov(local_genes_1');
                local_genes_cov_2 = cov(local_genes_2');

                n = size(local_genes_cov_1, 1);
                inv_local_genes_cov_2 = inv(local_genes_cov_2);
                mean_error = (local_genes_mean_2 - local_genes_mean_1)';

                test_val = det(local_genes_cov_2) / det(local_genes_cov_1);
                if test_val > 0
                    local_KL = 0.5 * (trace(inv_local_genes_cov_2 * local_genes_cov_1) + (mean_error * inv_local_genes_cov_2 * mean_error') - n + log(test_val));

                    mean_error_2 = (local_genes_mean_1 - local_genes_mean_2)';
                    local_KL_2 = 0.5 * (trace(inv(local_genes_cov_1) * local_genes_cov_2) + mean_error_2 * inv(local_genes_cov_1) * mean_error_2' - n + log(det(local_genes_cov_1) / det(local_genes_cov_2)));
                    % eq (5) of thesis

                    local_KL = (local_KL + local_KL_2) / 2.0; % eq (6) of thesis

                    if local_KL ~= inf
                        genes_KL(gene_idx) = local_KL / sum(gene_neighbours);
                    end
                end
            catch
                fprintf('Error occurred\n');
            end
        end
    end
end
