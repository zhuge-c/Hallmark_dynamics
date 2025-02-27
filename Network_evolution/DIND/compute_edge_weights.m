%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is the corresponding DIN build process in DIND
% This is divided into 3 steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function edge_weights = compute_edge_weights(stage_data, stage_gene_pccs,...
    node_neighbour_counts, neighbour_count_thres, reference_count)
    
node_count = size(stage_data, 1);

edge_weights = zeros(node_count);

for node_i = 1:node_count - 1

    %% The common point z of x and y and the corresponding expression data are obtained
    stage_gene_pccs = abs(stage_gene_pccs);

    invalid_bool_idxs = ~ (node_neighbour_counts >= neighbour_count_thres);
    
    % get the index of z
    reference_idxs = get_reference_idxs(node_i, stage_gene_pccs, invalid_bool_idxs, reference_count);

    reference_data = reshape(stage_data(reference_idxs(:), :), size(reference_idxs, 1), size(reference_idxs, 2), []);

    %% (1) compute the weight 'w_x', 'w_y' 
    i_data = stage_data(node_i, :);

    j_data = stage_data(node_i + 1:end, :);

    w_denominator = sum(reference_data.^2, 3);

    % temp_data=zeros(node_count-1,size(reference_data,2));

     %% y
     rho_x=[];
    for j=1:size(reference_data,1)
       %% z
        temp_reference_data=reshape(reference_data(j,:,:),size(reference_data,2),size(reference_data,3));
        temp_z2=w_denominator(j,:);
        rho_x_y_z=[];
        for i=1:size(temp_reference_data,1)
            wx=sum(i_data.*temp_reference_data(i,:))/temp_z2(i);
            ex=i_data-wx.*temp_reference_data(i,:);
            wy=sum(j_data(j,:).*temp_reference_data(i,:))/temp_z2(i);
            ey=j_data(j,:)-wy.*temp_reference_data(i,:);
            n=size(ey,2);
            rho_up=n*sum(ex.*ey)-sum(ex)*sum(ey);
            rho_down=sqrt(n*sum(ex.^2)-sum(ex)^2)*sqrt(n*sum(ey.^2)-sum(ey)^2);
            rho=rho_up/rho_down;
            rho_x_y_z=[rho_x_y_z,rho];
        end
        rho_x_y=mean(rho_x_y_z);
        rho_x=[rho_x,rho_x_y];
    end

    edge_weights(node_i, node_i + 1:end) = rho_x;
end

end