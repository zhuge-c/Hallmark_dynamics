%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute DIND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DNB_score_dind,local_dnb]=Compute_DNB(Rand_All_sample,reference_num,edge_weights_threshold, figure_type,file_path)

[DNB_score_dind,local_dnb]=DNB_DIND(Rand_All_sample,8,reference_num,edge_weights_threshold);

figure
plot(1:size(DNB_score_dind,2),DNB_score_dind)
xlabel("Time")
ylabel("DIND")

xt=0:ceil(size(DNB_score_dind,2)/10):size(DNB_score_dind,2);
xt(1)=1;
xlab=1:size(DNB_score_dind,2);
xlab=xlab(xt);

xticks(xt)
xticklabels(xt)
title('DIND trends over time')

print([file_path,'DIND'],figure_type)

% save("DNB","DNB_score_dind","local_dnb")