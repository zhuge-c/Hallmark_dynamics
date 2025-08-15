function Plot_local_DNB(Rand_All_sample,figure_type)

load DNB % get from Compute_DNB

round_num=size(local_dnb,1);

time_num=size(local_dnb{1},1);

gene=size(local_dnb{1}{end},1);

local_DNB=zeros(round_num,time_num,gene);

for i=1:round_num
    for j=1:time_num
        if ~isempty(local_dnb{i}{j})
            local_DNB(i,j,:)=local_dnb{i}{j};
        end
    end
end
D=reshape(mean(local_DNB,1),time_num,gene);

save("local_DNB","D")

bar3(D,0.9)
xticks(1:0.01:10)

print(['fig/','local_dnb'],figure_type)
