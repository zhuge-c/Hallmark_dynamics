import scipy.io
import numpy as np

from PlotNetwork import PlotTimeNetwork
from sklearn.preprocessing import StandardScaler


if __name__=="__main__":
    # 加载 .mat 文件
    mat_file_path = 'Result2024-04-21-22-23-10.mat'
    mat_data = scipy.io.loadmat(mat_file_path)

    # nodes weight dimension: gene, time seq, sample
    # edge weight dimension: sample, time seq, gene, gene
    nodes_weight=mat_data['All_sample_X']
    edge_weight=mat_data['All_sample_V']

    # # 创建 StandardScaler 对象
    # scaler = StandardScaler()
    # # 对数据进行标准化
    # nodes_weight = scaler.fit_transform(nodes_weight)
    # edges_weight = scaler.fit_transform(edge_weight)

    time_point_id=mat_data['time_point_id'].flatten().tolist()

    mean_weight=np.mean(edge_weight,axis=0)
    mean_nodes=np.mean(nodes_weight,axis=2)
    mean_nodes=mean_nodes[:,time_point_id]

    node_size=1e6+1
    edge_width=10
    PlotTimeNetwork(mean_nodes,mean_weight,node_size,edge_width,list(range(1,100,10)),2,5,'time_evolution.jpg')



