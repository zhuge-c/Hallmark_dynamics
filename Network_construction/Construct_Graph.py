# 存储数据
import pickle

def Construct_Graph(HallmarksData,Graph_Data):
    # 从文件中读取多个变量
    with open(HallmarksData, "rb") as f:
        loaded_data = pickle.load(f)
        Hallmarks_TFS=loaded_data["Hallmarks_TFS"]
        Hallmarks_Gene=loaded_data["Hallmarks_Gene"]
        network=loaded_data["network"]

    Graph={}
    for hallmarks in Hallmarks_TFS.keys():
        temp_key=hallmarks
        temp_edge={}
        temp_tfs=Hallmarks_TFS[hallmarks] # 矩阵行因子
        # 提取子矩阵
        interset_tfs = list(set(network.index).intersection(set(temp_tfs)))
        matrix = network.loc[interset_tfs]
        for tedge in Hallmarks_TFS.keys():
            edge_gene=Hallmarks_Gene[tedge] # 矩阵列因子
            col_matrix=matrix[edge_gene] # 获取基因矩阵
            positive_weight=col_matrix[col_matrix>0].sum().sum()
            negative_weight=col_matrix[col_matrix<=0].sum().sum()
            temp_edge[tedge]=[positive_weight,negative_weight]
        Graph[temp_key]=temp_edge

    pickle_data={"Graph":Graph}

    with open(Graph_Data, "wb") as f:
        pickle.dump(pickle_data, f)

