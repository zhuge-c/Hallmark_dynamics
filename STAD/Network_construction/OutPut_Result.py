import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import numpy as np
import pickle

def GraphData_Process(Graph_Data, Num_time):
    with open(Graph_Data, "rb") as f:
        loaded_data = pickle.load(f)
        Graph = loaded_data["Graph"]

    for key in Graph.keys():
        edges = Graph[key]
        for edge in edges.keys():
            Graph[key][edge] = [x / Num_time for x in edges[edge]]

    pickle_data = {"Graph": Graph}

    with open(Graph_Data, "wb") as f:
        pickle.dump(pickle_data, f)



def Print_Graph_Adj(Graph_Data,Out_Graph_Adj):
    # 从文件中读取多个变量
    with open(Graph_Data, "rb") as f:
        loaded_data = pickle.load(f)
        Graph=loaded_data["Graph"]

    ################## 将数据输出为邻接表
    Graph_adj=pd.DataFrame(Graph)

    Graph_adj.to_csv(Out_Graph_Adj)

def Print_Negative_Positive_Adj(Graph_Data,Out_Negative_net, Out_Positive_net):
    # 从文件中读取多个变量
    with open(Graph_Data, "rb") as f:
        loaded_data = pickle.load(f)
        Graph = loaded_data["Graph"]


    ################################
    Graph_positive = {}
    Graph_negative = {}
    for key in Graph.keys():
        temp_edge = Graph[key]
        temp_n = {}
        temp_p = {}
        for edge in temp_edge.keys():
            temp_n[edge] = temp_edge[edge][1]
            temp_p[edge] = temp_edge[edge][0]
        Graph_negative[key] = temp_n
        Graph_positive[key] = temp_p
    Graph_positive = pd.DataFrame(Graph_positive)
    Graph_negative = pd.DataFrame(Graph_negative)
    ##### 获得数据均值
    Median_N = np.median(Graph_negative.values.ravel())
    Median_P = np.median(Graph_positive.values.ravel())

    Graph_negative.to_csv(Out_Negative_net)
    Graph_positive.to_csv(Out_Positive_net)

######################################
# 绘制网络
######################################
def Plot_NetWork(Graph_Data,Out_figure):
    # 从文件中读取多个变量
    with open(Graph_Data, "rb") as f:
        loaded_data = pickle.load(f)
        Graph = loaded_data["Graph"]
    G = nx.DiGraph()


    for hm in Graph.keys():
        temp_edge=Graph[hm]
        for eg in temp_edge.keys():
            G.add_edge(hm,eg,weight1=round(temp_edge[eg][0],5),weight2=round(temp_edge[eg][1],5))

    # 计算节点的布局
    pos = nx.spring_layout(G)

    plt.figure(figsize=(24,16)) # 长，高
    # 绘制图形
    nx.draw(G, pos, with_labels=True, node_size=1000, node_color='skyblue', font_size=12, font_color='black', edge_color='gray') # 绘制顶点并给与标签

    # 获取边的属性字典，并添加到边的标签中
    edge_labels = {}
    for u, v, attrs in G.edges(data=True):
        weight1 = attrs.get('weight1', '')
        weight2 = attrs.get('weight2', '')
        label = f'({weight1}, {weight2})'
        edge_labels[(u, v)] = label

    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels) # 绘制边标签

    plt.title('Graph with Multiple Weighted Edges')

    plt.show(block=False)

    plt.savefig(Out_figure)

