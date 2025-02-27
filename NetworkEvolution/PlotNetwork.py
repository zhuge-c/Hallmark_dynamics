import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

node_size_stand=100
edge_width_stand=40
edge_times=8
node_times=1

def PlotNetWork(G,nodemin,nodemax,edgemin,edgemax,ax,fig_title):
    # G is a graph, contain nodes, nodes weight, edges, edge weight
    # ax is the canvas

    pos = nx.circular_layout(G)  # Position nodes in a circle

    # Extract node and edge attributes for visual styles
    node_weights = np.array([G.nodes[n]['weight'] for n in G.nodes])
    edge_weights = np.array([G.edges[e]['weight'] for e in G.edges])

    # 可以使用 np.maximum 函数将所有小于等于零的值替换为一个很小的正数（例如，1e-10）
    temp_edges_log = (np.log(np.maximum(edge_weights, 1e-10))+ edge_width_stand)/edge_times
    temp_node_log = (np.log(np.maximum(node_weights, 1e-10))+node_size_stand)/node_times

    # Draw the graph
    nodes = nx.draw_networkx_nodes(G, pos, node_size=temp_node_log,
                                   node_color=node_weights, vmin=nodemin, vmax=nodemax,
                                   cmap=plt.cm.viridis, ax=ax)
    edges = nx.draw_networkx_edges(G, pos, ax=ax,
                                   edge_color=edge_weights, edge_cmap=plt.cm.plasma,
                                   edge_vmin=edgemin, edge_vmax=edgemax,
                                   width=temp_edges_log)

    # 将节点标签位置稍微偏离原节点位置
    scale_size=1.13
    pos = {key: [value[0]*scale_size, value[1]*scale_size] for key, value in pos.items()}  # 增加y坐标

    nx.draw_networkx_labels(G, pos, labels={node: node for node in G.nodes}, ax=ax,font_size=14)

    # Set the title for each subplot
    ax.set_title(fig_title,fontsize=16)
    ax.axis('off')

    return nodes, edges

def Create_Graph(nodes_weight,edge_weight):
    # nodes_weight
    # edge_weight
    G = nx.Graph()
    num_nodes=len(nodes_weight)
    for i in range(1,num_nodes+1):
        G.add_node(i, weight=nodes_weight[i-1])
        for j in range(1,num_nodes+1):
            if i==j:
                continue
            G.add_edge(i, j, weight=edge_weight[i-1,j-1])
    return G

def PlotTimeNetwork(nodes_weight,edges_weight,node_weight_time,edge_weight_time,times_point,fig_row,fig_col,Out_figure):
    fig, axes = plt.subplots(fig_row, fig_col, figsize=(16, 10), subplot_kw={'aspect': 'auto'})
    # figsize = (20, 5),
    G_list=[]
    nodemin, edgemin=np.inf, np.inf
    nodemax, edgemax=-np.inf, -np.inf

    for i in range(len(times_point)):
        t=times_point[i]
        if t<3:
            time_windows=list(range(0,t+3))
        elif t>nodes_weight.shape[1]:
            time_windows=list(range(t-3,nodes_weight.shape[1]))
        else:
            time_windows = list(range(t-3,t+3))
        temp_nodes=np.var(nodes_weight[:,time_windows],axis=1)
        temp_edges=np.var(edges_weight[time_windows,:,:],axis=0)
        nodemin=min(nodemin,np.min(temp_nodes))
        nodemax=max(nodemax,np.max(temp_nodes))
        edgemin=min(edgemin,np.min(temp_edges))
        edgemax=max(edgemax,np.max(temp_edges))
        G=Create_Graph(temp_nodes,
                       temp_edges)
        G_list.append(G)
        # print(f'{i//fig_col},{i%fig_col}')

    for i in range(len(times_point)):
        G=G_list[i]
        ax=axes[i//fig_col,i%fig_col]
        nodes, edges = PlotNetWork(G, nodemin,nodemax,edgemin,edgemax, ax, f't={times_point[i]}')

    # Layout adjustments for subplot to prevent overlap
    plt.subplots_adjust(bottom=0.2, top=0.8, left=0.05, right=0.90)

    # Add color bar for nodes
    cbar_pos = fig.add_axes([0.05, 0.15, 0.9, 0.03])  # x, y, width, height
    cbar = fig.colorbar(nodes, cax=cbar_pos, orientation='horizontal')
    cbar.set_label('Hallmark expression level',fontsize=16)

    # Add color bar for edges
    cbar_pos_edges = fig.add_axes([0.05, 0.05, 0.9, 0.03])  # x, y, width, height
    cbar_edges = fig.colorbar(edges, cax=cbar_pos_edges, orientation='horizontal')
    cbar_edges.set_label('Hallmark Interaction strength',fontsize=16)

    # Show the plot
    plt.show()
    plt.savefig(Out_figure)



if __name__=="__main__":
    # Number of nodes
    num_nodes = 10
    # Number of time points
    time_points = 5

    # Initialize a list of graphs
    graphs = []

    # Create random graph data for each time point
    for _ in range(time_points):
        G = nx.Graph()
        for i in range(1, num_nodes + 1):
            G.add_node(i, weight=np.random.rand())
        for i in range(1, num_nodes):
            G.add_edge(i, i + 1, weight=np.random.rand())
        graphs.append(G)

    # Set up the matplotlib figure and axes
    fig, axes = plt.subplots(1, time_points, figsize=(20, 5), subplot_kw={'aspect': 'auto'})

    # Draw each graph in its subplot
    for i, ax in enumerate(axes):
        G = graphs[i]
        nodes,edges=PlotNetWork(G,500,5,ax,f'Time {i + 1}')

    # Layout adjustments for subplot to prevent overlap
    plt.subplots_adjust(bottom=0.2, top=0.8, left=0.05, right=0.95)

    # Add color bar for nodes
    cbar_pos = fig.add_axes([0.05, 0.1, 0.9, 0.03])  # x, y, width, height
    cbar = fig.colorbar(nodes, cax=cbar_pos, orientation='horizontal')
    cbar.set_label('Node Weight')

    # Add color bar for edges
    cbar_pos_edges = fig.add_axes([0.05, 0.05, 0.9, 0.03])  # x, y, width, height
    cbar_edges = fig.colorbar(edges, cax=cbar_pos_edges, orientation='horizontal')
    cbar_edges.set_label('Edge Weight')

    # Show the plot
    plt.show()
