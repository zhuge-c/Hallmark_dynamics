import pandas as pd
import pickle
import os

from ID_Symbol_Correct import ID_Symbol_Correct
from Get_Hallmarks_gene import Get_Hallmarks_gene
from Construct_Graph import Construct_Graph
from OutPut_Result import Print_Graph_Adj
from OutPut_Result import Print_Negative_Positive_Adj
from OutPut_Result import Plot_NetWork
from OutPut_Result import GraphData_Process
from HallMark_Expression import *
import time

start_time=time.time()
flag=0
type_name=["normal","cancer"]
# ########################################################
Hallmark_gene_file="HG.pkl"
for tn in type_name:
    input_network_data=f"{tn}_network.csv"
    output_network_data=f"new_{tn}_network.csv"
    output_pkl_data=f"{tn}_data.pkl"
    output_graph_data=f"{tn}_Graph.pkl"
    output_network_fig=f"{tn}_network.jpg"
    output_negative_adj=f"{tn}_Negative.csv"
    output_positive_adj=f"{tn}_Positive.csv"
    ##################################
    # 输入基因-iD 对应文件，网络数据文件，
    # 输出新命名网络，需要太多时间
    ID_Symbol_Correct("./input/gene_names.csv","./input/"+input_network_data,"./output/"+output_network_data)
    ##################################
    # 输入Hallmarks 基因所在文件夹，对应网络数据
    # 输出pkl数据
    #  "Hallmarks_TFS":Hallmarks_TFS,
    #  "Hallmarks_Gene":Hallmarks_Gene,
    #  "network":network
    Get_Hallmarks_gene("./HallmarkGene/", "./output/"+output_network_data, "./output/"+Hallmark_gene_file, "./output/"+output_pkl_data)
    ##################################
    # 输入 Hallmarks pkl 数据
    # 输出 pkl 图数据
    Construct_Graph("./output/"+output_pkl_data,"./output/"+output_graph_data)
    ##################################
    # 如果图数据权重不一样需要进行处理
    if flag==1:
        GraphData_Process("./output/"+output_graph_data,5000.)
    ##################################
    # 输入网络图数据，绘制网络图
    Plot_NetWork("./output/"+output_graph_data,"./output/"+output_network_fig)
    ##################################
    # 输出正负邻接矩阵
    Print_Negative_Positive_Adj("./output/"+output_graph_data,"./output/"+output_negative_adj,"./output/"+output_positive_adj)

########################################################

######################################################
# output Hallmark expression
###############################
ty_nm="normal"
expression_file_path = f"./input/{ty_nm}_expression.csv"
gene_name_file = "./input/gene_names.csv"
gene_expression_file = f"./output/Gene_{ty_nm}_expression_data.csv"
Hallmark_data_file = f"./output/{ty_nm}_data.pkl"
Hallmark_Expression_file = f"./output/HM_{ty_nm}_expression.csv"

Hallmark_Expression(expression_file_path, gene_name_file, gene_expression_file, Hallmark_data_file, Hallmark_Expression_file)

######################### 计算运行时间
end_time=time.time()
execution_time = end_time - start_time

print(f"代码执行时间为: {execution_time} 秒")
print(f"代码执行时间为: {execution_time/60} 分钟")
