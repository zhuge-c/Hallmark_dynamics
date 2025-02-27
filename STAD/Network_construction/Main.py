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

########################################################
# 输入基因-iD 对应文件，网络数据文件，
# 输出新命名网络，需要太多时间

# ID_Symbol_Correct("./input/gene_names.csv","./input/cancer_network.csv","./output/new_cancer_network.csv")

########################################################
########################################################
# 输入Hallmarks 基因所在文件夹，对应网络数据
# 输出pkl数据
#  "Hallmarks_TFS":Hallmarks_TFS,
#  "Hallmarks_Gene":Hallmarks_Gene,
#  "network":network
########################################

# Get_Hallmarks_gene("./HallmarkGene/","./output/new_cancer_network.csv","./output/cancer_data.pkl")

# ########################################################
# ########################################################
# 输入 Hallmark

# Construct_Graph("./output/cancer_data.pkl","./output/cancer_Graph.pkl")

########################################################

# ########################################################
# # 如果图数据权重不一样需要进行处理
# #######  GraphData_Process("./output/Graph.pkl",50000.)
# Plot_NetWork("./output/cancer_Graph.pkl","./output/cancer_network.jpg")
#
# Print_Negative_Positive_Adj("./output/cancer_Graph.pkl","./output/cancer_Negative.csv","./output/cancer_Positive.csv")

# ########################################################
# # 正常基因
ID_Symbol_Correct("./input/gene_names.csv","./input/normal_network.csv","./output/new_normal_network.csv")

Get_Hallmarks_gene("./HallmarkGene/","./output/new_normal_network.csv","./output/normal_data.pkl")

Construct_Graph("./output/normal_data.pkl","./output/normal_Graph.pkl")

# GraphData_Process("./normal_Graph.pkl",5000.)

Plot_NetWork("./output/normal_Graph.pkl","./output/normal_network.jpg")

Print_Negative_Positive_Adj("./output/normal_Graph.pkl","./output/normal_Negative.csv","./output/normal_Positive.csv")
########################################################

######################################################
# output Hallmark expression
###############################
normal_file_path = "./input/normal_expression.csv"
gene_name_file = "./input/gene_names.csv"
normal_gene_expression_file = "./output/Gene_normal_expression_data.csv"
Hallmark_data_file = "./output/normal_data.pkl"
Hallmark_Expression_file = "./output/HM_normal_expression.csv"

Hallmark_Expression(normal_file_path, gene_name_file, normal_gene_expression_file, Hallmark_data_file, Hallmark_Expression_file)

######################### 计算运行时间
end_time=time.time()
execution_time = end_time - start_time

print(f"代码执行时间为: {execution_time} 秒")
print(f"代码执行时间为: {execution_time/60} 分钟")
