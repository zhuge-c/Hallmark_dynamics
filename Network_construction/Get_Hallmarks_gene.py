import os
import re
import pandas as pd

def Get_Hallmarks_gene(Hallmarks_path,network_file,output_data):
    ################# 获取文件夹下所有txt
    folder_path = Hallmarks_path  # 替换为你的文件夹路径
    txt_files = []  # 用于存储所有的 .txt 文件路径
    for file in os.listdir(folder_path):
        if file.endswith(".txt"):
            txt_files.append(file)
    ####################### 获取所有Hallmarks
    Hallmarks=[]
    pattern=r"\((.*?)\)"
    for f in txt_files:
        temp=re.search(pattern,f)
        Hallmarks.append(temp.group(1))
    ######################### 读取基因数据
    Total_Hallmarks_Gene={}
    for file in txt_files:
        index=txt_files.index(file)
        hm=Hallmarks[index]
        data=pd.read_csv(Hallmarks_path+file,header=None)
        Total_Hallmarks_Gene[hm]=list(data[data.columns[0]])
    ########################### 读取网络数据
    network=pd.read_csv(network_file)
    gene_tfs=network[network.columns[0]]
    network=network.drop(network.columns[0],axis=1)
    network=network.rename(index=gene_tfs)

    Hallmarks_Gene={}
    # 确定Hallmarks中的基因
    for key in Hallmarks:
        tt=Total_Hallmarks_Gene[key]
        temp1=set(tt)
        temp2=set(network.columns)
        temp_total_gene=temp1.intersection(temp2)
        Hallmarks_Gene[key]=list(temp_total_gene)

    # 确定Hallmarks 中的TF
    Hallmarks_TFS={}
    for key in Hallmarks:
        tt=Total_Hallmarks_Gene[key]
        temp1=set(tt)
        temp2=set(gene_tfs)
        temp_tfs_gene=temp1.intersection(temp2)
        Hallmarks_TFS[key]=list(temp_tfs_gene)

    # 存储数据
    import pickle
    pickle_data={"Hallmarks_TFS":Hallmarks_TFS,
                 "Hallmarks_Gene":Hallmarks_Gene,
                 "network":network}

    with open(output_data, "wb") as f:
        pickle.dump(pickle_data, f)