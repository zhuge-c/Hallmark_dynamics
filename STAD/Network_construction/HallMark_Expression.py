import os
import pandas as pd
import pickle
import warnings
warnings.filterwarnings("ignore")

######################
## 读取表达数据
#####################
def Hallmark_Expression(normal_file_path, gene_name_file, normal_gene_expression_file, Hallmark_data_file, Hallmark_Expression_file):

    expression_data=pd.read_csv(normal_file_path,index_col=0)
    gene_name=pd.read_csv(gene_name_file,index_col=0)

    expression_data_id=expression_data.index
    gene_name_id=gene_name.index

    intersetion_id=set(expression_data_id).intersection(set(gene_name_id))

    Gene_sym_Data=pd.DataFrame()
    it=0
    for id in intersetion_id:
        it=it+1
        temp_symbol=list(gene_name.loc[id])[0]
        expression_data.rename(index={id:temp_symbol},inplace=True) # 改变行名
        # Gene_sym_Data=pd.concat([Gene_sym_Data,expression_data.loc[temp_symbol]],axis=0)

    expression_data.to_csv(normal_gene_expression_file)
    # 该数据行为基因，列表达不同的样本


    ##################################
    #### 统计 Hallmark 表达量
    ##################################


    HallMarkData=Hallmark_data_file
    # 从文件中读取多个变量
    with open(HallMarkData, "rb") as f:
        loaded_data = pickle.load(f)
        Hallmarks_TFS = loaded_data["Hallmarks_TFS"]
        Hallmarks_Gene = loaded_data["Hallmarks_Gene"]
        network = loaded_data["network"]

    HM_expression=pd.DataFrame()
    Gene_list=[]
    for key in Hallmarks_Gene.keys():
        temp_all_gene=set(Hallmarks_Gene[key]).union(set(Hallmarks_TFS[key]))
        expression_gene=list(set(expression_data.index).intersection(set(temp_all_gene)))
        expression_gene_symbol=expression_data.loc[expression_gene]
        hm=expression_gene_symbol.sum(axis=0)
        # print(hm.mean())
        HM_expression[key]=[hm.mean()]

    HM_expression.to_csv(Hallmark_Expression_file,index=False)

