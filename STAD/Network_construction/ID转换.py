# 安装和导入所需的库
import pandas as pd
from biomart import BiomartServer

# 连接到 BioMart 服务器
server = BiomartServer("http://www.ensembl.org/biomart")

# 选择数据库
database = server['ENSEMBL_MART_ENSEMBL']

# 选择数据集
dataset = database['hsapiens_gene_ensembl']

# 读取CSV文件
csv_data = pd.read_csv("./network.csv", header=None)

# 获取ensembl_ids
ensembl_ids = csv_data.iloc[0, 1:].values

# 查询基因名
response = dataset.query(attributes=["ensembl_gene_id", "hgnc_symbol"],
                         filters={"ensembl_gene_id": ensembl_ids})

# 将结果转换为DataFrame
result_df = pd.DataFrame(response)

# 将结果保存为CSV文件
result_df.to_csv("./gene_names.csv", index=False)
