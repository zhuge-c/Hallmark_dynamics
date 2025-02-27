import pandas as pd

def ID_Symbol_Correct(gene_name_file,network_file,new_network_file):
    gene_name= pd.read_csv(gene_name_file)
    gene_name.columns=["Gene_id","Gene_symbol"]
    for r in gene_name.index:
        temp_row=gene_name.loc[r]
        if str(temp_row["Gene_symbol"])=="nan":
            gene_name.loc[r]["Gene_symbol"]=temp_row["Gene_id"]

    gene_name.to_csv("./input/new_gene_names.csv")

    network=pd.read_csv(network_file)
    for i in range(1,len(network.columns)):
        if gene_name["Gene_id"].isin([str(network.columns[i])]).any():
            nr=gene_name[gene_name["Gene_id"].isin([str(network.columns[i])])].index[0]
            network.rename(columns={network.columns[i]:gene_name.loc[nr]["Gene_symbol"]},inplace=True)

    gene_tfs=network[network.columns[0]]

    network=network.drop(network.columns[0],axis=1)

    network=network.rename(index=gene_tfs)

    network.to_csv(new_network_file)
