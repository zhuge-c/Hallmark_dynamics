rm(list = ls())
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# library(BiocManager)

library(biomaRt) # 多种数据库

my_mart <- useMart("ensembl")
datasets <- listDatasets(my_mart)
my_dataset <- useDataset("hsapiens_gene_ensembl",
                  mart = my_mart)

library(readr)
csv_data <- read_csv("./normal_network.csv", col_names = FALSE)

ensembl_ids <- t(csv_data[1, -1])

# 查询基因名
gene_names <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = ensembl_ids,
                    mart = my_dataset)

write.csv(gene_names, "./gene_names.csv", row.names = FALSE)
