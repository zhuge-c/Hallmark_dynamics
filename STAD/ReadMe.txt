To complish the project, we should do many work, Here is a short deciption.

1. Get Data
we should get data, including normal tissues network, normal tissues expression, cancer network, from GRAND(https://grand.networkmedicine.org/cancers/CESC_cancer/#cardtissuetcga).

GRAND -> Network Construction
normal tissues adj/network -> normal_network
normal tissues expression -> normal_expression
cancer adj/network -> cancer_network

2. Convert Ensembl ID to Gene ID

we need use "ensembl id转换.R" program to convert ensembl ID to gene ID to get "gene_namme.csv".

Input is the "normal_network.csv"
Output is the "gene_name.csv"

3. Construction Network

The Input Data of "Network construction" are:
normal_network
normal_expression
cancer_network
gene_name

The output are:
cancer_Positive
normal_Positive
HM_normal_expression

4. Evolution Network

The Input Data of "Network Evolution" is from the Output of "Network construction"