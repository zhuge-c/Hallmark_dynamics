ensembl id转换.R

该程序用来将基因ID转换为基因Symbol，输入数据为网络邻接矩阵数据，输出为基因ID-Symbol数据

ID_Symbol_Correct.py

该程序用来将基因symbol补充完全（因为ID转换symbol后，会有新的基因缺少，所以需要以基因ID补全）。然后将补全的基因ID应用到network数据中，将列名替换，重新输出network数据。

NetworkConstruct.py

该程序用来统计各个Hallmarks所拥有的基因和转录因子，最后将获得的结果通过.pkl数据存储

NetworkConstruct2.py

