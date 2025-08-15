library(R.matlab) # 用于读取 .mat 文件
library(tidyverse) # 用于数据整理和绘图 (包含了 ggplot2, dplyr, tidyr)
library(ggpubr)    # 用于添加显著性标记

setwd("I:/ProjectResult/STAD")

mat_data <- readMat("BoxData.mat")