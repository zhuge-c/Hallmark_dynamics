import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import friedmanchisquare

# --- 1. 加载和准备数据 ---
file_path = './cancer_time_t1t2.csv'
try:
    # 将第一列'Var1'作为索引（癌症名称）
    df = pd.read_csv(file_path, index_col='Var1')
except FileNotFoundError:
    print(f"错误：找不到文件 '{file_path}'。请检查文件名和路径是否正确。")
    exit()

# 将列名修改为更简洁的 H1, H2...
# 假设 Cancer_Time_diference_1 对应 H1, 以此类推
hallmark_codes = [f'H{i}' for i in range(1, 11)]
df.columns = hallmark_codes

# --- 2. 统计分析：弗里德曼检验 ---
# a) 将原始数据转换为排序数据 (1-10)
# axis=1 表示我们对每一行（每一种癌症）独立进行排序
# ascending=False 表示值越大（时间差越大），排名越靠后（排名数字越大）
# 注意：这里我们假设 t2-t1 的值越大，代表越“不重要”或越“晚发生”，因此排名靠后。
# 如果您的逻辑相反，请将 ascending 改为 True。
rank_matrix = df.rank(axis=1, method='min', ascending=False)

print("\n--- 排序矩阵 (前5行) ---")
print(rank_matrix.head())

# b) 执行弗里德曼检验
stat, p_value = friedmanchisquare(*[rank_matrix[col] for col in rank_matrix.columns])

print("\n--- 弗里德曼检验结果 ---")
print(f"Chi-squared statistic: {stat:.4f}")
print(f"P-value: {p_value:.2e}")

if p_value < 0.05:
    print("\n结果解释: P值显著，表明在15种癌症中，Hallmark的t2-t1时间差排序存在一致的模式。")
else:
    print("\n结果解释: P值不显著，没有足够的证据证明排序存在一致的模式。")


# --- 3. 可视化：绘制平均排名条形图 ---
# a) 计算每个Hallmark的平均排名
average_ranks = rank_matrix.mean(axis=0)

# b) 为了更清晰地展示，按平均排名对Hallmark进行排序
sorted_ranks = average_ranks.sort_values(ascending=True)

# c) 绘制条形图
plt.figure(figsize=(9.5, 5.4))
sns.set_theme(style="whitegrid", font_scale=1.2)
ax = sorted_ranks.plot(kind='bar', color=sns.color_palette('viridis', 10))

# d) 美化图表
ax.set_title(r'Average Rank of Time Difference ($\Delta t = t_2 - t_1$) Across 15 Cancer Types', fontsize=16, pad=20)
ax.set_ylabel('Average Rank (Lower is Earlier)', fontsize=16, labelpad=15)
ax.set_xlabel('Hallmarks of Cancer', fontsize=16, labelpad=15)
plt.xticks(rotation=45, ha="right", fontsize=14)
plt.yticks(fontsize=14)
ax.grid(axis='x') # 只保留Y轴的网格线

# 在图上标注弗里德曼检验的p值
text_str = f'Friedman Test\np < {p_value:.2e}' if p_value > 0 else 'Friedman Test\np < 2.2e-16'
plt.text(0.2, 0.90, text_str, transform=ax.transAxes, fontsize=14,
         verticalalignment='top', horizontalalignment='right',
         bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

plt.tight_layout()

# e) 保存图表
barchart_filename = 'Time_Difference_Rank_BarChart.png'
plt.savefig(barchart_filename, dpi=300)
print(f"\n条形图已保存为 '{barchart_filename}'")
plt.show()

