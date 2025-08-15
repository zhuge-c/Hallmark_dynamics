import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import friedmanchisquare

font_size = 10

label_pad =5
def overlap_plot(file_path, result_path):
    # --- 1. 加载和准备数据 ---

    try:
        # 将第一列'Var1'作为索引（癌症名称）
        df = pd.read_csv(file_path, index_col='Var1')
    except FileNotFoundError:
        print(f"错误：找不到文件 '{file_path}'。请检查文件名和路径是否正确。")
        exit()

    # 将列名修改为更简洁的 H1, H2...
    # 假设 Cancer_overleap_1 对应 H1, 以此类推
    hallmark_codes = [f'H{i}' for i in range(1, 11)]
    df.columns = hallmark_codes

    # --- 2. 统计分析：弗里德曼检验 ---
    # a) 将原始数据转换为排序数据 (1-10)
    # axis=1 表示我们对每一行（每一种癌症）独立进行排序
    # ascending=False 表示值越大，排名越靠前（排名数字越小）
    rank_matrix = df.rank(axis=1, method='min', ascending=False)

    print("\n--- 排序矩阵 (前5行) ---")
    print(rank_matrix.head())

    # b) 执行弗里德曼检验
    # friedmanchisquare 需要接收每一组（Hallmark）的数据作为独立的参数
    # 所以我们将数据框的每一列传入
    stat, p_value = friedmanchisquare(*[rank_matrix[col] for col in rank_matrix.columns])

    print("\n--- 弗里德曼检验结果 ---")
    print(f"Chi-squared statistic: {stat:.4f}")
    print(f"P-value: {p_value:.2e}")

    friedman_results_df = pd.DataFrame({
        'Statistic': ['Chi-squared', 'p-value'],
        'Value': [stat, p_value]
    })

    if p_value < 0.05:
        print("\n结果解释: P值显著，表明在15种癌症中，Hallmark的排序存在一致的模式。")
    else:
        print("\n结果解释: P值不显著，没有足够的证据证明Hallmark的排序存在一致的模式。")


    # --- 3. 可视化：绘制平均排名条形图 ---
    # a) 计算每个Hallmark的平均排名
    average_ranks = rank_matrix.mean(axis=0)

    # b) 为了更清晰地展示，按平均排名对Hallmark进行排序
    sorted_ranks = average_ranks.sort_values(ascending=True)

    excel_filename = result_path + 'Overleap_Hallmark_Analysis_Results.xlsx'
    try:
        with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Original_Data')

            rank_matrix.to_excel(writer, sheet_name='Rank_Matrix')

            sorted_ranks.to_frame().to_excel(writer, sheet_name='Average_Ranks_Sorted')

            friedman_results_df.to_excel(writer, sheet_name='Friedman_Test_Summary', index=False)

        print(f"\n所有分析結果已成功匯出到 '{excel_filename}'")

    except Exception as e:
        print(f"\n匯出到 Excel 時發生錯誤: {e}")
        print("請確保您已安裝 'openpyxl' 函式庫 (pip install openpyxl)")

    # c) 绘制条形图
    width_cm = 10
    height_cm = 6
    width_inch = width_cm / 2.54
    height_inch = height_cm / 2.54

    plt.figure(figsize=(width_inch, height_inch))

    sns.set_theme(style="whitegrid", font_scale=1.2)

    ax = sorted_ranks.plot(kind='bar', color=sns.color_palette('viridis', 10))

    # 在每个条形图上方添加数值标签
    for p in ax.patches:
        ax.annotate(format(p.get_height(), '.2f'),
                    (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center',
                    xytext=(0, 9),
                    textcoords='offset points',
                    fontsize=font_size - 2, color='black')

    # d) 美化图表
    # ax.set_title('Average Rank of Hallmark Level Divergence Across 15 Cancer Types', fontsize=16, pad=20)
    ax.set_ylabel('Average Rank', fontsize=font_size, labelpad=label_pad)
    ax.set_xlabel('Hallmarks of Cancer', fontsize=font_size, labelpad=label_pad)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    ax.grid(axis='x') # 只保留Y轴的网格线

    ax.spines["top"].set_visible(False)

    # 在图上标注弗里德曼检验的p值
    text_str = f'p < {p_value:.2e}' if p_value > 0 else 'p < 2.2e-16'
    plt.text(0.3, 0.90, text_str, transform=ax.transAxes, fontsize=font_size,
             verticalalignment='top', horizontalalignment='right')

    plt.tight_layout()

    # e) 保存图表
    barchart_filename = result_path + 'Hallmark_Average_Rank_BarChart.png'
    plt.savefig(barchart_filename, dpi=1200)
    print(f"\n条形图已保存为 '{barchart_filename}'")
    plt.show()

def time_diff_plot(file_path, result_path):

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

    friedman_results_df = pd.DataFrame({
        'Statistic': ['Chi-squared', 'p-value'],
        'Value': [stat, p_value]
    })

    if p_value < 0.05:
        print("\n结果解释: P值显著，表明在15种癌症中，Hallmark的t2-t1时间差排序存在一致的模式。")
    else:
        print("\n结果解释: P值不显著，没有足够的证据证明排序存在一致的模式。")

    # --- 3. 可视化：绘制平均排名条形图 ---
    # a) 计算每个Hallmark的平均排名
    average_ranks = rank_matrix.mean(axis=0)

    # b) 为了更清晰地展示，按平均排名对Hallmark进行排序
    sorted_ranks = average_ranks.sort_values(ascending=True)

    excel_filename = result_path + 'Time_Hallmark_Analysis_Results.xlsx'
    try:
        with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
            df.to_excel(writer, sheet_name='Original_Data')

            rank_matrix.to_excel(writer, sheet_name='Rank_Matrix')

            sorted_ranks.to_frame().to_excel(writer, sheet_name='Average_Ranks_Sorted')

            friedman_results_df.to_excel(writer, sheet_name='Friedman_Test_Summary', index=False)

        print(f"\n所有分析結果已成功匯出到 '{excel_filename}'")

    except Exception as e:
        print(f"\n匯出到 Excel 時發生錯誤: {e}")
        print("請確保您已安裝 'openpyxl' 函式庫 (pip install openpyxl)")

    # c) 绘制条形图
    width_cm = 10
    height_cm = 6
    width_inch = width_cm / 2.54
    height_inch = height_cm / 2.54

    plt.figure(figsize=(width_inch, height_inch))
    sns.set_theme(style="whitegrid", font_scale=1.2)
    ax = sorted_ranks.plot(kind='bar', color=sns.color_palette('viridis', 10))

    # 在每个条形图上方添加数值标签
    for p in ax.patches:
        ax.annotate(format(p.get_height(), '.2f'),
                    (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='center',
                    xytext=(0, 9),
                    textcoords='offset points',
                    fontsize=font_size - 2, color='black')

    # d) 美化图表
    # ax.set_title(r'Average Rank of Time Difference ($\Delta t = t_2 - t_1$) Across 15 Cancer Types', fontsize=16, pad=20)
    ax.set_ylabel('Average Rank', fontsize=font_size, labelpad=label_pad)
    ax.set_xlabel('Hallmarks of Cancer', fontsize=font_size, labelpad=label_pad)
    plt.xticks(fontsize=font_size)
    plt.yticks(fontsize=font_size)
    ax.grid(axis='x')  # 只保留Y轴的网格线

    ax.spines["top"].set_visible(False)
    # 在图上标注弗里德曼检验的p值
    text_str = f'p < {p_value:.2e}' if p_value > 0 else 'p < 2.2e-16'
    plt.text(0.3, 0.90, text_str, transform=ax.transAxes, fontsize=font_size,
             verticalalignment='top', horizontalalignment='right')
    # , bbox = dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5)

    plt.tight_layout()

    # e) 保存图表
    barchart_filename = result_path + 'Time_Difference_Rank_BarChart.png'
    plt.savefig(barchart_filename, dpi=1200)
    print(f"\n条形图已保存为 '{barchart_filename}'")
    plt.show()

