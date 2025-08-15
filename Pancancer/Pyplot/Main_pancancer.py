from fig_plot import overlap_plot, time_diff_plot


file_path = './input/cancer_overleap.csv'
result_path = "./output/"

overlap_plot(file_path, result_path)

file_path2 = './input/cancer_time_t1t2.csv'

time_diff_plot(file_path2, result_path)