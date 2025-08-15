function test_output(p, tbl, stats, outfile)

% 打开写入文件
fid = fopen(outfile, 'w');

% 写入 p 值
fprintf(fid, '=== Friedman Test Result ===\n\n');
fprintf(fid, 'p-value: %.6f\n\n', p);

% 写入表格 tbl（逐行输出 cell 内容）
fprintf(fid, '--- ANOVA Table ---\n');
for i = 1:size(tbl, 1)
    for j = 1:size(tbl, 2)
        entry = tbl{i,j};
        if isnumeric(entry)
            fprintf(fid, '%-15.4f', entry);
        elseif ischar(entry)
            fprintf(fid, '%-15s', entry);
        else
            fprintf(fid, '%-15s', '');  % 空白处理
        end
    end
    fprintf(fid, '\n');
end
fprintf(fid, '\n');

% 写入 stats 结构体
fprintf(fid, '--- Stats Structure ---\n');
stats_fields = fieldnames(stats);
for i = 1:length(stats_fields)
    key = stats_fields{i};
    val = stats.(key);
    if isnumeric(val)
        fprintf(fid, '%s: %.6f\n', key, val);
    else
        fprintf(fid, '%s: [non-numeric data]\n', key);
    end
end

% 关闭文件
fclose(fid);

disp("file output")