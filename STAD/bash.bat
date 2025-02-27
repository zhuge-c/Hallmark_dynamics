set start_time=%TIME%
echo The start time: %start_time%

cd  ./Network_construction
python Main.py

copy /Y .\output\cancer_Positive.csv ..\Network_evolution\orginal_data
copy /Y .\output\normal_Positive.csv ..\Network_evolution\orginal_data
copy /Y .\output\HM_normal_expression.csv ..\Network_evolution\orginal_data

cd ..\Network_evolution
matlab -r "run('main.m'); run('main_analysis.m')"

set end_time=%TIME%
echo The end time: %end_time%

cmd
