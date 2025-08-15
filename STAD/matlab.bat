cd  ./Hallmark_Network

copy /Y .\output\cancer_Positive.csv ..\Network_evolution\orginal_data
copy /Y .\output\normal_Positive.csv ..\Network_evolution\orginal_data
copy /Y .\output\HM_normal_expression.csv ..\Network_evolution\orginal_data


cd ..\Network_evolution 

matlab -r "run('main.m'); run('Result.m')" || goto error

echo Successfully accomplish
goto end


:error
echo Something is wrong
goto end

: end

set end_time=%TIME%
echo The end time: %end_time%

cmd
