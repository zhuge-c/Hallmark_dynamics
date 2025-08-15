set start_time=%TIME%
echo The start time: %start_time%

cd  ./Hallmark_Network

:: cd ./input
:: Rscript ensembl_id_convert.R || goto error
:: cd ../

python Main.py  || goto error


goto end


:error
echo Something is wrong
goto end

: end

set end_time=%TIME%
echo The end time: %end_time%

cmd
