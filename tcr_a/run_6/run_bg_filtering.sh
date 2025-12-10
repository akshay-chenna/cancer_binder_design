source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

awk -F , '{print $4 }' out_1/final_ranked_designs/final_designs_metrics_1000.csv | sort -u | nl -w 4 -n rz | sed 's/^/trca_bg_mb_ac_/' > trca_bg_mb_ac.csv
awk -F , '{print $4 }' out_2/final_ranked_designs/final_designs_metrics_500.csv | sort -u | nl -w 4 -n rz | sed 's/^/trca_bg_mb_nocys1_ac_/' >> trca_bg_mb_ac.csv
awk -F , '{print $4 }' out_3/final_ranked_designs/final_designs_metrics_500.csv | sort -u | nl -w 4 -n rz | sed 's/^/trca_bg_mb_nocys2_ac_/' >> trca_bg_mb_ac.csv

wait

