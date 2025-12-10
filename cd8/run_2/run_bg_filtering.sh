source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

boltzgen run ./cd8_mb.yaml  --output ./out_1 --protocol protein-anything --steps filtering --budget 1500
boltzgen run ./cd8_mb.yaml  --output ./out_2 --protocol protein-anything --steps filtering --budget 1500

awk -F , '{print $4 }' out_1/final_ranked_designs/final_designs_metrics_1500.csv | sort -u | nl -w 4 -n rz | sed 's/^/cd8_bg_mb_ac_/' >> cd8_bg_mb_ac.csv
awk -F , '{print $4 }' out_2/final_ranked_designs/final_designs_metrics_1500.csv | sort -u | nl -w 4 -n rz | sed 's/^/cd8_bg_mb_nocys_ac_/' >> cd8_bg_mb_ac.csv

wait

