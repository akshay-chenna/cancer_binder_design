source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

FILE=hla_mb.yaml
BUDGET=5000
DATE=2511121535


cp ./out/boltzgen/${DATE}/combined/intermediate_designs_inverse_folded/aggregate_metrics_analyze.csv ./out/boltzgen/${DATE}/combined/intermediate_designs_inverse_folded/old.csv
awk 'BEGIN{FS=OFS=","} NR > 1 {$2 = $2 ".cif"} 1' ./out/boltzgen/${DATE}/combined/intermediate_designs_inverse_folded/old.csv > ./out/boltzgen/${DATE}/combined/intermediate_designs_inverse_folded/aggregate_metrics_analyze.csv
boltzgen run ./${FILE}  --output ./out/boltzgen/${DATE}/combined --protocol protein-anything --budget ${BUDGET} --steps filtering

awk -F , '{print $4}' out/boltzgen/${DATE}/combined/final_ranked_designs/final_designs_metrics_${BUDGET}.csv | sort -u | nl -w 4 -n rz | head -${BUDGET} | sed 's/^/a0201tp53r175h_bg_nb_ac_/'
