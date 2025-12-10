source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

FILE=cd3e_mb.yaml
DATE=2511021909

for i in {0..24} 
do
a=`wc -l ./out/boltzgen/${DATE}/unfiltered/job_${i}/intermediate_designs_inverse_folded/aggregate_metrics_analyze.csv | awk '{print $1}'`
b=$((a/10))
boltzgen run ./${FILE}  --output ./out/boltzgen/${DATE}/unfiltered/job_${i} --protocol protein-anything --budget $b --steps filtering

cat ./out/boltzgen/${DATE}/unfiltered/job_${i}/final_ranked_designs/final_designs_metrics_*.csv >> top_designs.csv

done

sort -u top_designs.csv > top_designs_10percent.csv
rm top_designs.csv
awk -F , '{print $NF }' top_designs_10percent.csv | sort -u | nl -w 4 -n rz | sed 's/^/cd3e_bg_mb_ac_/' > cd3e_bg_mb_ac.csv
