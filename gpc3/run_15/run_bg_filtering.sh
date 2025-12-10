source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

FILE=gpc3_nb2.yaml
DATE=2511020059

for i in {0..24} 
do
a=`wc -l ./out/boltzgen/${DATE}/unfiltered/job_${i}/intermediate_designs_inverse_folded/aggregate_metrics_analyze.csv | awk '{print $1}'`
b=$((a/10))
boltzgen run ./${FILE}  --output ./out/boltzgen/${DATE}/unfiltered/job_${i} --protocol protein-anything --budget $b --steps filtering

cat ./out/boltzgen/${DATE}/unfiltered/job_${i}/final_ranked_designs/final_designs_metrics_*.csv >> top_designs.csv

done

sort -u top_designs.csv > top_designs_10percent.csv
rm top_designs.csv
awk -F , '{print $4 }' top_designs_10percent.csv | sort -u | nl -w 4 -n rz | sed 's/^/gpc3_bg_nb_h2_ac_/' >> ../gpc3_bg_nb_ac.csv
awk -F , '{print $4 }' top_designs_10percent.csv | sort -u | nl -w 4 -n rz | sed 's/^/gpc3_bg_nb_h2_ac_/' >> orders_instance_gpc3_run_15_dec2025.csv
