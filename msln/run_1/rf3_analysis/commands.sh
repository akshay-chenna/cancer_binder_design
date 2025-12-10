# merge all accepted score files
cat outputs/*csv | awk '!seen[$0]++' >> bindcraft_passed-rf3_metrics.csv 

cat ../out_hotspots_*/final_design_stats.csv | awk '!seen[$0]++' >> bindcraft_passed-bindcraft_metrics.csv
