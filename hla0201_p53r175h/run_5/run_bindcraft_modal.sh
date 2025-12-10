VOLUME_NAME=a0201_p53r175h_run_5 GPU=H100 TIMEOUT_HOURS=24 uvx --python 3.12 --with pandas modal run modal_parallel_bindcraft_volume.py \
    --settings hla0201_p53r175h_h2.json \
    --filters relaxed_filters.json \
    --advanced hotspot10_4stage_multimer_hardtarget.json \
    --pdb-in 6W51_AC_cleaned.pdb \
    --num-gpus 25

uvx modal volume get $VOLUME_NAME / volume_data --force
cd volume_data
mkdir Accepted
cp job_*/Accepted/*pdb Accepted/.
cat job_*/final_design_stats.csv | sort -u | tac > final_design_stats.csv
for i in Accepted/*pdb ; do python find_interface.py $i >> hotspot_overlap.txt ; done
