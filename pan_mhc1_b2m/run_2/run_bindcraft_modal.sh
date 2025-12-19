VOLUME_NAME=panmhc1b2m_run2 GPU=H100 TIMEOUT_HOURS=24 uvx --python 3.12 --with pandas modal run modal_parallel_bindcraft_volume.py \
    --settings panmhc1_b2m.json \
    --filters default_filters.json \
    --advanced default_4stage_multimer.json \
    --pdb-in 7uc5_cleaned_renum.pdb \
    --num-gpus 2

#uvx modal volume get $VOLUME_NAME / volume_data --force
#cd volume_data
#mkdir Accepted
#cp job_*/Accepted/*pdb Accepted/.
#cat job_*/final_design_stats.csv | sort -u | tac > final_design_stats.csv
#for i in Accepted/*pdb ; do python find_interface.py $i >> hotspot_overlap.txt ; done
