GPU=H100 TIMEOUT_HOURS=1 uvx --python 3.12 --with pandas modal run modal_parallel_bindcraft.py \
    --settings PDL1.json \
    --filters no_filters.json \
    --advanced default_4stage_multimer.json \
    --pdb-in PDL1.pdb \
    --num-gpus 3 
