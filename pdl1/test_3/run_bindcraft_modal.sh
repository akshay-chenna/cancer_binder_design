VOLUME_NAME=actest GPU=H100 TIMEOUT_HOURS=24 uvx --python 3.12 --with pandas modal run modal_parallel_bindcraft_volume.py \
    --settings PDL1.json \
    --filters relaxed_filters.json \
    --advanced hotspot10_4stage_multimer_hardtarget.json \
    --pdb-in PDL1.pdb \
    --num-gpus 3
