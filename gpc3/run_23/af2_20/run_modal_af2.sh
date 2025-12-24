uvx modal volume put af2_akshay gpc3_run_23_af2_20 ./
VOLUME_NAME=af2_akshay GPU=A100 TIMEOUT_MINUTES=1440 uvx modal run modal_binder_af2.py --path ./gpc3_run_23_af2_20 --out-path ./gpc3_run_23_af2_20_outputs --binder-chain "A" --target-chain "B" --num-gpus 100
mkdir volume_data
uvx modal volume get af2_akshay ./gpc3_run_23_af2_20_outputs volume_data --force
