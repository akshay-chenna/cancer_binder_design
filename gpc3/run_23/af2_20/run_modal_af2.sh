#uvx modal volume delete af2_akshay -y
#uvx modal volume create --version=2 af2_akshay
#uvx modal volume put af2_akshay reordered_modified_gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb

#VOLUME_NAME=af2_akshay GPU=L40S TIMEOUT_MINUTES=1440 uvx --with numpy modal run modal_binder_af2.py --path ./ --out-path ./outputs --binder-chain "A" --target-chain "B" --num-gpus 1
VOLUME_NAME=af2_akshay GPU=L40S TIMEOUT_MINUTES=1440 uvx --with numpy modal run modal_binder_af2.py --path ./ --out-path ./outputs_2 --binder-chain "A" --target-chain "B" --num-gpus 1 --no-use-multimer
#rm -rf volume_data
#mkdir volume_data
uvx modal volume get af2_akshay . volume_data --force
