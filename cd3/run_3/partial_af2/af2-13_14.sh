source /paperspace/apps/source_conda.sh
conda activate af2_binder_design_1

export CUDA_VISIBLE_DEVICES=2

export XLA_PYTHON_CLIENT_MEM_FRACTION=.48

for i in {1..8}
do
	/paperspace/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_top_13-${i}.silent -outsilent af2_top_13-${i}.silent -scorefilename af2_top_13-${i}.sc -checkpoint_name af2_top_13-${i}.point 
done &

for i in {1..8}
do
	/paperspace/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_top_14-${i}.silent -outsilent af2_top_14-${i}.silent -scorefilename af2_top_14-${i}.sc -checkpoint_name af2_top_14-${i}.point 
done &
wait
