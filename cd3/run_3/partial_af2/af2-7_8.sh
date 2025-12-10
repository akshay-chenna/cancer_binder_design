source /paperspace/apps/source_conda.sh
conda activate af2_binder_design_1

export CUDA_VISIBLE_DEVICES=3

export XLA_PYTHON_CLIENT_MEM_FRACTION=.48

for i in {1..8}
do
	/paperspace/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_top_7-${i}.silent -outsilent af2_top_7-${i}.silent -scorefilename af2_top_7-${i}.sc -checkpoint_name af2_top_7-${i}.point 
done &

for i in {1..8}
do
	/paperspace/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_top_8-${i}.silent -outsilent af2_top_8-${i}.silent -scorefilename af2_top_8-${i}.sc -checkpoint_name af2_top_8-${i}.point 
done &
wait
