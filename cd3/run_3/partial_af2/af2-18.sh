source /paperspace/apps/source_conda.sh
conda activate af2_binder_design_1

export XLA_PYTHON_CLIENT_MEM_FRACTION=.48

for i in {1..2}
do
	CUDA_VISIBLE_DEVICES=1 /paperspace/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_top_18-${i}.silent -outsilent af2_top_18-${i}.silent -scorefilename af2_top_18-${i}.sc -checkpoint_name af2_top_18-${i}.point 
done &

for i in {3..4}
do
	CUDA_VISIBLE_DEVICES=1 /paperspace/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_top_18-${i}.silent -outsilent af2_top_18-${i}.silent -scorefilename af2_top_18-${i}.sc -checkpoint_name af2_top_18-${i}.point 
done &

for i in {5..6}
do
	CUDA_VISIBLE_DEVICES=3 /paperspace/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_top_18-${i}.silent -outsilent af2_top_18-${i}.silent -scorefilename af2_top_18-${i}.sc -checkpoint_name af2_top_18-${i}.point 
done &

for i in {7..8}
do
	CUDA_VISIBLE_DEVICES=3 /paperspace/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_top_18-${i}.silent -outsilent af2_top_18-${i}.silent -scorefilename af2_top_18-${i}.sc -checkpoint_name af2_top_18-${i}.point 
done &
wait
