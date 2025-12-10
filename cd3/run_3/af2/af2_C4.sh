source ~/scripts/source_conda.sh
conda activate af2_binder_design

export CUDA_VISIBLE_DEVICES=1

export XLA_PYTHON_CLIENT_MEM_FRACTION=.11

for i in {1..8}
do
        ~/apps/dl_binder_design/af2_initial_guess/predict.py -silent mpnn_C4_${i}.silent -outsilent af2_C4_${i}.silent -scorefilename af2_C4_${i}.sc -checkpoint_name af2_C4_${i}.point &
done
wait
