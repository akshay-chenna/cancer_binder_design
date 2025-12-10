source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
mkdir out_{1..3}
CUDA_VISIBLE_DEVICES=2 boltzgen run ./tcr_mb.yaml  --output ./out_1 --protocol protein-anything --num_designs 10000 --budget 1000 --reuse &
CUDA_VISIBLE_DEVICES=3 boltzgen run ./tcr_mb.yaml  --output ./out_2 --protocol protein-anything --inverse_fold_avoid 'C' --num_designs 5000 --budget 500 &
CUDA_VISIBLE_DEVICES=4 boltzgen run ./tcr_mb.yaml  --output ./out_3 --protocol protein-anything --inverse_fold_avoid 'C' --num_designs 5000 --budget 500 &
wait
