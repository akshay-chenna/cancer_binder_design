source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
mkdir out_{1..1}
CUDA_VISIBLE_DEVICES=6 boltzgen run ./gpc3_mb2.yaml  --output ./out_1 --protocol protein-anything --num_designs 15000 --budget 1500 &
#CUDA_VISIBLE_DEVICES=3 boltzgen run ./gpc3_mb.yaml  --output ./out_2 --protocol protein-anything --inverse_fold_avoid 'C' --num_designs 15000 --budget 1500 &
wait
