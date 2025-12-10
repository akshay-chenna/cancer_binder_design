source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
#mkdir out_{1..2}
#CUDA_VISIBLE_DEVICES=0 boltzgen run ./cd8_mb.yaml  --output ./out_1 --protocol protein-anything --num_designs 15000 --budget 1500 &
#CUDA_VISIBLE_DEVICES=1 boltzgen run ./cd8_mb.yaml  --output ./out_2 --protocol protein-anything --inverse_fold_avoid 'C' --num_designs 15000 --budget 1500 &

CUDA_VISIBLE_DEVICES=0 boltzgen run ./cd8_mb.yaml  --output ./out_1 --protocol protein-anything --steps design_folding analysis filtering &
CUDA_VISIBLE_DEVICES=1 boltzgen run ./cd8_mb.yaml  --output ./out_2 --protocol protein-anything --steps design_folding analysis filtering &

wait
