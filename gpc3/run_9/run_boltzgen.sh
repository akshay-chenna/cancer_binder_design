source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
mkdir out_{1..2}
CUDA_VISIBLE_DEVICES=4 boltzgen run ./gpc3_nb.yaml  --output ./out_1 --protocol protein-anything --num_designs 15000 --budget 750 &
CUDA_VISIBLE_DEVICES=5 boltzgen run ./gpc3_nb.yaml  --output ./out_2 --protocol protein-anything --num_designs 15000 --budget 750 &
wait
