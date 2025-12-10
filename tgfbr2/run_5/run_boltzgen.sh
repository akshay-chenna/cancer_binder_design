source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
mkdir out_1
CUDA_VISIBLE_DEVICES=3 boltzgen run ./tgfbr2_nb.yaml  --output ./out_1 --protocol protein-anything --num_designs 10000 --budget 1000
