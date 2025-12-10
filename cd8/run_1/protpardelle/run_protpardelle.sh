source /paperspace/apps/source_conda.sh
conda deactivate
conda activate protpardelle

CUDA_VISIBLE_DEVICES=0 python3 -m protpardelle.sample run_1.yaml --num-samples 25 --motif-dir ./ --num-mpnn-seqs 0
