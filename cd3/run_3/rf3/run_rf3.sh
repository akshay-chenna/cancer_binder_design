source ~/apps/source_conda.sh
conda deactivate
conda activate rf3

for i in *.pdb
do
	CUDA_VISIBLE_DEVICES=0 time rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/home/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2'
done
