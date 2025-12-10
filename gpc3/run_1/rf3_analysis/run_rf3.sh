<< 'END'
## without template
while read -r l
do
	CUDA_VISIBLE_DEVICES=0 rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2'
done < xaa &

while read -r l
do
	CUDA_VISIBLE_DEVICES=1 rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2'
done < xab &

while read -r l
do
	CUDA_VISIBLE_DEVICES=2 rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2'
done < xac &

while read -r l
do
	CUDA_VISIBLE_DEVICES=3 rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2'
done < xad &
wait

END

## with template
while read -r l
do
	CUDA_VISIBLE_DEVICES=0 rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2' template_selection="[A]"
done < xaa &

while read -r l
do
	CUDA_VISIBLE_DEVICES=1 rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2' template_selection="[A]"
done < xab &

while read -r l
do
	CUDA_VISIBLE_DEVICES=2 rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2' template_selection="[A]"
done < xac &

while read -r l
do
	CUDA_VISIBLE_DEVICES=3 rf3 fold inference_engine=rf3 inputs=$i ckpt_path='/paperspace/apps/modelforge/rf3_latest.pt' num_steps='200' early_stopping_plddt_threshold='0.2' template_selection="[A]"
done < xad &
wait
