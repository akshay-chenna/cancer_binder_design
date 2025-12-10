source ~/apps/source_conda.sh
conda deactivate
conda activate FoldCraft

export XLA_PYTHON_CLIENT_MEM_FRACTION=.23
for i in {1..4}
do

CUDA_VISIBLE_DEVICES=0 python /home/paperspace/apps/FoldCraft/FoldCraft.py --output_folder 7dv4vhh_${i} --binder_template 7dv4_vhh_serialnumbered.pdb --target_template cd3e_cl_00.pdb --target_hotspots '33-38,50-56,77-81'  --binder_mask '26-33,51-57,96-108' --binder_chain 'H'  --num_designs 500 & 

done
wait
