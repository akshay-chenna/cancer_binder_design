source ~/apps/source_conda.sh
conda deactivate
conda activate proteinhunter

#export CUDA_VISIBLE_DEVICES=0

python boltz_ph/design.py --num_designs 1000 --num_cycles 20 --protein_ids A:B \
--template_path 6W51_AC_cleaned.pdb --template_chain_id A:B \
--protein_msas "" --gpu_id 0 \
--name r175h_p53_0201 \
--percent_X 30 --min_design_protein_length 50 --max_design_protein_length 150 \
--high_iptm_threshold 0.7 \
--contact_residues 73,76,77,80,81,84,143,146,147,150,152,187,188,189 --no_potentials False \
--work_dir ~/cancer/hla0201_p53r175h/run_1 \
--plot