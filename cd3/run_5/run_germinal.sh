source /paperspace/apps/source_conda.sh
conda deactivate
conda activate germinal

for i in {0..3}
do

CUDA_VISIBLE_DEVICES=$i python run_germinal.py run=vhh experiment_name="vhh_cd3e_${i}" target.target_hotspots=\"A4,A15,A18,A21,A22,A46,A50,A76\" target.target_name="cd3e" target.target_pdb_path="pdbs/cd3e_cl_00_renum.pdb" target.target_chain="A" target.binder_chain="B" target.length=89 &

done
wait
