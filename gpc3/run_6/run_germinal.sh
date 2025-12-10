source /paperspace/apps/source_conda.sh
conda deactivate
conda activate germinal

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
for i in {2..3}
do

CUDA_VISIBLE_DEVICES=$i python run_germinal.py experiment_name="scfv_gpc3_hotspots_2_${i}" target.target_hotspots=\"A22,A73,A77,A42,A49\" target.target_name="gpc3_h2" target.target_pdb_path="pdbs/correct_uniprot_gpc3_hotspot2_renum.pdb" target.target_chain="A" target.binder_chain="B" max_trajectories=50000 max_hallucinated_trajectories=10000 max_passing_designs=1000 &

done
wait
