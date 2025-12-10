source /paperspace/apps/source_conda.sh
conda deactivate
conda activate germinal

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
for i in {0..0}
do

CUDA_VISIBLE_DEVICES=$i python run_germinal.py experiment_name="vhh_gpc3_hotspots_1_${i}" target.target_hotspots=\"A301,A467,A471,A393,A400\" target.target_name="gpc3_h1" target.target_pdb_path="pdbs/correct_uniprot_gpc3_hotspot1_renum.pdb" target.target_chain="A" target.binder_chain="B" max_trajectories=50000 max_hallucinated_trajectories=10000 max_passing_designs=1000 &

done
wait
