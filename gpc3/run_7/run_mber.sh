source /paperspace/apps/source_conda.sh
conda deactivate
conda activate mber

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
for i in {0..7}
do
	CUDA_VISIBLE_DEVICES=${i} mber-vhh --input-pdb ./correct_uniprot_gpc3_hotspot2.pdb --output-dir ./run_${i} --target-name gpc3 --chains A --hotspots 'A301,A467,A471,A393,A400' --num-accepted 1000  --max-trajectories 100000 --min-iptm 0.75 --min-plddt 0.7  --no-animations  --no-pickle --no-png &

done
wait
