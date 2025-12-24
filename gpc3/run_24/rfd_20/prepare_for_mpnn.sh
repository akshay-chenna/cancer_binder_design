source /paperspace/apps/source_conda.sh
conda deactivate
conda activate proteinmpnn_binder_design

for i in {0..7}
do	
	cd mb_${i}
	/paperspace/apps/silent_tools/silentfrompdbs mb_${i}_{0..124}.pdb > rfd_mb_${i}.silent && rm *.pdb ; rm *.trb ; rm -rf traj &
	cd ..
done
