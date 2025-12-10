source /paperspace/apps/source_conda.sh
conda activate proteinmpnn_binder_design

/paperspace/apps/silent_tools/silentfrompdbs *.pdb > rfd_C3.silent && rm *.pdb ; rm *.trb ; rm -rf traj ; rm nohup.out
