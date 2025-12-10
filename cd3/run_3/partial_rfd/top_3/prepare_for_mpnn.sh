source ~/scripts/source_conda.sh
conda activate proteinmpnn_binder_design

/home/paperspace/apps/silent_tools/silentfrompdbs top_3_{0..999}.pdb > rfd_top_3.silent && rm *.pdb ; rm *.trb ; rm nohup.out ; rm -rf traj
