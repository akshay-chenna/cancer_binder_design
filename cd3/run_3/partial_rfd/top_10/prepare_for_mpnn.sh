source ~/scripts/source_conda.sh
conda activate proteinmpnn_binder_design

/home/paperspace/apps/silent_tools/silentfrompdbs top_10_{0..999}.pdb > rfd_top_10.silent && rm *.pdb ; rm *.trb ; rm nohup.out ; rm -rf traj
