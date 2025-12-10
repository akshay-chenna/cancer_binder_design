source /paperspace/apps/source_conda.sh
conda activate proteinmpnn_binder_design

/paperspace/apps/silent_tools/silentfrompdbs top_18_{0..999}.pdb > rfd_top_18.silent && rm *.pdb ; rm *.trb ; rm -rf traj ; rm nohup.out
