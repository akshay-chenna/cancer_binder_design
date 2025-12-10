source /paperspace/apps/source_conda.sh
conda activate proteinmpnn_binder_design

/paperspace/apps/silent_tools/silentfrompdbs top_6_{0..999}.pdb > rfd_top_6.silent && rm *.pdb ; rm *.trb ; rm -rf traj ; rm nohup.out
