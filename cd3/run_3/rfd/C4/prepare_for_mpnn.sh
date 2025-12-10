source ~/scripts/source_conda.sh
conda activate proteinmpnn_binder_design

/home/paperspace/apps/silent_tools/silentfrompdbs C4_{0..2500}.pdb > rfd_C4.silent
