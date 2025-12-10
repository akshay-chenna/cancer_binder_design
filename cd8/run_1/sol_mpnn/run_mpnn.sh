for i in {3..4}
do
	bash sol_mpnn_C${i}.sh &
done
wait

<< 'END'
source /paperspace/apps/source_conda.sh
conda activate proteinmpnn_binder_design

export CUDA_VISIBLE_DEVICES=-1

/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_13.silent -outsilent mpnn_top_13-3.silent -temperature 0.0001 -checkpoint_name mpnn_top_13-3.point &

/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_13.silent -outsilent mpnn_top_13-6.silent -temperature 0.0001 -checkpoint_name mpnn_top_13-6.point &

/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_13.silent -outsilent mpnn_top_13-8.silent -temperature 0.0001 -checkpoint_name mpnn_top_13-8.point &

/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_13.silent -outsilent mpnn_top_13-7.silent -temperature 0.0001 -checkpoint_name mpnn_top_13-7.point &

/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_15.silent -outsilent mpnn_top_15-8.silent -temperature 0.0001 -checkpoint_name mpnn_top_15-8.point &

/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_15.silent -outsilent mpnn_top_15-7.silent -temperature 0.0001 -checkpoint_name mpnn_top_15-7.point &

/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_4.silent -outsilent mpnn_top_4-8.silent -temperature 0.0001 -checkpoint_name mpnn_top_4-8.point &

/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_5.silent -outsilent mpnn_top_5-8.silent -temperature 0.0001 -checkpoint_name mpnn_top_5-8.point &
wait
END
