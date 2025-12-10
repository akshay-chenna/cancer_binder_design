source /paperspace/apps/source_conda.sh
conda activate proteinmpnn_binder_design

export CUDA_VISIBLE_DEVICES=-1

mpnn(){
	a=$1
	/paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_top_4.silent -outsilent mpnn_top_4-${a}.silent -temperature 0.0001 -checkpoint_name mpnn_top_4-${a}.point
	sleep 2
}

for i in {1..8}
do
	mpnn $i
done
wait
