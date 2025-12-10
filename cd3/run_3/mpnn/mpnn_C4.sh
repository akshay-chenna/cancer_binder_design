source ~/scripts/source_conda.sh
conda activate proteinmpnn_binder_design

export CUDA_VISIBLE_DEVICES=1

mpnn(){
	a=$1
	~/apps/dl_binder_design/mpnn_fr/dl_interface_design.py -silent rfd_C4.silent -outsilent mpnn_C4_${a}.silent -temperature 0.0001 -checkpoint_name mpnn_C4_${a}.point
	sleep 2
}

for i in {1..8}
do
	mpnn $i &
done
wait

