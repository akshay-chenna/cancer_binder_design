source ~/scripts/source_conda.sh
conda activate proteinmpnn_binder_design

export CUDA_VISIBLE_DEVICES=1

mpnn(){
	a=$1
	~/apps/dl_binder_design/mpnn_fr/dl_interface_design.py -silent rfd_C2.silent -outsilent mpnn_C2_${a}.silent -temperature 0.0001 -checkpoint_name mpnn_C2_${a}.point
	sleep 2
}

for i in {1..8}
do
	mpnn $i &
done
wait

