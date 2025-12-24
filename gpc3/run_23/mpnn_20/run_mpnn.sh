source /paperspace/apps/source_conda.sh
conda deactivate
conda activate proteinmpnn_binder_design

mpnn(){
	a=$1
	b=$2
	CUDA_VISIBLE_DEVICES=$b /paperspace/apps/dl_binder_design/mpnn_fr/dl_interface_design_soluble.py -silent rfd_mb_${b}.silent -outsilent mpnn_mb_${b}-${a}.silent -temperature 0.1 -checkpoint_name mpnn_mb_${b}-${a}.point
	sleep 2
}

for i in {0..7}
do
	for j in {0..7}
	do
		mpnn $j $i &
	done
done
wait
