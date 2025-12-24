source ~/apps/source_conda.sh
conda deactivate
conda activate proteinmpnn_binder_design

for j in {1..8}
do
for i in -${j}.silent
do
	silentls $i | sed "s/$/_${j}/g" | silentrename $i > renamed_$i
done
done

for i in renamed_*silent
do
	silentextract $i &
done
wait

mkdir gpc3_run_23_af2_20
mv mb*pdb gpc3_run_23_af2_20
