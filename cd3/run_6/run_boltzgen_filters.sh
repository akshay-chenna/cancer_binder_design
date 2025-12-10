source /paperspace/apps/source_conda.sh
conda deactivate
conda activate bg

cd ./out/boltzgen/2511021909/combined/intermediate_designs_inverse_folded

NAME=cd3e_mb

#Rename files
rename "s/${NAME}_(\\d+)\\.cif/sprintf(\"${NAME}_%04d.cif\", \$1)/e" ${NAME}_*.cif
rename "s/${NAME}_(\\d+)\\.npz/sprintf(\"${NAME}_%04d.npz\", \$1)/e" ${NAME}_*.npz

cd fold_out_design_npz
rename "s/${NAME}_(\\d+)\\.npz/sprintf(\"${NAME}_%04d.npz\", \$1)/e" ${NAME}_*.npz
cd ..

cd refold_design_cif
rename "s/${NAME}_(\\d+)\\.cif/sprintf(\"${NAME}_%04d.cif\", \$1)/e" ${NAME}_*.cif
cd ..

cd fold_out_npz
rename "s/${NAME}_(\\d+)\\.npz/sprintf(\"${NAME}_%04d.npz\", \$1)/e" ${NAME}_*.npz
cd ..

cd refold_cif
rename "s/${NAME}_(\\d+)\\.cif/sprintf(\"${NAME}_%04d.cif\", \$1)/e" ${NAME}_*.cif
cd ..

cd ~/cancer/cd3/run_6

CUDA_VISIBLE_DEVICES=0 boltzgen run ./cd3e_mb.yaml  --output ./out/boltzgen/2511021909/combined --protocol protein-anything --budget 1000 --steps analysis filtering
