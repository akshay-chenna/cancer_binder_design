source ~/apps/source_conda.sh
conda deactivate
conda activate BindCraft

export XLA_PYTHON_CLIENT_MEM_FRACTION=.18
for i in {1..5}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './peptide_relaxed_filters.json' --advanced './peptide_3stage_multimer_weightsconinter10.json' &
done
wait
