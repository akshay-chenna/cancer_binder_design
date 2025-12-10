source ~/apps/source_conda.sh
conda deactivate
conda activate BindCraft

export XLA_PYTHON_CLIENT_MEM_FRACTION=.23
for i in {1..4}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hla_peptide_${i}.json" --filters './default_filters.json' --advanced './hotspots10w_4stage_multimer_flexible_hardtarget.json' &
done
wait

