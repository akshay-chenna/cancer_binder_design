source ~/apps/source_conda.sh
conda deactivate
conda activate BindCraft

export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
for i in {1..2}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "hla0201_p53r175h.json" --filters './relaxed_filters.json' --advanced './default_4stage_multimer_hardtarget.json' &
done
wait

