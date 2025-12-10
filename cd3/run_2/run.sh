source ~/apps/source_conda.sh
conda deactivate
conda activate BindCraft

<< 'END'
export XLA_PYTHON_CLIENT_MEM_FRACTION=.23
for i in {1..4}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './peptide_filters.json' --advanced './peptide_3stage_multimer_weightsconinter5.json' &
done

for i in {5..8}
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './peptide_filters.json' --advanced './peptide_3stage_multimer_weightsconinter20.json' &
done
wait

export XLA_PYTHON_CLIENT_MEM_FRACTION=.23
for i in {9..12}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './peptide_relaxed_filters.json' --advanced './peptide_custom_1.json' &
done

for i in {13..16}
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './peptide_relaxed_filters.json' --advanced './peptide_custom_2.json' &
done
wait
END

export XLA_PYTHON_CLIENT_MEM_FRACTION=.23
for i in {17..20}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './peptide_relaxed_filters.json' --advanced './peptide_custom_3.json' &
done

for i in {21..24}
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './peptide_relaxed_filters.json' --advanced './peptide_custom_4.json' &
done
wait

