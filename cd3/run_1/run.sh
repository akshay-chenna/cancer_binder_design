source ~/apps/source_conda.sh
conda deactivate
conda activate BindCraft

<<'END'
python ~/apps/BindCraft/bindcraft.py --settings './run1.json' --filters './default_filters.json' --advanced './default_4stage_multimer.json'

export XLA_PYTHON_CLIENT_MEM_FRACTION=.20
for i in {1..4}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './default_4stage_multimer.json' &
done

for i in {5..8}
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './default_4stage_multimer.json' &
done
wait

export XLA_PYTHON_CLIENT_MEM_FRACTION=.22
for i in {9..12}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './hotspots10_4stage_multimer.json' &
done

for i in {13..16}
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './hotspots10_4stage_multimer.json' &
done
wait
END

export XLA_PYTHON_CLIENT_MEM_FRACTION=.22
for i in {17..20}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './hotspots20_4stage_multimer.json' &
done

for i in {21..24}
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './hotspots50_4stage_multimer.json' &
done
wait
