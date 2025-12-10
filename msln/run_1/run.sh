source ~/apps/source_conda.sh
conda deactivate
conda activate BindCraft

<< 'END'
export XLA_PYTHON_CLIENT_MEM_FRACTION=.18
for i in {1..5}
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './default_4stage_multimer.json' &
done
wait

export XLA_PYTHON_CLIENT_MEM_FRACTION=.20
for i in {6..9}
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './hotspots10_4stage_multimer.json' &
done
wait

export XLA_PYTHON_CLIENT_MEM_FRACTION=.30
for i in 4 5 10
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './hotspots10_4stage_multimer.json' &
done
wait
END

export XLA_PYTHON_CLIENT_MEM_FRACTION=.30
for i in 6 7
do
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./hotspots_${i}.json" --filters './default_filters.json' --advanced './hotspots10_4stage_multimer.json' &
done
wait
