source ~/apps/source_conda.sh
conda deactivate
conda activate BindCraft

export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
for i in {0..1}
do
	for j in {1..2}
	do
		CUDA_VISIBLE_DEVICES=$i python ~/apps/BindCraft/bindcraft.py --settings "./gpc3_1.json" --filters './default_filters.json' --advanced './default_4stage_multimer.json' &
	done
done
wait

<< 'END'
for i in {1..1}
do
CUDA_VISIBLE_DEVICES=0 python ~/apps/BindCraft/bindcraft.py --settings "./gpc3_2.json" --filters './default_filters.json' --advanced './default_4stage_multimer_hardtarget.json' &
CUDA_VISIBLE_DEVICES=1 python ~/apps/BindCraft/bindcraft.py --settings "./gpc3_2.json" --filters './default_filters.json' --advanced './default_4stage_multimer_hardtarget.json' &
done
wait
END
