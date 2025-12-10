source /paperspace/apps/source_conda.sh
conda deactivate
conda activate BindCraft

export XLA_PYTHON_CLIENT_MEM_FRACTION=.47
for i in {4..4}
do
	for j in {2..2}
	do
		CUDA_VISIBLE_DEVICES=$i python /paperspace/apps/BindCraft/bindcraft.py --settings "./trac_${j}.json" --filters './default_filters.json' --advanced './hotspots10_4stage_multimer.json' &
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
