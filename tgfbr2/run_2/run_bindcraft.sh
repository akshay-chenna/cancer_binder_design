source /paperspace/apps/source_conda.sh
conda deactivate
conda activate BindCraft

#export XLA_PYTHON_CLIENT_MEM_FRACTION=.48
for i in {1..1}
do
	for j in {0..0}
	do
		CUDA_VISIBLE_DEVICES=$i python /paperspace/apps/BindCraft/bindcraft.py --settings "./tgfbr2_1.json" --filters './default_filters.json' --advanced './betasheet_4stage_multimer.json' &
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
