mkdir outputs_aa

for i in dir_1/*
do
	CUDA_VISIBLE_DEVICES=0 python /paperspace/apps/Protenix/runner/inference.py --dump_dir ./outputs_aa --input_json_path $i --need_atom_confidence true 
done &

for i in dir_2/*
do
	CUDA_VISIBLE_DEVICES=1 python /paperspace/apps/Protenix/runner/inference.py --dump_dir ./outputs_aa --input_json_path $i --need_atom_confidence true 
done &

for i in dir_3/*
do
	CUDA_VISIBLE_DEVICES=2 python /paperspace/apps/Protenix/runner/inference.py --dump_dir ./outputs_aa --input_json_path $i --need_atom_confidence true 
done &

for i in dir_4/*
do
	CUDA_VISIBLE_DEVICES=3 python /paperspace/apps/Protenix/runner/inference.py --dump_dir ./outputs_aa --input_json_path $i --need_atom_confidence true 
done &
wait
mv outputs_aa ../.
