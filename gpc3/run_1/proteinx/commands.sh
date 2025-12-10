# Get the continuous chain
pdb_tidy fold_glypican_3_model_0.pdb | pdb_selres -59:477 | pdb_tofasta >> gpc3_59to477.fasta ## Residue 435 missing !!!!
protenix msa --input gpc3_59to477.fasta --out_dir .
mv 0.a3m gpc3_59to477.a3m
mv 0 msa_gpc3_59to477

# Prepare jsons
cp bindcraft_passed/*pdb bindcraft_passed_forprotenix/
cd bindcraft_passed_forprotenix
for i in *.pdb ; do a=`echo $i | cut -d _ -f 4-6` ; mv $i ${a}.pdb ; done
for i in *pdb ; do protenix tojson --input $i --out_dir input_jsons ; done

for i in *.json ; do a=`echo $i | cut -d - -f 1` ; mv $i ${a}.json ; done

# Make msa
for i in *.json
do
	a=`jq -r '.[0].sequences[1].proteinChain.sequence' ${i}`
	b=`echo ${i} | cut -d . -f 1`
	mkdir msa_binder_${b}
	echo -e ">query\n${a}" >> msa_binder_${b}/pairing.a3m
	echo -e ">query\n${a}" >> msa_binder_${b}/non_pairing.a3m
done

#Modify json
for i in *json
do
	a=`echo $i | cut -d . -f 1`
	jq --arg binder_dir "../msas/msa_binder_${a}" '.[0].sequences[0].proteinChain += {"msa": {"precomputed_msa_dir": "../msas/msa_gpc3_59to477/", "pairing_db": "uniref100"}} | .[0].sequences[1].proteinChain += {"msa": {"precomputed_msa_dir": $binder_dir, "pairing_db": "uniref100"}}' $i > msa_${i}
done
mv msa_l*json ../../input_jsons_withmsa/.
cd input_jsons_withmsa
sed -i "s/${b}/${a}/g" *json


# Run protenix with AA files

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

# Compute ipSAE, actifipTM, LIS
for i in l*
do
	cd ${i}/seed_101/predictions/
	python /paperspace/apps/ipsae_modified.py *data_sample_0.json *sample_0.cif 10 10 "protenix"
	cd -
done

grep max l*/seed_101/predictions/*10.txt | awk '{print $6, $9, $13, $NF}'  >> protenix_ipsaemax_actifiptm_lis.txt


# Compute rosetta_metrics 
python /paperspace/apps/de_novo_binder_scoring/scripts/compute_rosetta_metrics.py --folder px:outputs_aa/ --out-csv rosetta_px_metrics.csv
