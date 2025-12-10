for i in l*
do
	cd ${i}/seed_101/predictions/
	python /paperspace/apps/ipsae_modified.py *data_sample_0.json *sample_0.cif 10 10 "protenix"
	cd -
done

grep max l*/seed_101/predictions/*10.txt | awk '{print $6, $9, $13, $NF}'  >> protenix_ipsaemax_actifiptm_lis.txt
