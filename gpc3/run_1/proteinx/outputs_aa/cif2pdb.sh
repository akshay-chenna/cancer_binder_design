for i in l*
do
	cd ${i}/seed_101/predictions/
	for i in *sample_0.cif ; do a=`echo $i | cut -d . -f 1` ; obabel -icif $i -opdb -O ${a}.pdb ; done
	cd -
done


