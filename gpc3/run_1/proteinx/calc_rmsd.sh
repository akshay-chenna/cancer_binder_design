for i in *.pdb
do
	a=`echo $i | cut -d . -f 1`
	python calc_all_rmsds.py $i outputs_aa/${a}_sample_0.pdb >> rmsds.txt
done
