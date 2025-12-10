for i in *.pdb
do
	a=`echo $i | cut -d . -f 1`
	python calc_all_rmsds.py $i outputs_template/${a}_model_0.pdb >> rmsds.txt
done
