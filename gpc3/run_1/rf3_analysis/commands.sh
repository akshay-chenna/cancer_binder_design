# merge all accepted score files
cat outputs_template/*csv | awk '!seen[$0]++' >> bindcraft_passed-rf3template_metrics.csv 

cat ../out_1/final_design_stats.csv | awk '!seen[$0]++' >> bindcraft_passed-bindcraft_metrics.csv

#convert cif.gz to pdb
for i in *model_0.cif.gz ; do a=`echo $i | cut -d . -f 1` ; obabel -icif $i -opdb -O ${a}.pdb ; done

bash calc_rmsd.sh
python /paperspace/apps/de_novo_binder_scoring/scripts/compute_rosetta_metrics.py --folder rf3:outputs_template/ --out-csv rosetta_rf3_metrics.csv
python /paperspace/apps/de_novo_binder_scoring/scripts/compute_rosetta_metrics.py --folder af2:bindcraft_passed/ --out-csv rosetta_af2_metrics.csv
