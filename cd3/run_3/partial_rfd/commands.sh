n=1
while read -r l 
do
pdb_selchain -B ${l}.pdb | pdb_reres -33 | pdb_keepcoord > tmpB.pdb
pdb_selchain -A ${l}.pdb | pdb_keepcoord > tmpA.pdb 
pdb_merge tmpA.pdb tmpB.pdb | pdb_reatom -1 | pdb_tidy > top_${n}.pdb
n=$((n+1))
done < top_visually_selected.txt

for i in {1..18}
do
	cp rfd_binder_partial.sh rfd_binder_partial_${i}.sh
	a=`pdb_selchain -A top_${i}.pdb | pdb_wc | sed -n 3p | awk '{print $3}'`
	sed -i "s/top_1/top_${i}/g" rfd_binder_partial_${i}.sh
	sed -i "s/71/${a}/g" rfd_binder_partial_${i}.sh
done	

for i in {9..18} ; do sed -i 's/CUDA_VISIBLE_DEVICES=0/CUDA_VISIBLE_DEVICES=1/g' rfd_binder_partial_${i}.sh ; done
