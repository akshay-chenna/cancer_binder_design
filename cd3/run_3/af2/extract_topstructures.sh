cat *.sc | sort -rnk 9 | awk '$9 < 3.5' > top_targetaligned.txt
awk '{print $NF}' top_targetaligned.txt  | sort -u > top_structures.txt

cat top_structures.txt | while read -r l 
do
	a=`grep $l *sc | sort -nk 9 | head -1 | awk '{print $1}' | cut -d . -f 1`
	echo $a >> top_structures_sc.txt
done	

exec 3< top_structures.txt
exec 4< top_structures_sc.txt

while read -u 3 l && read -u 4 m
do
	echo $l
	echo $m
	silentextractspecific ${m}.silent $l
done

mkdir top_structures
mv *pdb top_structures
