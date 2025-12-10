#for j in {1..18} ; do for i in {1..8} ; do silentsequence af2_top_${j}-${i}.silent >> af2_top_${j}-${i}.sequences  & done ; done
#wait
#for j in {1..18} ; do for i in {1..8} ; do sed -i "s/$/_${i}/g" af2_top_${j}-${i}.sequences  & done ; done
#wait
#cat *.sequences >> all.sequences
for j in {1..18}; do for i in {1..8}; do sed -i "s/$/_${i}/g" af2_top_${j}-${i}.sc; done; done
cat *sc >> all.sc
