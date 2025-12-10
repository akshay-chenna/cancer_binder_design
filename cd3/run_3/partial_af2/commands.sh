sort -rnk 9 all.sc | awk '$9 < 3.5'  > top_designs.sc
awk '{print $NF}' top_designs.sc | while read -r l ; do grep $l all.sequences | awk '{print $3,$1}' >> top_sequences.csv ; done
mv top_sequences.csv cd3_rfd_mb_nov2025_1.csv
