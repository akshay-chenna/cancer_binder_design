sort -t ',' -nk 7 *csv | awk -F , '{print $1, $7}' | awk '!seen[$0]++' >> iptm.txt
