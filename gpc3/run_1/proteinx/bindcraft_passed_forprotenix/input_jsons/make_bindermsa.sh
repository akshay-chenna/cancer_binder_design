for i in *.json
do
	a=`jq -r '.[0].sequences[1].proteinChain.sequence' ${i}`
	b=`echo ${i} | cut -d . -f 1`
	mkdir msa_binder_${b}
	echo -e ">query\n${a}" >> msa_binder_${b}/pairing.a3m
	echo -e ">query\n${a}" >> msa_binder_${b}/non_pairing.a3m
done
