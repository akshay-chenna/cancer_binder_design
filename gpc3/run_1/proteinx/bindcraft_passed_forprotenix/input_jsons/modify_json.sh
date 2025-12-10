for i in *json
do 
	a=`echo $i | cut -d . -f 1`	
	jq --arg binder_dir "../msas/msa_binder_${a}" '.[0].sequences[0].proteinChain += {"msa": {"precomputed_msa_dir": "../msas/msa_gpc3_59to477/", "pairing_db": "uniref100"}} | .[0].sequences[1].proteinChain += {"msa": {"precomputed_msa_dir": $binder_dir, "pairing_db": "uniref100"}}' $i > msa_${i}
done
