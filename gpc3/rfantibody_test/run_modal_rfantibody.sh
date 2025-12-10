TIMEOUT_MINUTES=1400 GPU=H100 uvx modal run modal_rfantibody.py \
  	--contigs '[L1:9-14,L2:1-5,L3:7-13,H1:5-10,H2:5-12,H3:5-13]' \
  	--run-name gpc3_run20 \
  	--num-seqs 8 \
  	--num-gpus 2 \
  	--num-designs 4 \
  	--input-framework-filepath ./yp7_HLT.pdb \
  	--input-target-filepath ./correct_uniprot_gpc3_59-477_renum_rfa.pdb \
	--hotspot-res '[T67,T71,T147,T153,T156,T275]'
#	--n-iterations-per-design INTEGER \

#uvx modal volume get $VOLUME_NAME / volume_data --force
