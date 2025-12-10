VOLUME_NAME=tp53_run6

VOLUME_NAME=${VOLUME_NAME} TIMEOUT_MINUTES=1400 GPU=H100 uvx modal run modal_proteinhunter.py \
	--protein-seqs HMTEVVRHC,MGSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTVQRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQ \
	--protein-ids C,A \
	--protein-msas "" \
	--template-path 6w51_cleaned.cif,6w51_cleaned.cif \
	--template-cif-chain-id C,A \
	--binder-chain B \
	--min-design-protein-length 50 \
	--max-design-protein-length 150 \
	--high-iptm-threshold 0.7 \
	--contact-residues "7,8,9" \
	--use-potentials \
	--num-designs 1000 \
	--num-cycles 20 \
	--num-gpus 50 \
	--name ${VOLUME_NAME}

uvx modal volume get $VOLUME_NAME / volume_data --force
