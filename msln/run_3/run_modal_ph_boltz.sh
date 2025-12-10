VOLUME_NAME=msln_run3

VOLUME_NAME=${VOLUME_NAME} TIMEOUT_MINUTES=1400 GPU=H100 uvx modal run modal_proteinhunter.py \
	--protein-seqs KTACPSGKKAREIDESLIFYKKWELEACVDAALLATQMDRVNAIPFTYEQLDVLKHKLDELYPQGYPESVIQHLGYLFLKMSPEDIRKWNVTSLETLKALLEVNKGHEMSPQVATLIDRFVKGRGQLDKDTLDTLTAFYPGYLCSLSPEELSSVPPSSIWAVRPQDLDTCDPRQLDVLYPKARLAFQNMNGSEYFVKIQSFLGGAPTEDLKALSQQNVSMDLATFMKLRTDAVLPLTVAEVQKLL \
	--protein-ids B \
	--protein-msas "" \
	--template-path af3_MSLN.cif \
	--template-cif-chain-id A \
	--binder-chain A \
	--min-design-protein-length 50 \
	--max-design-protein-length 150 \
	--high-iptm-threshold 0.7 \
        --contact-residues "18,19,23,47,76,79,109,138,160,184,189,222,226,234,245" \
        --use-potentials \
	--num-designs 1000 \
	--num-cycles 20 \
	--num-gpus 50 \
	--name ${VOLUME_NAME}

uvx modal volume get $VOLUME_NAME / volume_data --force
