1. run_1: Bindcraft Minibinder design against the CD3e epitopes
	run1.json; outputs: 4 hotspots used: 36,47,54,82
	hotspot{1..8}.json: single hotspots used -- 50,78,53,108,36,47,54,82
	hotspot{9..16}.json: single hotspots used -- 50,78,53,108,36,47,54,82 -- higher hotspots weights (weights_con_inter) of 10: hotspots10_4stage_multimer.json
	hotspots{17..20}.json: All hotspots used -- 36,47,50,53,54,78,82,108 -- higher hotspots weights of 20: hotspots20_4stage_multimer.json
	hotspots{21..24}.json: All hotspots used -- 36,47,50,53,54,78,82,108 -- higher hotspots weights of 50: hotspots50_4stage_multimer.json
2. run_2: Bindcraft peptide design
	hotspots{1..4}.json: Peptide design lengths=10-30. Hotspots= 36,50,53,78. weights_con_inter=5 (peptide_3stage_multimer_weightsconinter5.json). Filters=peptide_filters.json
	hotspots{5..8}.json: Peptide design lengths=10-30. Hotspots= 36,50,53,78. weights_con_inter=20 (peptide_3stage_multimer_weightsconinter20.json). Filters=peptide_filters.json
	hotspots{9..12}.json: Peptide design lengths=10-30. Hotspots= 36,50,53,78. .Settings=peptide_custom_1.json. Filters=peptide_relaxed_filters.json
		<     "weights_con_intra": 0.5,
		<     "weights_con_inter": 0.5,
		---
		>     "weights_con_intra": 0.0,
		>     "weights_con_inter": 50,
		29,30c29,30
		<     "weights_helicity": 0.95,
		<     "random_helicity": false,
		---
		>     "weights_helicity": 0.0,
		>     "random_helicity": true,
	hotspots{13..16}.json: Peptide design lengths=10-30. Hotspots= 36,50,53,78. .Settings=peptide_custom_2.json. Filters=peptide_relaxed_filters.json
		<     "weights_con_intra": 0.5,
		<     "weights_con_inter": 0.5,
		---
		>     "weights_con_intra": 0.0,
		>     "weights_con_inter": 50,
		29,30c29,30
		<     "weights_helicity": 0.95,
		<     "random_helicity": false,
		---
		>     "weights_helicity": 0.0,
		>     "random_helicity": true,
		<     "use_termini_distance_loss": false,
		---
		>     "use_termini_distance_loss": true,
	hotspots{17..20}.json
		24c24
		<     "weights_con_inter": 0.5,
		---
		>     "weights_con_inter": 50,
		29c29
		<     "weights_helicity": 0.95,
		---
		>     "weights_helicity": 2,
		44c44
		<     "mpnn_weights": "soluble",
		---
		>     "mpnn_weights": "original",
	hotspots{21..24}.json
		24c24
		<     "weights_con_inter": 0.5,
		---
		>     "weights_con_inter": 50,
		29,30c29,30
		<     "weights_helicity": 0.95,
		<     "random_helicity": false,
		---
		>     "weights_helicity": 0.0,
		>     "random_helicity": true,
		44c44
		<     "mpnn_weights": "soluble",
		---
		>     "mpnn_weights": "original",
3. run_3: RFDiffusion binder design
4. run_4: FoldCraft binder design
5. run_5: Germinal VHH
6. run_6: Modal boltzgen CD3e minibinder design using boltzgen with hotspots from run_3. cd3e_cl_00.pdb renumbered to cd3e_cl_00_renum.pdb
7. run_7: Modal boltzgen CD3e nanobody design using boltzgen with hotspots from run_3. cd3e_cl_00.pdb renumbered to cd3e_cl_00_renum.pdb
