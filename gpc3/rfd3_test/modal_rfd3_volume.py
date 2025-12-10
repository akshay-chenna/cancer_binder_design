"""
Modal design for RFD3

VOLUME_NAME=test GPU=A100 TIMEOUT_HOURS=24 modal run modal_rfd3_volume.py \
    --pdb-in \
    --contig \
    --hotspots \
    --n-mpnn \
    --fixed-chains \
    --designable-chains \
    --num-designs \
    --num-gpus # number of GPUs to use in parallel

modal volume get $VOLUME_NAME / out_1 # download modal volume
"""

import os
from pathlib import Path 
from modal import App, Image, Volume 

GPU = os.environ.get("GPU", "A100")
TIMEOUT = int(os.environ.get("TIMEOUT_HOURS", "24")) * 60 * 60  # in seconds
VOLUME_NAME = os.environ.get("VOLUME_NAME", "rfd3_volume")
print(f"Using GPU {GPU} for a maximum of {TIMEOUT} seconds")

image = (
    Image.debian_slim(python_version="3.12")
    .apt_install("git", "wget", "pip")
    .env({"CCD_MIRROR_PATH": "", "PDB_MIRROR_PATH": ""})
    .run_commands("pip install -q 'rc-foundry[all]'")
    .run_commands("foundry install rfd3 proteinmpnn rf3")
)

with image.imports():
    import json
    import subprocess
    import numpy as np
    import pandas as pd
    from pathlib import Path

    from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
    from mpnn.inference_engines.mpnn import MPNNInferenceEngine
    from rf3.inference_engines.rf3 import RF3InferenceEngine
    from rf3.utils.inference import InferenceInput
    
    from biotite.structure import get_residue_starts
    from biotite.sequence import ProteinSequence
    from biotite.structure import rmsd, superimpose

    from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES
    from atomworks.io.utils.io_utils import to_cif_file


app = App("RFD3", image=image)
volume = Volume.from_name(VOLUME_NAME, create_if_missing=True)

@app.function(gpu=GPU, timeout=TIMEOUT, volumes={"/data": volume})
def run_parallel(
    pdb: str,
    contig: str,
    hotspots: str,
    n_mpnn: int,
    fixed_chains: str,
    designable_chains: str,
    num_designs_per_job: int,
    job_idx: int):

    in_dir_path = Path(f"/data/job_{job_idx}")
    # Ensure the directory exists
    in_dir_path.mkdir(exist_ok=True, parents=True)
    
    pdb_path = in_dir_path / "pdb.pdb"
    pdb_path.write_text(pdb)

    def rfd3_mpnn_rf3(
        pdb_path,
        contig,
        hotspots,
        n_mpnn,
        fixed_chains,
        designable_chains,
        struct_idx,
        job_idx):

        idx = f"job_{job_idx}-struct_{struct_idx}"
        folder_path = f"/data/job_{job_idx}"

        # Run RFD3
        config = RFD3InferenceConfig(
            specification={
                "dialect": 2,
                "infer_ori_strategy": "hotspots",
                "input": str(pdb_path),
                "contig": contig,
                "select_hotspots": hotspots,
                "extra": {},
            },
            diffusion_batch_size=1,
        )

        model = RFD3InferenceEngine(**config)
        outputs = model.run(
            inputs=None,
            out_dir=None,
            n_batches=1,
        )

        first_key = next(iter(outputs.keys()))
        atom_array = outputs[first_key][0].atom_array

        # Run MPNN
        engine_config = {
            'model_type': "protein_mpnn",
            "is_legacy_weights": True,
            "out_directory": "new",
            "write_structures": True,
            "write_fasta": False,
        }

        input_configs = [
            {
                "batch_size": n_mpnn,
                "remove_waters": True,
                "fixed_chains": [fixed_chains],
                "temperature": 0.001,
            }
        ]

        model = MPNNInferenceEngine(**engine_config)
        mpnn_outputs = model.run(input_dicts=input_configs, atom_arrays=[atom_array])

        aa_generated = atom_array              # Original RFD3 backbone (from Section 1)
        chA_generated = aa_generated[np.isin(aa_generated.chain_id,designable_chains)]
        bb_chA_generated = chA_generated[np.isin(chA_generated.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)]
        to_cif_file(aa_generated, f"{folder_path}/{idx}_generated.cif")

        #Run RF3
        def run_rf3(mpnn_atom_array,mpnn_id):
            inference_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)

            input_structure = InferenceInput.from_atom_array(mpnn_atom_array, example_id="binder", template_selection=[fixed_chains])
            rf3_outputs = inference_engine.run(inputs=input_structure)

            rf3_output = rf3_outputs["binder"][0] #Picks the best structure

            #Calculate RMSD
            aa_refolded = rf3_output.atom_array    # RF3-predicted structure
            
            # Export structures to CIF format for visualization in PyMOL/ChimeraX
            to_cif_file(aa_refolded, f"{folder_path}/{idx}_mpnn_{mpnn_id}_refolded.cif")

            chA_refolded = aa_refolded[np.isin(aa_refolded.chain_id,designable_chains)]
            bb_chA_refolded = chA_refolded[np.isin(chA_refolded.atom_name, PROTEIN_BACKBONE_ATOM_NAMES)]

            bb_refolded_fitted, _ = superimpose(bb_chA_generated, bb_chA_refolded)
            rmsd_value = rmsd(bb_chA_generated, bb_chA_refolded)

            res_starts = get_residue_starts(chA_generated)
            seq_1letter = ''.join(
                ProteinSequence.convert_letter_3to1(res_name)
                for res_name in chA_generated.res_name[res_starts]
            )

            rf3_output.summary_confidences['bb_sc_rmsd'] = rmsd_value.item()
            rf3_output.summary_confidences['binder_sequence'] = seq_1letter
            rf3_output.summary_confidences['sequence_length'] = len(seq_1letter)
            rf3_output.summary_confidences['name'] = f"{idx}_mpnn_{mpnn_id}"

            #Save metric file
            df = pd.DataFrame([rf3_output.summary_confidences])
            df.to_csv(f"{folder_path}/{idx}_metrics.csv", index=False,mode='a')

        for mpnn_id in range(n_mpnn):
            run_rf3(mpnn_outputs[mpnn_id].atom_array,mpnn_id)

    for struct_idx in range(num_designs_per_job):
        rfd3_mpnn_rf3(
            pdb_path,
            contig,
            hotspots,
            n_mpnn,
            fixed_chains,
            designable_chains,
            struct_idx,
            job_idx,
        )

@app.local_entrypoint()
def main(
    pdb_in: str,
    contig: str,
    hotspots: str | None = None,
    n_mpnn: int = 8,
    fixed_chains: str = "B",
    designable_chains: str = "A",
    num_designs: int = 1000,
    num_gpus: int = 1,
):
    
    pdb_path = Path(pdb_in)
    pdb_str = pdb_path.read_text()
    num_jobs = num_gpus
    num_designs_per_job = num_designs // num_jobs

    for outputs in run_parallel.starmap(
        [
            (
                pdb_str,
                contig,
                hotspots,
                n_mpnn,
                fixed_chains,
                designable_chains,
                num_designs_per_job,
                job_idx,
            )
            for job_idx in range(num_jobs)
        ],
        order_outputs=False,
    ):
        print("Completed RFD3, MPNN, and RF3 job(s).")