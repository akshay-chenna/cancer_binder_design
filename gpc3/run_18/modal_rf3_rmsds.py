"""
Modal RMSD calculations using biotite

VOLUME_NAME=test TIMEOUT_HOURS=24 modal run modal_rfd3_rmsd.py \
    --n-mpnn \
    --fixed-chains \
    --designable-chains \
    --num-designs \
    --num-jobs # number of jobs to use in parallel

modal volume get $VOLUME_NAME / out_1 # download modal volume
"""

import os
from pathlib import Path 
from modal import App, Image, Volume 

TIMEOUT = int(os.environ.get("TIMEOUT_HOURS", "24")) * 60 * 60  # in seconds
VOLUME_NAME = os.environ.get("VOLUME_NAME","rmsd_volume")

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
    
    from biotite.structure.io import load_structure
    from biotite.structure import rmsd, superimpose

app = App("RMSD", image=image)
volume = Volume.from_name(VOLUME_NAME, create_if_missing=False)

@app.function(timeout=TIMEOUT, volumes={"/data": volume})
def run_parallel(
    n_mpnn: int,
    fixed_chains: str,
    designable_chains: str,
    num_designs_per_job: int,
    job_idx: int):

    def rmsd_calc(
        fixed_chains,
        designable_chains,
        struct_idx,
        job_idx,
        mpnn_idx):

        idx = f"job_{job_idx}-struct_{struct_idx}"
        folder_path = f"/data/job_{job_idx}"

        aa_generated = load_structure(f"{folder_path}/{idx}_generated.cif")
        aa_refolded = load_structure(f"{folder_path}/{idx}_mpnn_{mpnn_idx}_refolded.cif")

        bb_generated = aa_generated[np.isin(aa_generated.atom_name,["N","CA","C","O"])]
        bb_refolded = aa_refolded[np.isin(aa_refolded.atom_name,["N","CA","C","O"])]

        binder_generated = aa_generated[np.isin(aa_generated.chain_id,designable_chains)]
        binder_refolded = aa_refolded[np.isin(aa_refolded.chain_id,designable_chains)]

        bb_binder_generated = binder_generated[np.isin(binder_generated.atom_name, ["N","CA","C","O"])]
        bb_binder_refolded = binder_refolded[np.isin(binder_refolded.atom_name, ["N","CA","C","O"]) ]
        
        bb_target_generated = bb_generated[np.isin(bb_generated.chain_id,fixed_chains)]
        bb_target_refolded = bb_refolded[np.isin(bb_refolded.chain_id,fixed_chains)]

        bb_binder_refolded_binder_aligned, _ = superimpose(bb_binder_generated, bb_binder_refolded)
        bb_binder_sc_rmsd = rmsd(bb_binder_generated, bb_binder_refolded_binder_aligned)

        bb_target_refolded_target_aligned, transformation = superimpose(bb_target_generated, bb_target_refolded)
        bb_binder_refolded_target_aligned = transformation.apply(bb_binder_refolded)
        bb_binder_target_aligned_rmsd = rmsd(bb_binder_generated, bb_binder_refolded_target_aligned)

        df = pd.DataFrame([{'bb_sc_rmsd': bb_binder_sc_rmsd, 'bb_target_aligned_rmsd': bb_binder_target_aligned_rmsd}])
        df.to_csv(f"{folder_path}/{idx}_rmsd.csv", index=False,mode='a')


    for struct_idx in range(num_designs_per_job):
        for mpnn_idx in range(n_mpnn):
            rmsd_calc(
                fixed_chains,
                designable_chains,
                struct_idx,
                job_idx,
                mpnn_idx)

@app.local_entrypoint()
def main(
    n_mpnn: int = 8,
    fixed_chains: str = "B",
    designable_chains: str = "A",
    num_designs: int = 1000,
    num_gpus: int = 1,
):
    
    num_jobs = num_gpus
    num_designs_per_job = num_designs // num_jobs

    for outputs in run_parallel.starmap(
        [
            (
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
        print("Completed RMSD job(s).")