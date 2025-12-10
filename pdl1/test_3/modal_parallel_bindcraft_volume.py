"""
Modal script to run BindCraft with parallelization.

VOLUME_NAME=test GPU=H100 TIMEOUT_HOURS=24 modal run modal_parallel_bindcraft_volume.py \
    --settings \
    --filters \
    --advanced \
    --pdb-in \
    --num-gpus # number of GPUs to use in parallel
"""

import os
from pathlib import Path
import modal
from modal import App, Image, Volume
import sys

GPU = os.environ.get("GPU", "A100")
TIMEOUT = int(os.environ.get("TIMEOUT_HOURS", "24")) * 60 * 60  # in seconds
VOLUME_NAME = os.environ.get("VOLUME_NAME", "bindcraft_volume")
print(f"Using GPU {GPU} for a maximum of {TIMEOUT} seconds")

image = (
    Image.debian_slim(python_version="3.11")
    .apt_install("git", "wget", "aria2", "ffmpeg")
    .pip_install("numpy<2.0")  # Pin NumPy FIRST before any dependencies
    .pip_install(
        "pdb-tools==2.4.8", "ffmpeg-python==0.2.0", "plotly==5.18.0", "kaleido==0.2.1"
    )
    .pip_install("git+https://github.com/sokrypton/ColabDesign.git")
    .run_commands(
        "git clone https://github.com/martinpacesa/BindCraft /root/bindcraft",
        "cd /root/bindcraft && git checkout c0a48d595d4976694aa979438712ac94c16620bb",
        "chmod +x /root/bindcraft/functions/dssp",
        "chmod +x /root/bindcraft/functions/DAlphaBall.gcc",
    )
    .run_commands(
        "ln -s /usr/local/lib/python3.*/dist-packages/colabdesign colabdesign && mkdir /params"
    )
    .run_commands(
        "aria2c -q -x 16 https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"
        " && mkdir -p /root/bindcraft/params"
        " && tar -xf alphafold_params_2022-12-06.tar -C /root/bindcraft/params"
    )
    .pip_install(
        "numpy<2.0",  # Re-enforce after pyrosetta (which may upgrade it)
        "jax[cuda]<0.7.0",  # Pin to avoid 'wraps' removal in JAX 0.7.0
        "matplotlib==3.8.1",  # https://github.com/martinpacesa/BindCraft/issues/4
    )
	.apt_install("wget")
    .run_commands("wget https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Debug.python311.linux.wheel/pyrosetta-2024.39+release.59628fb-cp311-cp311-linux_x86_64.whl")
    .run_commands("pip install pyrosetta-2024.39+release.59628fb-cp311-cp311-linux_x86_64.whl")
)

with image.imports():
    import json
    import subprocess

app = App("Bindcraft", image=image)
volume = Volume.from_name(VOLUME_NAME, create_if_missing=True)

@app.function(gpu=GPU, timeout=TIMEOUT, volumes={"/data":volume}) 
def run_bindcraft(
    settings_str: str,
    filters_str: str,
    advanced_str: str,
    setting_name: str,
    filter_name: str,
    advanced_name: str,
    pdb_str: str,
    pdb_name: str,
    job_idx: int,
):
    from pathlib import Path

    in_dir_path = Path(f"/data/job_{job_idx}")
    # Ensure the directory exists
    in_dir_path.mkdir(exist_ok=True, parents=True)

    settings_path = in_dir_path / setting_name
    filters_path = in_dir_path / filter_name
    advanced_path = in_dir_path / advanced_name
    pdb_path = in_dir_path / pdb_name
    
    settings_json = json.loads(settings_str)
    settings_json["starting_pdb"] = str(pdb_path)
    settings_json["design_path"] = str(in_dir_path)
    settings_str = json.dumps(settings_json)
    
    settings_path.write_text(settings_str)
    filters_path.write_text(filters_str)
    advanced_path.write_text(advanced_str)
    pdb_path.write_text(pdb_str)
    print(f"Wrote PDB to {pdb_path}")

    cmd = [
        "python", "-u", "/root/bindcraft/bindcraft.py",
        "--settings", str(settings_path),
        "--filters", str(filters_path),
        "--advanced", str(advanced_path),
    ]
    subprocess.run(cmd, check=True)
    
@app.local_entrypoint()
def main(
    settings: str,
    filters: str,
    advanced: str,
    pdb_in: str,
    num_gpus: int = 1,
    out_dir: str = "./out/bindcraft",
):

    settings_path = Path(settings)
    settings_str = settings_path.read_text()
    filters_path = Path(filters)
    filters_str = filters_path.read_text()
    advanced_path = Path(advanced)
    advanced_str = advanced_path.read_text()
    pdb_path = Path(pdb_in)
    pdb_str = pdb_path.read_text()
    
    num_jobs = num_gpus
    
    for outputs in run_bindcraft.starmap(
        [ 
         (settings_str,
          filters_str,
          advanced_str,
          settings_path.name,
          filters_path.name,
          advanced_path.name,
          pdb_str,
          pdb_path.name,
          i,
           ) 
         for i in range(num_jobs)
         ],
        order_outputs = False,
    ):
        print("Completed BindCraft job(s).")
    print(f"All BindCraft jobs completed. Outputs saved to the Volume")