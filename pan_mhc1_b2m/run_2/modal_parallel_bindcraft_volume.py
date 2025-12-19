"""
Modal script to run BindCraft with parallelization.

VOLUME_NAME=test GPU=H100 TIMEOUT_HOURS=24 modal run modal_parallel_bindcraft_volume.py \
    --settings \
    --filters \
    --advanced \
    --pdb-in \
    --num-gpus # number of GPUs to use in parallel

modal volume get $VOLUME_NAME / out_1 # download modal volume
"""

import os
from pathlib import Path
from modal import App, Image, Volume

GPU = os.environ.get("GPU", "A100")
TIMEOUT = int(os.environ.get("TIMEOUT_HOURS", "24")) * 60 * 60  # in seconds
VOLUME_NAME = os.environ.get("VOLUME_NAME", "bindcraft_volume")
print(f"Using GPU {GPU} for a maximum of {TIMEOUT} seconds")

image = (
    Image.debian_slim()
    .micromamba()
    .apt_install(["git","wget","mlocate","libquadmath0"])
    .run_commands("git clone https://github.com/martinpacesa/BindCraft.git")
    .workdir("BindCraft")
    .add_local_file("install_bindcraft.sh", remote_path="/tmp/BindCraft/install_bindcraft.sh", copy=True)
    .run_commands("bash install_bindcraft.sh --cuda '12.4' --pkg_manager 'micromamba'", gpu="A100")
    .apt_install("libgfortran5")
    .workdir("/root")
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
        "micromamba", "run", "-n", "BindCraft",
        "python", "-u", "/tmp/BindCraft/bindcraft.py",
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