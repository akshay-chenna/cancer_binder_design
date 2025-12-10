from datetime import datetime

import math

import modal

import os

from pathlib import Path

import re

import shutil


VOLUME_NAME = os.environ.get("VOLUME_NAME")
VOLUME_PATH = "/results"
GPU = os.environ.get("GPU", "H100")
TIMEOUT_MINUTES = int(os.environ.get("TIMEOUT_MINUTES", 60))
RFANTIBODY_REPO_PATH = "/root/RFantibody"

if VOLUME_NAME is not None:
    vol = modal.Volume.from_name(VOLUME_NAME)
else:
    vol = None


def download_weights():
    import subprocess
    urls = [
        "https://files.ipd.uw.edu/pub/RFantibody/RFdiffusion_Ab.pt",
        "https://files.ipd.uw.edu/pub/RFantibody/ProteinMPNN_v48_noise_0.2.pt",
        "https://files.ipd.uw.edu/pub/RFantibody/RF2_ab.pt",
    ]
    processes = [
        subprocess.Popen(
            [
                "aria2c",
                "-x",
                "8",
                url,
            ],
        )
        for url in urls
    ]
    for process in processes:
        process.wait()
    print("RFantibody weights downloaded.")


image = (
    modal.Image.debian_slim(python_version = "3.12")
    .apt_install("git", "aria2")
    .run_function(download_weights)
    .uv_pip_install("dgl", find_links = "https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html", gpu = "A10G")
    .uv_pip_install("torch", "hydra-core", "pyrsistent", "opt-einsum", "e3nn", "icecream", gpu = "A10G")
    .run_commands(
        f"git clone https://github.com/RosettaCommons/RFantibody {RFANTIBODY_REPO_PATH}",
        f"pip install --no-deps {RFANTIBODY_REPO_PATH}/include/SE3Transformer",
        # Cannot use uv here because RFantibody uses poetry.
        f"pip install --no-deps -e {RFANTIBODY_REPO_PATH}", 
        f"g++ -static -O3 -ffast-math -lm -o {RFANTIBODY_REPO_PATH}/include/USalign/USalign {RFANTIBODY_REPO_PATH}/include/USalign/USalign.cpp",
    )
)


with image.imports():
    import os
    import subprocess
    import tempfile


app = modal.App("RFantibody", image = image)


@app.function(timeout = TIMEOUT_MINUTES * 60, volumes = {VOLUME_PATH: vol} if vol is not None else {})
def process_and_merge_data(all_data: list[tuple[str, bytes]], out_path: str, designs_per_job: int) -> None:
    out_path = Path(out_path)
    out_path.mkdir(exist_ok = True, parents = True)

    for filepath_str, content in all_data:
        local_filepath = out_path / Path(filepath_str)
        local_filepath.parent.mkdir(parents = True, exist_ok = True)
        local_filepath.write_bytes(content)

    job_dirs = sorted([d for d in out_path.glob('job_*') if d.is_dir()], key = lambda d: int(d.name.split("_")[1]))
    job_ids = [int(d.name.split("_")[1]) for d in job_dirs]
    
    if not job_dirs:
        print(f"No job directories found in {out_path}")
        return
    
    print(f"Found {len(job_dirs)} job directories to process")

    # Deal with structure related files
    structure_file_re = re.compile(r"design_(?P<num>\d+)(?:_(?P<suffix>.+))?\.(?P<ext>.+)")
    for job_id, job_dir in zip(job_ids, job_dirs):
        print(f"Processing {job_dir.name}")
        structure_files = [
            f for ext in ("pdb", )
            for f in job_dir.glob(f"**/*.{ext}")
        ]
        for structure_file in structure_files:
            parent_dir_name = structure_file.parent.name
            match_components = structure_file_re.match(structure_file.name)
            old_num = match_components.group("num")
            new_num = job_id * designs_per_job + int(old_num)
            _suffix = match_components.groupdict().get("suffix")
            if _suffix is not None:
                out_name = f"design_{new_num}_{_suffix}.{match_components.group('ext')}"
            else:
                out_name = f"design_{new_num}.{match_components.group('ext')}"
            out_dir = out_path.joinpath(structure_file.relative_to(job_dir).parent)
            if not out_dir.exists():
                out_dir.mkdir(parents = True)
            shutil.copy(structure_file, out_dir.joinpath(out_name))
    print(f"\nSuccessfully merged all files to {out_path}.")


@app.function(
    gpu = GPU,
    timeout = TIMEOUT_MINUTES * 60,
    volumes = {VOLUME_PATH: vol} if vol is not None else {},
)
def run(
    input_target_data: bytes,
    input_framework_data: bytes,
    job_id: int,
    run_name: str,
    num_designs: int = 500,
    num_seqs: int = 8,
    contigs: str | None = None,
    n_iterations_per_design: int = 50,
    hotspot_res: str | None = None,
) -> list[tuple[str, bytes]]:
    tmp_d_used = False
    if vol is None:
        print("No volume specified. Using temporary directory.")
        out_path = Path(tempfile.mkdtemp()) / f"{run_name}"
        tmp_d_used = True
    else:
        out_path = Path(VOLUME_PATH) / f"{run_name}"
    if not out_path.exists():
        out_path.mkdir(parents = True, exist_ok = True)

    out_path.joinpath("input_target_file.pdb").write_bytes(input_target_data)
    out_path.joinpath("input_framework_file.pdb").write_bytes(input_framework_data)

    out_path.joinpath(f"job_{job_id}", "mpnn").mkdir(parents = True, exist_ok = True)
    out_path.joinpath(f"job_{job_id}", "rf2ab").mkdir(parents = True, exist_ok = True)

    def _backbone_design() -> None:
        cmd = [
            "python", f"{RFANTIBODY_REPO_PATH}/scripts/rfdiffusion_inference.py",
            "--config-path", f"{RFANTIBODY_REPO_PATH}/src/rfantibody/rfdiffusion/config/inference",
            "--config-name", "antibody",
            f"antibody.target_pdb={out_path}/input_target_file.pdb",
            f"antibody.framework_pdb={out_path}/input_framework_file.pdb",
            f"antibody.design_loops={contigs}",
            f"inference.output_prefix={out_path}/job_{job_id}/design",
            f"inference.ckpt_override_path=RFdiffusion_Ab.pt",
            f"inference.num_designs={num_designs}",
            f"diffuser.T={n_iterations_per_design}",
        ]
        if hotspot_res is not None:
            cmd.append(f"ppi.hotspot_res={hotspot_res}")
        print(f"Running job {job_id} backbone design")
        # Using popen is important because
        # we need to pass in the new environment.
        _env = os.environ.copy()
        _env["PYTHONPATH"] = f"{RFANTIBODY_REPO_PATH}/src/rfantibody/rfdiffusion:{_env['PYTHONPATH']}"
        subprocess.Popen(cmd, env = _env).wait()

    def _sequence_design() -> None:
        cmd = [
            "python", f"{RFANTIBODY_REPO_PATH}/scripts/proteinmpnn_interface_design.py",
            f"-pdbdir", f"{out_path}/job_{job_id}/",
            f"-outpdbdir", f"{out_path}/job_{job_id}/mpnn",
            f"-checkpoint_path", f"ProteinMPNN_v48_noise_0.2.pt",
            f"-seqs_per_struct", f"{num_seqs}",
        ]
        print(f"Running job {job_id} sequence design")
        subprocess.check_call(cmd)

    def _structure_prediction() -> None:
        cmd = [
            "python", f"{RFANTIBODY_REPO_PATH}/scripts/rf2_predict.py",
            f"--config-path={RFANTIBODY_REPO_PATH}/src/rfantibody/rf2/config",
            f"input.pdb_dir={out_path}/job_{job_id}/mpnn",
            f"output.pdb_dir={out_path}/job_{job_id}/rf2ab",
            f"model.model_weights=RF2_ab.pt",
        ]
        print(f"Running job {job_id} structure prediction")
        subprocess.check_call(cmd)

    _backbone_design()
    _sequence_design()
    _structure_prediction()

    ret = [
        (str(out_file.relative_to(out_path)), out_file.read_bytes())
        for out_file in out_path.rglob("*")
        if out_file.is_file() and out_file.name not in (
            "input_target_file.pdb",
            "input_framework_file.pdb",
        )
    ]

    if tmp_d_used:
        # rmtree
        for root, dirs, files in out_path.walk(top_down = False):
            for fn in files:
                (root / fn).unlink()
            for dn in dirs:
                (root / dn).rmdir()
    return ret


@app.local_entrypoint()
def main(
    input_target_filepath: str,
    input_framework_filepath: str,
    num_designs: int = 500,
    num_gpus: int = 1,
    num_seqs: int = 8,
    run_name: str | None = None,
    contigs: str | None = None,
    n_iterations_per_design: int = 50,
    hotspot_res: str | None = None,
) -> int:
    print(f"Running {num_gpus} parallel jobs")
    if num_gpus > 1:
        num_designs_per_gpu = math.ceil(num_designs / num_gpus)
    else:
        num_designs_per_gpu = num_designs
    all_data = []
    input_target_data = Path(input_target_filepath).read_bytes()
    input_framework_data = Path(input_framework_filepath).read_bytes()
    if run_name is None:
        run_name = f"{datetime.now().strftime("%Y%m%d%H%M")}_rfantibody_run"
    for data in run.starmap(
        [
            (
                input_target_data,
                input_framework_data,
                i,
                run_name,
                num_designs_per_gpu,
                num_seqs,
                contigs,
                n_iterations_per_design,
                hotspot_res,
            )
            for i in range(num_gpus)
        ],
        order_outputs = False,
    ):
        all_data.extend(data)
    
    if vol is None:
        process_and_merge_data.local(all_data, f"./out/{run_name}", num_designs_per_gpu)
    else:
        process_and_merge_data.remote(all_data, f"{VOLUME_PATH}/{run_name}", num_designs_per_gpu)

    return 0
