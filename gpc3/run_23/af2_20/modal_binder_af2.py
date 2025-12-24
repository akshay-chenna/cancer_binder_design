import itertools

import math

import modal

import os

from pathlib import Path


VOLUME_NAME = os.environ.get("VOLUME_NAME", "af2-binder-eval")
GPU = os.environ.get("GPU", "A100")
TIMEOUT_MINUTES = int(os.environ.get("TIMEOUT_MINUTES", 60))
PARAMS_DIR = "/root/.cache/params"
VOLUME_MOUNT = "/vol"
IPAE_NORMALIZATION_FACTOR = 31

vol = modal.Volume.from_name(VOLUME_NAME, create_if_missing = False)


def download_af2_params():
    import subprocess
    path = Path(PARAMS_DIR)
    path.mkdir(parents = True, exist_ok = True)
    url = "https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"
    print(f"Downloading AF2 parameters to {path}...")
    subprocess.check_call(["aria2c", "-x", "8", url, "-d", str(path)])
    subprocess.check_call(["tar", "xvf", str(path / "alphafold_params_2022-12-06.tar"), "-C", str(path)])
    (path / "alphafold_params_2022-12-06.tar").unlink(missing_ok = True)
    print("AF2 parameters downloaded and extracted.")


image = (
    modal.Image.debian_slim(python_version = "3.12")
    .apt_install("git", "aria2", "wget", "build-essential")
    .uv_pip_install(
        "deltalake",
        "flax",
        "gemmi",
        "jax[cuda12]<0.6.0",
        "numpy",
        "polars",
        "rich",
        "colabdesign @ git+https://github.com/sokrypton/ColabDesign.git",
        gpu = GPU,
    )
    .run_function(download_af2_params)
)


app = modal.App("af2-binder-eval", image = image)


with image.imports():
    from colabdesign import mk_afdesign_model
    from colabdesign.shared.utils import clear_mem
    import gemmi
    import numpy as np
    import polars as pl


def get_binder_seq(structure_path: Path, chain: str) -> str:
    structure = gemmi.read_structure(str(structure_path))
    if "," in chain:
        return "".join(
            gemmi.one_letter_code([res.name for res in structure[0][c]])
            for c in chain.split(",")
        )
    return gemmi.one_letter_code([res.name for res in structure[0][chain]])


@app.function(volumes = {VOLUME_MOUNT: vol})
def list_pdbs(path: str, glob_str: str = "*.pdb") -> list[str]:
    full_path = Path(VOLUME_MOUNT) / path.lstrip("/")
    if not full_path.exists():
        raise FileNotFoundError(f"Path {full_path} does not exist on volume.")
    return [str(f.relative_to(VOLUME_MOUNT)) for f in full_path.glob(glob_str)]


def calc_min_ipae(model) -> float:
    pae = model.aux["pae"]
    target_len = model._target_len
    min_target_binder_pae = pae[: target_len, target_len :].min(-1).mean()
    min_binder_target_pae = pae[target_len :, : target_len].min(-1).mean()
    return np.mean([min_target_binder_pae, min_binder_target_pae]).item()


@app.function(gpu = GPU, timeout = TIMEOUT_MINUTES * 60, volumes = {VOLUME_MOUNT: vol})
def worker(
    pdb_rel_paths: list[str],
    binder_chain: str,
    target_chain: str,
    out_rel_dir: str,
    use_multimer: bool = True,
    use_initial_guess: bool = True,
    use_templates: bool = True,
    use_binder_template: bool = False,
    rm_template_ic: bool = True,
    num_recycles: int = 5,
    model_nums: list[int] = [0, 1],
) -> list[dict]:
    clear_mem()
    
    out_dir = Path(VOLUME_MOUNT) / out_rel_dir
    out_dir.mkdir(parents = True, exist_ok = True)

    af2_binder_model = mk_afdesign_model(
        protocol = "binder",
        data_dir = PARAMS_DIR,
        use_multimer = use_multimer,
        use_initial_guess = use_initial_guess,
        use_templates = use_templates,
    )
    
    results = []
    for pdb_rel_path in pdb_rel_paths:
        pdb_path = Path(VOLUME_MOUNT) / pdb_rel_path
        
        af2_binder_model.prep_inputs(
            pdb_filename = str(pdb_path),
            chain = target_chain,
            binder_chain = binder_chain,
            use_binder_template = use_binder_template,
            rm_template_ic = rm_template_ic,
        )
        
        binder_seq = get_binder_seq(pdb_path, binder_chain)
        
        af2_binder_model.predict(
            binder_seq,
            num_recycles = num_recycles,
            model_nums = model_nums,
            verbose = False,
        )
        
        pred_pdb_name = f"{pdb_path.stem}.pdb"
        pred_pdb_path = out_dir / pred_pdb_name
        af2_binder_model.save_pdb(str(pred_pdb_path))

        result = {
            "name": pdb_path.stem,
            "plddt": float(np.mean(af2_binder_model.aux["log"]["plddt"])),
            "pae": float(np.mean(af2_binder_model.aux["log"]["pae"])),
            "ptm": float(af2_binder_model.aux["log"]["ptm"]),
            "rmsd": float(af2_binder_model.aux["log"]["rmsd"]),
            "ipae": float(np.mean(af2_binder_model.aux["log"]["i_pae"])) * IPAE_NORMALIZATION_FACTOR,
            "min_ipae": calc_min_ipae(af2_binder_model),
            "iptm": float(af2_binder_model.aux["log"]["i_ptm"]),
        }
        results.append(result)
        af2_binder_model.restart()
    print(f"Finished batch of {len(pdb_rel_paths)} file(s).")
    metrics_save_filepath = Path(VOLUME_MOUNT).joinpath(out_rel_dir, "metrics.delta")
    pl.LazyFrame(results).collect().write_delta(metrics_save_filepath, mode = "append")
    vol.commit()


@app.local_entrypoint()
def main(
    path: str,
    out_path: str = None,
    glob_str: str = "*.pdb",
    binder_chain: str = "B",
    target_chain: str = "A",
    num_gpus: int = 1,
    use_multimer: bool = True,
    use_initial_guess: bool = True,
    use_templates: bool = True,
    use_binder_template: bool = False,
    rm_template_ic: bool = True,
    num_recycles: int = 5,
    model_nums = [0, 1],
) -> None:
    pdb_paths = list_pdbs.remote(path, glob_str)
    num_files_to_process = len(pdb_paths)
    
    if not pdb_paths:
        print(f"No PDB files found in {path} on volume '{VOLUME_NAME}'.")
        return

    if out_path is None:
        out_rel_dir = str(Path(path) / "af2_binder_eval")
    else:
        out_rel_dir = out_path

    batch_size = math.ceil(num_files_to_process / num_gpus)

    print(f"Found {num_files_to_process} files. Processing {batch_size} file(s) on each of {num_gpus} GPU(s)...")
    
    batched_tasks = [
        (
            batch,
            binder_chain,
            target_chain,
            out_rel_dir,
            use_multimer,
            use_initial_guess,
            use_templates,
            use_binder_template,
            rm_template_ic,
            num_recycles,
            model_nums,
        )
        for batch in itertools.batched(pdb_paths, batch_size)
    ]

    list(worker.starmap(batched_tasks, order_outputs = False))

    return 0