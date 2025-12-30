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
    from tempfile import NamedTemporaryFile


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

@app.function(gpu = GPU, timeout = TIMEOUT_MINUTES * 60, volumes = {VOLUME_MOUNT: vol})
def worker(
    structure_rel_paths: list[str],
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

    def calc_min_ipae(pae:np.ndarray, L: int) -> float:
        min_target_binder_pae = np.min(pae[: L, L:])
        min_binder_target_pae = np.min(pae[L:, : L])
        return np.min([min_target_binder_pae, min_binder_target_pae])

    af2_binder_model = mk_afdesign_model(
        protocol = "binder",
        data_dir = PARAMS_DIR,
        use_multimer = use_multimer,
        use_initial_guess = use_initial_guess,
        use_templates = use_templates,
    )

    results = []
    for structure_rel_path in structure_rel_paths:
        structure_path = Path(VOLUME_MOUNT) / structure_rel_path
        filepath = structure_path
        if "cif" in filepath.suffix:
            with NamedTemporaryFile(suffix = ".pdb", delete = False) as tmp_file:
                s = gemmi.read_structure(str(filepath))
                s.setup_entities()
                s.write_pdb(tmp_file.name)
            filepath = Path(tmp_file.name)
        
        af2_binder_model.prep_inputs(
            pdb_filename = str(filepath),
            chain = target_chain,
            binder_chain = binder_chain,
            use_binder_template = use_binder_template,
            rm_template_ic = rm_template_ic,
        )
        
        binder_seq = get_binder_seq(filepath, binder_chain)
        
        af2_binder_model.predict(
            binder_seq,
            num_recycles = num_recycles,
            model_nums = model_nums,
            verbose = False,
        )
        
        pred_pdb_name = f"{structure_path.stem}.pdb"
        pred_pdb_path = out_dir / pred_pdb_name
        af2_binder_model.save_pdb(str(pred_pdb_path))

        result = {
            "name": structure_path.stem,
            "plddt": float(np.mean(af2_binder_model.aux["log"]["plddt"])),
            "pae": float(np.mean(af2_binder_model.aux["log"]["pae"])),
            "ptm": float(af2_binder_model.aux["log"]["ptm"]),
            "rmsd": float(af2_binder_model.aux["log"]["rmsd"]),
            "ipae": float(np.mean(af2_binder_model.aux["log"]["i_pae"])) * IPAE_NORMALIZATION_FACTOR,
            "min_ipae": calc_min_ipae(af2_binder_model.aux["pae"], af2_binder_model._target_len),
            "iptm": float(af2_binder_model.aux["log"]["i_ptm"]),
        }
        results.append(result)
        af2_binder_model.restart()
    print(f"Finished batch of {len(structure_rel_paths)} file(s).")
    return results


@app.function(volumes = {VOLUME_MOUNT: vol}, timeout = 60 * 60 * 2)
def save_metrics(results: list[dict], out_rel_dir: Path) -> None:
    metrics_save_filepath = Path(VOLUME_MOUNT).joinpath(out_rel_dir, "metrics.delta")
    pl.LazyFrame(results).collect().write_delta(metrics_save_filepath, mode = "append")


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
    structure_paths = list_pdbs.remote(path, glob_str)
    num_files_to_process = len(structure_paths)
    
    if not structure_paths:
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
        for batch in itertools.batched(structure_paths, batch_size)
    ]

    results = []
    for r in worker.starmap(batched_tasks, order_outputs = False):
        results.extend(r)
    save_metrics.remote(results, out_rel_dir)

    return 0
