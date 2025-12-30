"""
VOLUME_NAME=gpc3_test

uvx modal volume delete $VOLUME_NAME -y
uvx modal volume create --version=2 $VOLUME_NAME

VOLUME_NAME=$VOLUME_NAME TIMEOUT_MINUTES=1400 GPU=L40S uvx --with numpy modal run modal_rfd1-solmpnn-af2.py \
        --pdb-in reordered_modified_gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb \
        --hotspot-res "B36,B48,B59,B81,B86,B89,B99" \
        --partialt 5 \
        --noise-scale-ca 1 \
        --noise-scale-frame 1 \
        --guiding-potentials "type:binder_ROG,weight:1,min_dist:5" \
        --guide-scale 2 \
        --guide-decay "quadratic" \
        --n-mpnn 2 \
        --mpnn-t 0.1 \
        --fixed-chains "B" \
        --designable-chains "A" \
        --no-use-binder-template \
        --num-designs 4 \
        --num-gpus 2
rm -rf volume_data
mkdir volume_data
uvx modal volume get $VOLUME_NAME / volume_data --force
"""

import os
import math
from pathlib import Path
import itertools
from modal import App, Image, Volume

GPU = os.environ.get("GPU","L40S")
TIMEOUT = int(os.environ.get("TIMEOUT_HOURS","24")) * 60 * 60
VOLUME_NAME = os.environ.get("VOLUME_NAME", "rfd1_volume")
RFDIFFUSION_REPO_PATH = "/root/RFdiffusion"
IPAE_NORMALIZATION_FACTOR = 31
PARAMS_DIR = "/root/.cache/params"

volume = Volume.from_name(VOLUME_NAME, create_if_missing=True, version = 2)

def download_weights():
    import subprocess
    urls = [
        "http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt",
    ]
    processes = [
        subprocess.Popen(
            [
                "aria2c",
                "-x",
                "16",
                url,
            ],
        )
        for url in urls
    ]
    for process in processes:
        process.wait()
    print("RFdiffusion weights downloaded.")


image_rfd = (
    Image.debian_slim(python_version = "3.12")
    .apt_install("git", "aria2")
    .uv_pip_install("torch", "hydra-core", "pyrsistent", "opt-einsum", "e3nn", gpu = "A10G")
    .uv_pip_install("dgl", find_links = "https://data.dgl.ai/wheels/torch-2.4/cu124/repo.html", gpu = "A10G")
    .run_commands(f"git clone https://github.com/RosettaCommons/RFdiffusion {RFDIFFUSION_REPO_PATH}")
    .pip_install("uv")
    .run_commands(
        f"uv pip install --system --compile-bytecode -e {RFDIFFUSION_REPO_PATH}/env/SE3Transformer",
        f"uv pip install --system --compile-bytecode -e {RFDIFFUSION_REPO_PATH}",
    )
    .pip_install(["biotite"])
    .run_function(download_weights)
)

with image_rfd.imports():
    import subprocess
    from pathlib import Path
    import biotite.structure as struc
    from biotite.structure.io import load_structure

app = App("RFD1", image=image_rfd)

@app.function(gpu=GPU, timeout=TIMEOUT, volumes={"/data": volume})
def run_rfd_parallel(
    pdb,
    contig,
    hotspot_res,
    partialt,
    noise_scale_ca,
    noise_scale_frame,
    ckpt_override_path,
    guiding_potentials,
    guide_scale,
    guide_decay,
    job_idx,
    num_designs_per_job,  
):

    pdb_path = Path("pdb.pdb")
    pdb_path.write_text(pdb)

    if contig is None:
        atom_array = load_structure(pdb_path)
        chain_ids = struc.get_chains(atom_array)
        start_stop = []
        for chain in chain_ids:
            chain_mask = (atom_array.chain_id == chain)
            chain_res_ids = atom_array.res_id[chain_mask]
            start_stop.append(str(chain_res_ids[0]))
            start_stop.append(str(chain_res_ids[-1]))
        
        contig = [f"{start_stop[1]}-{start_stop[1]}/0 "]

        for i, item in enumerate(start_stop):
            if i//2 == 0:
                continue
            if i%2 == 0:
                contig.append(chain_ids[i//2])
                contig.append(item)
                contig.append("-")
            else:
                contig.append(item)
                contig.append("/0 ")

        contig = "".join(contig)

    cmd = [
        "python", f"{RFDIFFUSION_REPO_PATH}/scripts/run_inference.py",
        f"inference.output_prefix=/data/rfd/job_{job_idx}",
        f"inference.input_pdb={pdb_path}",
        f"inference.num_designs={num_designs_per_job}",
        f"denoiser.noise_scale_ca={noise_scale_ca}",
        f"denoiser.noise_scale_frame={noise_scale_frame}",
        f"inference.ckpt_override_path={ckpt_override_path}",
        f"contigmap.contigs=[{contig}]",
    ]
    
    if hotspot_res is not None:
        cmd.append(f"ppi.hotspot_res=[{hotspot_res}]")
    
    if partialt is not None:
        cmd.append(f"diffuser.partial_T={partialt}")

    if guiding_potentials is not None:
        cmd.append(f'potentials.guiding_potentials=["{guiding_potentials}"]')
        cmd.append(f"potentials.guide_scale={guide_scale}")
        cmd.append(f"potentials.guide_decay={guide_decay}")

    subprocess.run(cmd, check=True)

@app.function(volumes={"/data": volume})
def list_pdbs(path:str, glob_str: str="*.pdb") -> list[str]:
    full_path = Path("/data") / path.lstrip("/")
    if not full_path.exists():
        raise FileNotFoundError(f"Path {full_path} does not exist on volume")
    return [str(f.relative_to("/data")) for f in full_path.glob(glob_str)]

@app.function(volumes={"/data": volume})
def clean_nan_from_files(path: str, glob_str: str = "*.cif"):
    full_path = Path("/data") / path.lstrip("/")
    if not full_path.exists():
        return
    
    for file_path in full_path.glob(glob_str):
        if file_path.is_file():
            lines = file_path.read_text().splitlines()
            clean_lines = [line for line in lines if "nan" not in line.lower()]
            if len(lines) != len(clean_lines):
                file_path.write_text("\n".join(clean_lines) + "\n")
                print(f"Cleaned 'nan' from {file_path.name}")

image_mpnn = (
    Image.debian_slim(python_version="3.12")
    .apt_install("git", "wget", "pip")
    .env({"CCD_MIRROR_PATH": "", "PDB_MIRROR_PATH": ""})
    .run_commands("pip install -q 'rc-foundry[all]'")
    .run_commands("foundry install rfd3 proteinmpnn solublempnn rf3")
)

with image_mpnn.imports():
    from pathlib import Path

    from biotite.structure.io import load_structure
    from mpnn.inference_engines.mpnn import MPNNInferenceEngine
    
    import gc
    import torch

app_mpnn = App("Solube MPNN", image=image_mpnn)

@app_mpnn.function(gpu="A10G", timeout=TIMEOUT, volumes={"/data": volume})
def run_mpnn_parallel(
    batch,
    n_mpnn,
    mpnn_t,
    fixed_chains,
    designable_chains,
):
    # Explicitly create output directory on the volume
    out_dir = Path("/data/mpnn")
    out_dir.mkdir(parents=True, exist_ok=True)

    for pdb_path in batch:
        pdb_path = Path("/data") / pdb_path
        atom_array = load_structure(pdb_path)
        
        engine_config = {
            'model_type': "protein_mpnn",
            "is_legacy_weights": True,
            "checkpoint_path": "/root/.foundry/checkpoints/solublempnn_v_48_020.pt",
            "out_directory": "/data/mpnn",
            "write_structures": True,
            "write_fasta": True,
        }

        input_configs = [
            {
                "batch_size": n_mpnn,
                "remove_waters": True,
                "fixed_chains": fixed_chains,
                "temperature": mpnn_t,
                "name": f"mpnn_{pdb_path.stem}",
            }
        ]

        model = MPNNInferenceEngine(**engine_config)
        mpnn_outputs = model.run(input_dicts=input_configs, atom_arrays=[atom_array])

        gc.collect()
        torch.cuda.empty_cache()

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

image_af2 = (
    Image.debian_slim(python_version = "3.12")
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

app_af2 = App("af2-binder-eval", image = image_af2)
with image_af2.imports():
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


@app_af2.function(gpu = GPU, timeout = TIMEOUT, volumes = {"/data": volume})
def run_af2_parallel(
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

    out_dir = Path("/data") / out_rel_dir
    out_dir.mkdir(parents = True, exist_ok = True)
    pae_dir = out_dir / "paes"
    pae_dir.mkdir(parents = True, exist_ok = True)

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
            structure_path = Path("/data") / structure_rel_path
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

            pae_file = f"{structure_path.stem}.npy"
            pae_file_path = pae_dir / pae_file
            np.save(str(pae_file_path), af2_binder_model.aux["pae"])

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


@app_af2.function(volumes = {"/data": volume}, timeout = 60 * 60 * 2)
def save_metrics(results: list[dict], out_rel_dir: Path) -> None:
    metrics_save_filepath = Path("/data").joinpath(out_rel_dir, "metrics.delta")
    pl.LazyFrame(results).collect().write_delta(metrics_save_filepath, mode = "append")


@app.local_entrypoint()
def main(
    pdb_in: str,
    contig: str | None = None,
    hotspot_res: str | None = None,
    partialt: int | None = None, # Partial T
    noise_scale_ca: float = 0.5,
    noise_scale_frame: float = 0.5,
    ckpt_override_path: str = "Complex_base_ckpt.pt",
    guiding_potentials: str | None = None,
    guide_scale: float | None = None,
    guide_decay: str | None = None,
    n_mpnn: int = 4,
    mpnn_t: float = 0.1,
    fixed_chains: str = "B",
    designable_chains: str = "A",
    use_multimer: bool = True,
    use_initial_guess: bool = True,
    use_templates: bool = True,
    use_binder_template: bool = False,
    rm_template_ic: bool = True,
    num_recycles: int = 5,
    model_nums = [0, 1],
    num_designs: int = 10,
    num_gpus: int = 1,
) -> None:

    pdb_path = Path(pdb_in)
    pdb_str = pdb_path.read_text()
    
    num_jobs = num_gpus
    num_designs_per_job = num_designs // num_jobs

    binder_chain = designable_chains
    target_chain = fixed_chains

    # Execute RFD1
    for outputs in run_rfd_parallel.starmap(
        [
            (
                pdb_str,
                contig,
                hotspot_res,
                partialt,
                noise_scale_ca,
                noise_scale_frame,
                ckpt_override_path,
                guiding_potentials,
                guide_scale,
                guide_decay,
                job_idx,
                num_designs_per_job,
            )
            for job_idx in range(num_jobs)
        ],
        order_outputs = False,
    ):
        print("Completed RFD1")

    # Execute MPNN
    structure_paths = list_pdbs.remote(path="rfd", glob_str="*.pdb")
    fixed_chains = fixed_chains.split(",")
    designable_chains = designable_chains.split(",")

    with app_mpnn.run():
        for mpnn_outputs in run_mpnn_parallel.starmap(
            [
                (
                    batch,
                    n_mpnn,
                    mpnn_t,
                    fixed_chains,
                    designable_chains,
                )
                for batch in itertools.batched(structure_paths,num_designs_per_job)
            ],
            order_outputs = False,
        ):
            print("Completed MPNN")

    # Clean MPNN outputs (remove 'nan' lines)
    clean_nan_from_files.remote(path="mpnn", glob_str="*.cif")
    
    # Execute AF2
    structure_paths = list_pdbs.remote(path="mpnn", glob_str="*.cif")
    num_files_to_process = len(structure_paths)
    batch_size = math.ceil(num_files_to_process / num_gpus)

    out_rel_dir = "af2"
    
    with app_af2.run():
        results = []
        for af2_outputs in run_af2_parallel.starmap(
            [
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
                for batch in itertools.batched(structure_paths,batch_size)
            ],
            order_outputs = False,
        ):
            results.extend(af2_outputs)
            print("Completed AF2")  
        
    with app_af2.run():
        save_metrics.remote(results, out_rel_dir)

    print("Completed all jobs")
