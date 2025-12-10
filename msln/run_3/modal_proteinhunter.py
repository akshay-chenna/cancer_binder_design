from dataclasses import asdict, dataclass

import math

import modal

import os

from pathlib import Path


VOLUME_NAME = os.environ.get("VOLUME_NAME", "modal-protein-hunter")
VOLUME_PATH = "/results"
GPU = os.environ.get("GPU", "H100")
TIMEOUT_MINUTES = int(os.environ.get("TIMEOUT_MINUTES", 600))
REPO_PATH = "/root/protein_hunter"


vol = modal.Volume.from_name(VOLUME_NAME, create_if_missing = True)


def download_boltz2_data() -> None:
    from boltz.main import download_boltz2
    p = Path.home().joinpath(".boltz")
    p.mkdir(parents = True, exist_ok = True)
    download_boltz2(p)


image: modal.Image = (
    modal.Image.debian_slim(python_version = "3.10")
    .apt_install("git", "wget")
    .run_commands(
        f"git clone https://github.com/yehlincho/Protein-Hunter {REPO_PATH}",
    )
    .uv_pip_install(
        "matplotlib",
        "ml-collections",
        "tqdm",
        "requests",
        "logmd==0.1.45",
        "py2Dmol",
        "py3Dmol",
        "pypdb",
        "PyYAML",
        "prody",
        "seaborn",
        f"{REPO_PATH}/boltz_ph",
        gpu = "A10G",
    )
    .run_commands(
        f"bash {REPO_PATH}/LigandMPNN/get_model_params.sh {REPO_PATH}/LigandMPNN/model_params",
    )
    .run_function(download_boltz2_data)
)


with image.imports():
    from datetime import datetime
    import shutil
    import subprocess
    import tempfile


app = modal.App("Protein-Hunter", image = image)


@dataclass
class ProteinHunterConfig:
    num_designs: int
    num_cycles: int
    protein_seqs: str
    protein_ids: str
    protein_msas: str = ""
    gpu_id: int = 0
    name: str | None = None
    min_design_protein_length: int = 90
    max_design_protein_length: int = 150
    high_iptm_threshold: float = 0.7
    template_path: str = ""
    name: str | None = None
    cyclic: bool = False
    mode: str = "binder"
    binder_chain: str = "A"
    ligand_id: str = "B"
    ligand_smiles: str = ""
    ligand_ccd: str = ""
    nucleic_type: str = "dna"
    nucleic_id: str = "B"
    nucleic_seq: str = ""
    template_path: str = ""
    template_chain_id: str = ""
    template_cif_chain_id: str = ""
    no_potentials: bool = True
    diffuse_steps: int = 200
    recycling_steps: int = 3
    boltz_model_version: str = "boltz2"
    boltz_model_path: str = "~/.boltz/boltz2_conf.ckpt"
    ccd_path: str = "~/.boltz/mols"
    randomly_kill_helix_feature: bool = False
    negative_helix_constant: float = 0.2
    logmd: bool = False
    save_dir: str = ""
    add_constraints: bool = False
    contact_residues: str = ""
    omit_AA: str = "C"
    exclude_P: bool = False
    percent_X: int = 80
    plot: bool = False
    constraint_target_chain: str = ""
    no_contact_filter: bool = False
    max_contact_filter_retries: int = 6
    contact_cutoff: float = 15.0
    work_dir: str = ""
    temperature_start: float = 0.05
    temperature_end: float = 0.001
    alanine_bias_start: float = -0.5
    alanine_bias_end: float = -0.2
    alanine_bias: float = False
    high_iptm_threshold: float = 0.7


@app.function(timeout = TIMEOUT_MINUTES * 60, volumes = {VOLUME_PATH: vol} if vol is not None else {})
def process_and_merge_data(all_data: list[tuple[str, bytes]], out_path: str, designs_per_job: int) -> None:
    out_path = Path(out_path)
    out_path.mkdir(exist_ok = True, parents = True)
    for p, data in all_data:
        local_filepath = out_path / p
        local_filepath.parent.mkdir(parents = True, exist_ok = True)
        local_filepath.write_bytes(data)


@app.function(
    timeout = TIMEOUT_MINUTES * 60,
    volumes = {VOLUME_PATH: vol} if vol is not None else {},
    gpu = GPU,
)
def run(
    config: ProteinHunterConfig,
    job_id: int,
    template_data: list[tuple[str, bytes]],
) -> list[tuple[str, bytes]]:
    tmp_d_used = False
    if vol is None:
        print("No volume specified. Using temporary directory.")
        out_path = Path(tempfile.mkdtemp()) / f"{config.name}"
        tmp_d_used = True
    else:
        out_path = Path(VOLUME_PATH) / f"{config.name}"
    if not out_path.exists():
        out_path.mkdir(parents = True, exist_ok = True)

    if len(template_data) != 0:
        for tn, td in template_data:
            out_path.joinpath(tn).write_bytes(td)
        config.template_path = ",".join(
            str(out_path.joinpath(tn)) for tn, _ in template_data
        )
    config.save_dir = out_path / f"job_{job_id}"

    _flags = {
        "cyclic",
        "randomly_kill_helix_feature",
        "logmd",
        "add_constraints",
        "exclude_P",
        "plot",
        "no_contact_filter",
        "alanine_bias",
    }
    cmd_args = (
        f"--{k}={v}" if k not in _flags else (f"--{k}" if v else None)
        for k, v in asdict(config).items()
    )
    cmd = [
        "uv",
        "run",
        "--directory", f"{REPO_PATH}",
        f"{REPO_PATH}/boltz_ph/design.py",
        *(
            arg for arg in cmd_args if arg is not None
        )
    ]
    print(f"Running job {job_id}")
    subprocess.check_call(cmd)

    ret = [
        (str(out_file.relative_to(out_path)), out_file.read_bytes())
        for out_file in out_path.rglob("*")
        if out_file.is_file()
    ]
    if tmp_d_used:
        shutil.rmtree(out_path)
    return ret


@app.local_entrypoint()
def main(
    num_designs: int,
    num_cycles: int,
    protein_seqs: str,
    protein_ids: str,
    num_gpus: int,
    protein_msas: str = "",
    gpu_id: int = 0,
    name: str | None = None,
    min_design_protein_length: int = 90,
    max_design_protein_length: int = 150,
    high_iptm_threshold: float = 0.7,
    cyclic: bool = False,
    mode: str = "binder",
    binder_chain: str = "A",
    ligand_id: str = "B",
    ligand_smiles: str = "",
    ligand_ccd: str = "",
    nucleic_type: str = "dna",
    nucleic_id: str = "B",
    nucleic_seq: str = "",
    template_path: str = "",
    template_chain_id: str = "",
    template_cif_chain_id: str = "",
    use_potentials: bool = False,
    diffuse_steps: int = 200,
    recycling_steps: int = 3,
    boltz_model_version: str = "boltz2",
    boltz_model_path: str = "~/.boltz/boltz2_conf.ckpt",
    ccd_path: str = "~/.boltz/mols",
    randomly_kill_helix_feature: bool = False,
    negative_helix_constant: float = 0.2,
    logmd: bool = False,
    save_dir: str = "",
    add_constraints: bool = False,
    contact_residues: str = "",
    omit_aa: str = "C",
    exclude_p: bool = False,
    percent_x: int = 80,
    plot: bool = False,
    constraint_target_chain: str = "",
    no_contact_filter: bool = False,
    max_contact_filter_retries: int = 6,
    contact_cutoff: float = 15.0,
    work_dir: str = "",
    temperature_start: float = 0.05,
    temperature_end: float = 0.001,
    alanine_bias_start: float = -0.5,
    alanine_bias_end: float = -0.2,
    alanine_bias: float = False,
) -> int:
    print(f"Running {num_gpus} parallel jobs")
    if num_gpus > 1:
        num_designs_per_gpu = math.ceil(num_designs / num_gpus)
    else:
        num_designs_per_gpu = num_designs

    run_name = f"{datetime.now().strftime('%Y%m%d%H%M')}_proteinhunter_run" if name is None else name

    template_data = []
    template_data_set = set()
    if template_path:
        delim = ","
        for tp in template_path.split(delim):
            tpp = Path(tp)
            if tpp.name in template_data_set:
                continue
            template_data.append(
                (
                    tpp.name,
                    tpp.read_bytes(),
                )
            )
            template_data_set.add(tpp.name)

    protein_hunter_config = ProteinHunterConfig(
        num_designs = num_designs,
        num_cycles = num_cycles,
        protein_seqs = protein_seqs,
        protein_ids = protein_ids,
        protein_msas = protein_msas,
        gpu_id = gpu_id,
        name = run_name,
        percent_X = percent_x,
        min_design_protein_length = min_design_protein_length,
        max_design_protein_length = max_design_protein_length,
        cyclic = cyclic,
        mode = mode,
        binder_chain = binder_chain,
        ligand_id = ligand_id,
        ligand_smiles = ligand_smiles,
        ligand_ccd = ligand_ccd,
        nucleic_type = nucleic_type,
        nucleic_id = nucleic_id,
        nucleic_seq = nucleic_seq,
        template_path = template_path,
        template_chain_id = template_chain_id,
        template_cif_chain_id = template_cif_chain_id,
        no_potentials = not use_potentials,
        diffuse_steps = diffuse_steps,
        recycling_steps = recycling_steps,
        boltz_model_version = boltz_model_version,
        boltz_model_path = boltz_model_path,
        ccd_path = ccd_path,
        randomly_kill_helix_feature = randomly_kill_helix_feature,
        negative_helix_constant = negative_helix_constant,
        logmd = logmd,
        save_dir = save_dir,
        add_constraints = add_constraints,
        contact_residues = contact_residues,
        omit_AA = omit_aa,
        exclude_P = exclude_p,
        plot = plot,
        constraint_target_chain = constraint_target_chain,
        no_contact_filter = no_contact_filter,
        max_contact_filter_retries = max_contact_filter_retries,
        contact_cutoff = contact_cutoff,
        work_dir = work_dir,
        temperature_start = temperature_start,
        temperature_end = temperature_end,
        alanine_bias_start = alanine_bias_start,
        alanine_bias_end = alanine_bias_end,
        alanine_bias = alanine_bias,
        high_iptm_threshold = high_iptm_threshold,
    )

    task_args = []
    for job_id in range(num_gpus):
        protein_hunter_config.num_designs = num_designs_per_gpu
        task_args.append(
            (
                protein_hunter_config,
                job_id,
                template_data,
            )
        )

    output = []
    for data in run.starmap(
        task_args,
        order_outputs = False,
    ):
        output.extend(data)

    # if vol is None:
    #     process_and_merge_data.local(output, f"./out/{run_name}", num_designs_per_gpu)
    # else:
    #     process_and_merge_data.remote(output, f"{VOLUME_PATH}/{run_name}", num_designs_per_gpu)

    return 0
