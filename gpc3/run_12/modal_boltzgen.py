"""
modal run modal_boltzgen.py \
  --input-yaml 1g13prot.yaml \
  --protocol protein-anything \
  --num-designs 8 \
  --designs-per-gpu 2
"""

import os
import re
from pathlib import Path

from modal import App, Image

GPU = os.environ.get("GPU", "B200")
TIMEOUT = int(os.environ.get("TIMEOUT", 18*60))


def download_boltzgen_models():
    """Download all boltzgen models during image build to avoid runtime timeouts."""
    import subprocess

    print("Downloading boltzgen models...")
    subprocess.run(["boltzgen", "download", "all"], check=True)
    print("Model download complete")


image = (
    Image.debian_slim()
    .apt_install("git", "wget", "build-essential")
    .pip_install("torch>=2.4.1", "pandas")
    .run_commands(
        "git clone https://github.com/HannesStark/boltzgen /root/boltzgen",
        "cd /root/boltzgen && git checkout 58c1eed2b07f00fd5263f78fe2821c80d6875699 && pip install -e .",
        gpu="a10g",
    )
    .run_function(download_boltzgen_models, gpu="a10g")
)

app = App("AC_boltzgen", image=image)


# File handling utilities
def should_renumber_file(file_path: str, base_name: str) -> bool:
    """Check if file should be renumbered based on design directories."""
    design_dirs = [
        "intermediate_designs/",
        "intermediate_designs_inverse_folded/",
        "intermediate_designs_inverse_folded/refold_design_cif/",
        "intermediate_designs_inverse_folded/refold_cif/",
        "intermediate_designs_inverse_folded/fold_out_design_npz/",
        "intermediate_designs_inverse_folded/fold_out_npz/",
        "intermediate_designs_inverse_folded/metrics_tmp/",
    ]
    
    if not any(file_path.startswith(d) for d in design_dirs):
        return False
    
    filename = file_path.split("/")[-1]
    pattern = rf'{re.escape(base_name)}_\d+\.(cif|npz)'
    if filename.startswith(("metrics_", "data_")):
        pattern = rf'(metrics|data)_{re.escape(base_name)}_\d+\.npz'
    
    return bool(re.search(pattern, filename))


def should_skip_file(file_path: str) -> bool:
    """Files to keep only from first job."""
    filename = file_path.split("/")[-1]
    return (
        file_path.startswith("config/") or
        "lightning_logs/" in file_path or
        (filename.endswith(".cif") and "_" not in filename) or
        filename == "steps.yaml"
    )


def should_merge_csv(file_path: str) -> bool:
    """CSV files that need merging across jobs."""
    filename = file_path.split("/")[-1]
    return filename in ["per_target_metrics_analyze.csv", "aggregate_metrics_analyze.csv"]


def should_merge_pkl(file_path: str) -> bool:
    """Pickle files that need merging across jobs."""
    return file_path.split("/")[-1] == "ca_coords_sequences.pkl.gz"


def merge_csv_files(csv_files_by_job: dict[int, bytes], base_name: str, designs_per_job: int) -> bytes:
    """Merge CSV files from multiple jobs, renumbering design IDs."""
    import pandas as pd
    from io import StringIO, BytesIO
    
    if not csv_files_by_job:
        return b""
    
    def renumber_field(field_str):
        if pd.isna(field_str):
            return field_str
        match = re.search(rf'{re.escape(base_name)}_(\d+)', str(field_str))
        if match:
            old_num = int(match.group(1))
            new_num = job_num * designs_per_job + old_num
            return field_str.replace(f'{base_name}_{old_num}', f'{base_name}_{new_num}')
        return field_str
    
    dfs = []
    for job_num in sorted(csv_files_by_job.keys()):
        df = pd.read_csv(StringIO(csv_files_by_job[job_num].decode('utf-8')))
        
        for col in ['id', 'file_name']:
            if col in df.columns:
                df[col] = df[col].apply(renumber_field)
        
        dfs.append(df)
    
    combined_df = pd.concat(dfs, ignore_index=True)
    output = BytesIO()
    combined_df.to_csv(output, index=False)
    return output.getvalue()


def merge_pkl_files(pkl_files_by_job: dict[int, bytes], base_name: str, designs_per_job: int) -> bytes:
    """Merge pickle files from multiple jobs, renumbering design keys/IDs."""
    import pickle
    import gzip
    from io import BytesIO
    import pandas as pd
    
    if not pkl_files_by_job:
        return b""
    
    merged_data = None
    
    for job_num in sorted(pkl_files_by_job.keys()):
        with gzip.open(BytesIO(pkl_files_by_job[job_num]), 'rb') as f:
            data = pickle.load(f)
        
        def renumber_id(id_str):
            if pd.isna(id_str):
                return id_str
            match = re.search(rf'{re.escape(base_name)}_(\d+)', str(id_str))
            if match:
                old_num = int(match.group(1))
                new_num = job_num * designs_per_job + old_num
                return id_str.replace(f'{base_name}_{old_num}', f'{base_name}_{new_num}')
            return id_str
        
        if isinstance(data, pd.DataFrame):
            if 'id' in data.columns:
                data['id'] = data['id'].apply(renumber_id)
            merged_data = data if merged_data is None else pd.concat([merged_data, data], ignore_index=True)
        
        elif isinstance(data, dict):
            if merged_data is None:
                merged_data = {}
            for key, value in data.items():
                new_key = renumber_id(key) if isinstance(key, str) else key
                merged_data[new_key] = value
        
        elif isinstance(data, list):
            merged_data = data if merged_data is None else merged_data.extend(data) or merged_data
        
        else:
            if merged_data is None:
                merged_data = data
    
    output = BytesIO()
    with gzip.open(output, 'wb') as f:
        pickle.dump(merged_data, f)
    return output.getvalue()


def process_and_merge_files(all_outputs, base_name, designs_per_job, output_dir):
    """Common logic for merging files from multiple jobs."""
    # Group files by job
    jobs_files = {}
    for file_path, content in all_outputs:
        parts = file_path.split("/")
        if parts[0].startswith("job_"):
            job_num = int(parts[0].split("_")[1])
            jobs_files.setdefault(job_num, []).append((file_path, content))
    
    print(f"Merging {len(jobs_files)} jobs into {output_dir}")
    
    written_skip_files = set()
    csv_files_to_merge = {}
    pkl_files_to_merge = {}
    
    # Process jobs in order
    for job_num in sorted(jobs_files.keys()):
        print(f"  Processing job_{job_num} with {len(jobs_files[job_num])} files")
        files_written = files_skipped = files_renumbered = 0
        
        for file_path, content in jobs_files[job_num]:
            clean_path = "/".join(file_path.split("/")[1:])
            
            # Handle CSV merging
            if should_merge_csv(clean_path):
                csv_files_to_merge.setdefault(clean_path, {})[job_num] = content
                continue
            
            # Handle PKL merging
            if should_merge_pkl(clean_path):
                pkl_files_to_merge.setdefault(clean_path, {})[job_num] = content
                continue
            
            # Handle skip files (write only once)
            if should_skip_file(clean_path):
                if clean_path in written_skip_files:
                    files_skipped += 1
                    continue
                written_skip_files.add(clean_path)
                output_path = output_dir / clean_path
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_bytes(content)
                files_written += 1
                continue
            
            # Renumber design files
            if should_renumber_file(clean_path, base_name):
                clean_path = re.sub(
                    rf'{re.escape(base_name)}_(\d+)\.(cif|npz)',
                    lambda m: f'{base_name}_{job_num * designs_per_job + int(m.group(1))}.{m.group(2)}',
                    clean_path
                )
                clean_path = re.sub(
                    rf'(metrics|data)_{re.escape(base_name)}_(\d+)\.npz',
                    lambda m: f'{m.group(1)}_{base_name}_{job_num * designs_per_job + int(m.group(2))}.npz',
                    clean_path
                )
                files_renumbered += 1
            
            output_path = output_dir / clean_path
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_bytes(content)
            files_written += 1
        
        print(f"    Job {job_num}: wrote {files_written}, renumbered {files_renumbered}, skipped {files_skipped}")
    
    # Merge CSV and PKL files
    for csv_path, csv_by_job in csv_files_to_merge.items():
        merged_csv = merge_csv_files(csv_by_job, base_name, designs_per_job)
        (output_dir / csv_path).parent.mkdir(parents=True, exist_ok=True)
        (output_dir / csv_path).write_bytes(merged_csv)
    
    for pkl_path, pkl_by_job in pkl_files_to_merge.items():
        merged_pkl = merge_pkl_files(pkl_by_job, base_name, designs_per_job)
        (output_dir / pkl_path).parent.mkdir(parents=True, exist_ok=True)
        (output_dir / pkl_path).write_bytes(merged_pkl)
    
    print(f"  Merged {len(csv_files_to_merge)} CSV and {len(pkl_files_to_merge)} PKL files")


@app.function(timeout=TIMEOUT * 60, gpu=GPU)
def boltzgen_run_without_filtering(
    yaml_str: str,
    yaml_name: str,
    additional_files: dict[str, bytes],
    protocol: str,
    designs_per_job: int,
    job_idx: int,
    steps: str | None = None,
    cache: str | None = None,
    devices: int | None = None,
    extra_args: str | None = None,
) -> list:
    """Run BoltzGen pipeline WITHOUT filtering for a batch of designs."""
    from subprocess import run
    from tempfile import TemporaryDirectory

    with TemporaryDirectory() as in_dir, TemporaryDirectory() as out_dir:
        yaml_path = Path(in_dir) / yaml_name
        yaml_path.write_text(yaml_str)

        for rel_path, content in additional_files.items():
            file_path = Path(in_dir) / rel_path
            file_path.parent.mkdir(parents=True, exist_ok=True)
            file_path.write_bytes(content)

        # Determine steps (exclude filtering)
        if steps:
            step_list = [s for s in steps.split() if s != "filtering"]
            steps_arg = " ".join(step_list)
        else:
            steps_arg = "design inverse_folding folding design_folding analysis"

        cmd = [
            "boltzgen", "run", str(yaml_path),
            "--output", out_dir,
            "--protocol", protocol,
            "--num_designs", str(designs_per_job),
            "--steps"
        ] + steps_arg.split()

        if cache:
            cmd.extend(["--cache", cache])
        if devices:
            cmd.extend(["--devices", str(devices)])
        if extra_args:
            cmd.extend(extra_args.split())

        print(f"Job {job_idx}: Running {designs_per_job} designs (without filtering)")
        run(cmd, check=True)

        return [
            (f"job_{job_idx}/{out_file.relative_to(out_dir)}", out_file.read_bytes())
            for out_file in Path(out_dir).rglob("*")
            if out_file.is_file()
        ]


@app.function(timeout=TIMEOUT * 60, gpu=GPU)
def boltzgen_filter(
    combined_files: list[tuple[str, bytes]],
    yaml_str: str,
    yaml_name: str,
    additional_files: dict[str, bytes],
    protocol: str,
    designs_per_job: int,
    cache: str | None = None,
    devices: int | None = None,
    extra_args: str | None = None,
) -> list:
    """Run only the filtering step on combined results."""
    from subprocess import run
    from tempfile import TemporaryDirectory

    with TemporaryDirectory() as in_dir, TemporaryDirectory() as work_dir:
        yaml_path = Path(in_dir) / yaml_name
        yaml_path.write_text(yaml_str)

        for rel_path, content in additional_files.items():
            file_path = Path(in_dir) / rel_path
            file_path.parent.mkdir(parents=True, exist_ok=True)
            file_path.write_bytes(content)

        work_path = Path(work_dir) / "combined_results"
        work_path.mkdir(parents=True, exist_ok=True)
        
        # Extract base name
        base_name = None
        for file_path, _ in combined_files:
            clean_path = "/".join(file_path.split("/")[1:])
            if "intermediate_designs/" in clean_path:
                match = re.search(r'/([^/]+)_\d+\.(cif|npz)$', clean_path)
                if match:
                    base_name = match.group(1)
                    break
        
        if not base_name:
            raise ValueError("Could not detect base name from design files")
        
        # Process and merge files
        process_and_merge_files(combined_files, base_name, designs_per_job, work_path)
        
        design_files = list(work_path.glob("intermediate_designs/*.cif"))
        print(f"Merged {len(design_files)} design files total")

        cmd = [
            "boltzgen", "run", str(yaml_path),
            "--output", str(work_path),
            "--protocol", protocol,
            "--steps", "filtering",
        ]

        if cache:
            cmd.extend(["--cache", cache])
        if devices:
            cmd.extend(["--devices", str(devices)])
        if extra_args:
            cmd.extend(extra_args.split())

        print("Running filtering on combined results")
        run(cmd, check=True)

        return [
            (str(out_file.relative_to(work_path)), out_file.read_bytes())
            for out_file in work_path.rglob("*")
            if out_file.is_file()
        ]


@app.local_entrypoint()
def main(
    input_yaml: str,
    protocol: str = "protein-anything",
    num_designs: int = 100,
    designs_per_gpu: int = 10,
    steps: str | None = None,
    cache: str | None = None,
    devices: int | None = None,
    extra_args: str | None = None,
    out_dir: str = "./out/boltzgen",
    run_name: str | None = None,
):
    """Parallelize BoltzGen by running all steps except filtering in parallel, then filter combined results."""
    import math
    from datetime import datetime

    yaml_path = Path(input_yaml)
    yaml_str = yaml_path.read_text()
    yaml_dir = yaml_path.parent


    # Collect additional files referenced in YAML
    additional_files = {}

    # First pass: Get files directly referenced in the main YAML
    for match in re.finditer(r"(?:path:\s*|^\s*-\s+)([^\s\n]+\.(?:yaml|cif|pdb))", yaml_str, re.MULTILINE):
        ref_file = match.group(1)
        ref_path = yaml_dir / ref_file
        if ref_path.exists():
            additional_files[ref_file] = ref_path.read_bytes()
            print(f"Added: {ref_file}")
        else:
            print(f"Warning: Referenced file not found: {ref_file}")

    # Second pass: Get files referenced in scaffold YAML files
    scaffold_files_to_check = [f for f in additional_files.keys() if f.endswith('.yaml')]
    for scaffold_file in scaffold_files_to_check:
        scaffold_content = additional_files[scaffold_file].decode('utf-8')
        scaffold_dir = yaml_dir / Path(scaffold_file).parent
        
        for match in re.finditer(r"path:\s*([^\s\n]+\.(?:cif|pdb))", scaffold_content, re.MULTILINE):
            nested_ref = match.group(1)
            # Try relative to scaffold file location
            nested_path = scaffold_dir / nested_ref
            
            if nested_path.exists():
                # Store with same directory structure
                nested_key = str(Path(scaffold_file).parent / nested_ref)
                if nested_key not in additional_files:
                    additional_files[nested_key] = nested_path.read_bytes()
                    print(f"Added from {scaffold_file}: {nested_key}")
            else:
                print(f"Warning: File referenced in {scaffold_file} not found: {nested_ref}")

    print(f"\nTotal files to upload: {len(additional_files)}")
    for f in additional_files.keys():
        print(f"  - {f}")

    num_jobs = math.ceil(num_designs / designs_per_gpu)
    
    print(f"Running {num_designs} designs across {num_jobs} parallel GPU jobs")
    print(f"Each job will generate {designs_per_gpu} designs")
    print("Phase 1: Parallel execution (all steps except filtering)")

    # Phase 1: Run all steps except filtering in parallel
    all_outputs = []
    for outputs in boltzgen_run_without_filtering.starmap(
        [
            (
                yaml_str,
                yaml_path.name,
                additional_files,
                protocol,
                min(designs_per_gpu, num_designs - i * designs_per_gpu),
                i,
                steps,
                cache,
                devices,
                extra_args,
            )
            for i in range(num_jobs)
        ],
        order_outputs=False,
    ):
        all_outputs.extend(outputs)
        print(f"Completed job with {len(outputs)} output files")

    # Setup output directories
    today = datetime.now().strftime("%Y%m%d%H%M")[2:]
    out_dir_full = Path(out_dir) / (run_name or today)
    unfiltered_dir = out_dir_full / "unfiltered"
    combined_dir = out_dir_full / "combined"
    filtered_dir = out_dir_full / "filtered"
    
    # Save unfiltered outputs
    print(f"\nSaving unfiltered outputs to: {unfiltered_dir}")
    for out_file, out_content in all_outputs:
        output_path = unfiltered_dir / out_file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(out_content)
    
    # Extract base name
    base_name = None
    for file_path, _ in all_outputs:
        clean_path = "/".join(file_path.split("/")[1:])
        if "intermediate_designs/" in clean_path:
            match = re.search(r'/([^/]+)_\d+\.(cif|npz)$', clean_path)
            if match:
                base_name = match.group(1)
                break
    
    if not base_name:
        raise ValueError("Could not detect base name from design files")
    
    # Save combined results
    print(f"\nPhase 2: Combining and saving to: {combined_dir}")
    process_and_merge_files(all_outputs, base_name, designs_per_gpu, combined_dir)
    
    design_files = list(combined_dir.glob("intermediate_designs/*.cif"))
    print(f"Saved {len(design_files)} combined design files")
    
    # Phase 3: Run filtering
    print("\nPhase 3: Running filtering on combined results")
    filtered_outputs = boltzgen_filter.remote(
        all_outputs,
        yaml_str,
        yaml_path.name,
        additional_files,
        protocol,
        designs_per_gpu,
        cache,
        devices,
        extra_args,
    )

    # Save filtered outputs
    print(f"Saving filtered outputs to: {filtered_dir}")
    for out_file, out_content in filtered_outputs:
        output_path = filtered_dir / out_file
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(out_content)

    print(f"\nResults saved to: {out_dir_full}")
    print(f"  - Unfiltered: {unfiltered_dir} ({len(all_outputs)} files)")
    print(f"  - Combined: {combined_dir} ({len(design_files)} design files)")
    print(f"  - Filtered: {filtered_dir} ({len(filtered_outputs)} files)")
    print(f"All {num_designs} designs processed!")
