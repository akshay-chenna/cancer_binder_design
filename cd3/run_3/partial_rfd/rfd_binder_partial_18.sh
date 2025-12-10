source /paperspace/apps/source_conda.sh
conda activate SE3nv

export CUDA_VISIBLE_DEVICES=2

file=top_18
length=71

# Running RFdiffusion
python /paperspace/apps/RFdiffusion/scripts/run_inference.py \
inference.output_prefix=${file}/${file} inference.input_pdb=${file}.pdb \
"contigmap.contigs=[${length}-${length}/0 B33-121]" \
'ppi.hotspot_res=[B36,B47,B50,B53,B78,B82]' \
inference.num_designs=1000 \
diffuser.partial_T=25 \
denoiser.noise_scale_ca=0.5 \
denoiser.noise_scale_frame=0.5 \
inference.ckpt_override_path=/paperspace/apps/RFdiffusion/models/Complex_base_ckpt.pt \
potentials.guiding_potentials=[\"type:binder_ROG,weight:2,min_dist:5\"] potentials.guide_scale=2 potentials.guide_decay='quadratic'
