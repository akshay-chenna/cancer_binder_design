source /paperspace/apps/source_conda.sh
conda activate SE3nv

export CUDA_VISIBLE_DEVICES=5

# Running RFdiffusion
python /paperspace/apps/RFdiffusion/scripts/run_inference.py \
inference.output_prefix=C2/C2 inference.input_pdb=af3_cd8ab_ecto_nosignalp_model_0_tidy.pdb \
'contigmap.contigs=[A1-111/0 B1-114/0 60-150]' \
inference.num_designs=2500 \
denoiser.noise_scale_ca=0.33 \
denoiser.noise_scale_frame=0.33 \
inference.ckpt_override_path=/paperspace/apps/RFdiffusion/models/Complex_base_ckpt.pt \
"potentials.guiding_potentials=['type:binder_ROG,weight:1,min_dist:5','type:interface_ncontacts,weight:1']" potentials.guide_scale=2 potentials.guide_decay='quadratic'
