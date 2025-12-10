source ~/scripts/source_conda.sh
conda activate SE3nv

export CUDA_VISIBLE_DEVICES=0

# Running RFdiffusion
python /home/paperspace/apps/RFdiffusion/scripts/run_inference.py \
inference.output_prefix=C3/C3 inference.input_pdb=cd3e_cl_00.pdb \
'contigmap.contigs=[A33-121/0 50-150]' \
'ppi.hotspot_res=[A36,A50,A53,A78]' \
inference.num_designs=2500 \
denoiser.noise_scale_ca=0.33 \
denoiser.noise_scale_frame=0.33 \
inference.ckpt_override_path=/home/paperspace/apps/RFdiffusion/models/Complex_base_ckpt.pt \
"potentials.guiding_potentials=['type:binder_ROG,weight:1,min_dist:5','type:interface_ncontacts,weight:1']" potentials.guide_scale=2 potentials.guide_decay='quadratic'
