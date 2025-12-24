source /paperspace/apps/source_conda.sh
conda deactivate
conda activate SE3nv

export CUDA_VISIBLE_DEVICES=1

file=mb_1
length=108

# Running RFdiffusion
python /paperspace/apps/RFdiffusion/scripts/run_inference.py \
inference.output_prefix=${file}/${file} inference.input_pdb=reordered_modified_gpc3_out_1_l108_s599292_mpnn3_model2_relaxed.pdb \
"contigmap.contigs=[${length}-${length}/0 B1-291]" \
'ppi.hotspot_res=[B36,B48,B59,B81,B86,B89,B99]' \
inference.num_designs=125 \
diffuser.partial_T=20 \
denoiser.noise_scale_ca=1 \
denoiser.noise_scale_frame=1 \
inference.ckpt_override_path=/paperspace/apps/RFdiffusion/models/Complex_base_ckpt.pt \
potentials.guiding_potentials=[\"type:binder_ROG,weight:1,min_dist:5\"] potentials.guide_scale=2 potentials.guide_decay='quadratic'
