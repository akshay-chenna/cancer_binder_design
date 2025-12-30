source ~/apps/source_conda.sh
conda deactivate
conda activate foundry

python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='demo_15.json' out_dir=./outputs_15 n_batches=5 diffusion_batch_size=4 dump_trajectories=True prevalidate_inputs=True
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='demo_1.json' out_dir=./outputs_1 n_batches=5 diffusion_batch_size=4 dump_trajectories=True prevalidate_inputs=True
