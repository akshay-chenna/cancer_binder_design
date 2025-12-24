source ~/apps/source_conda.sh
conda deactivate
conda activate foundry

#python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='demo.json' out_dir=./outputs n_batches=33 diffusion_batch_size=3
#python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='demo_2.json' out_dir=./outputs_2 n_batches=33 diffusion_batch_size=3
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='demo_new.json' out_dir=./outputs_new n_batches=33 diffusion_batch_size=3
