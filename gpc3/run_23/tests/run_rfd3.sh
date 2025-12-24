source ~/apps/source_conda.sh
conda deactivate
conda activate foundry

python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_1.json' out_dir=./outputs_1 n_batches=25 diffusion_batch_size=4
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_2.json' out_dir=./outputs_2 n_batches=25 diffusion_batch_size=4
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_3.json' out_dir=./outputs_3 n_batches=25 diffusion_batch_size=4
