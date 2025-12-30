source ~/apps/source_conda.sh
conda deactivate
conda activate foundry

python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_0.json' out_dir=./out_0 n_batches=5 diffusion_batch_size=4
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_1.json' out_dir=./out_1 n_batches=5 diffusion_batch_size=4
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_3.json' out_dir=./out_3 n_batches=5 diffusion_batch_size=4
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_5.json' out_dir=./out_5 n_batches=5 diffusion_batch_size=4
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_10.json' out_dir=./out_10 n_batches=5 diffusion_batch_size=4
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_15.json' out_dir=./out_15 n_batches=5 diffusion_batch_size=4
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_30.json' out_dir=./out_30 n_batches=5 diffusion_batch_size=4
