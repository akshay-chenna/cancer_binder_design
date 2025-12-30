source ~/apps/source_conda.sh
conda deactivate
conda activate foundry

#python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='pd_01.json' out_dir=./out_01 n_batches=5 diffusion_batch_size=4
#python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='pd_05.json' out_dir=./out_05 n_batches=5 diffusion_batch_size=4
#python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='pd_1.json' out_dir=./out_1 n_batches=5 diffusion_batch_size=4
#python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='pd_3.json' out_dir=./out_3 n_batches=5 diffusion_batch_size=4

rfd3 ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='pd_01.json' out_dir=./cli_out_01 n_batches=5 diffusion_batch_size=4
rfd3 ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='pd_05.json' out_dir=./cli_out_05 n_batches=5 diffusion_batch_size=4
rfd3 ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='pd_1.json' out_dir=./cli_out_1 n_batches=5 diffusion_batch_size=4
rfd3 ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='pd_3.json' out_dir=./cli_out_3 n_batches=5 diffusion_batch_size=4
