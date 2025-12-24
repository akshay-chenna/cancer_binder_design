source ~/apps/source_conda.sh
conda deactivate
conda activate foundry

#python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion.json' out_dir=./rfd3_outputs n_batches=200 diffusion_batch_size=3
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_1.json' out_dir=./rfd3_outputs_1 n_batches=33 diffusion_batch_size=3
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_2.json' out_dir=./rfd3_outputs_2 n_batches=33 diffusion_batch_size=3
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_3.json' out_dir=./rfd3_outputs_3 n_batches=33 diffusion_batch_size=3
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_4.json' out_dir=./rfd3_outputs_4 n_batches=33 diffusion_batch_size=3
python ~/apps/foundry/models/rfd3/src/rfd3/run_inference.py ckpt_path='/home/ubuntu/.foundry/checkpoints/rfd3_latest.ckpt' inputs='partial_diffusion_5.json' out_dir=./rfd3_outputs_5 n_batches=33 diffusion_batch_size=3
