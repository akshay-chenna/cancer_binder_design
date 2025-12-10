i=job_72-struct_108_generated.cif

rf3 fold inference_engine=rf3 inputs=$i ckpt_path='rf3' num_steps='200' early_stopping_plddt_threshold='0.2'
