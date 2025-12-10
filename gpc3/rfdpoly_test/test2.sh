#This script works with apptainer
APPTAINER_PATH="/paperspace/apps/RFDpoly/exec/SE3nv.sif"
MODEL_WEIGHTS_PATH="/paperspace/apps/RFDpoly/weights/train_session2024-07-08_1720455712_BFF_3.00.pt"
export RFDPOLY_CKPT_PATH=$MODEL_WEIGHTS_PATH
RFDPOLY_DIR="/paperspace/apps"
DESIGN_DIR="/paperspace/apps/RFDpoly/design_jobs"

$APPTAINER_PATH $RFDPOLY_DIR/RFDpoly/rf_diffusion/run_inference.py --config-name=multi_polymer \
diffuser.T=50 \
inference.ckpt_path=$MODEL_WEIGHTS_PATH \
inference.input_pdb=$RFDPOLY_DIR/RFDpoly/rf_diffusion/test_data/DBP035.pdb \
inference.num_designs=3 \
contigmap.contigs=[\'33\ 33\ 75\'] \
contigmap.polymer_chains=[\'dna\',\'rna\',\'protein\'] \
inference.output_prefix=$DESIGN_DIR/test_outputs/basic_uncond_test02
