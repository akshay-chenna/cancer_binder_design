uvx --python 3.12 --with pandas modal run modal_boltzgen.py \
  --input-yaml hla_mb.yaml \
  --protocol protein-anything \
  --num-designs 50000 \
  --designs-per-gpu 1000
