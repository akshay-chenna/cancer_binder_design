uvx --python 3.12 --with pandas modal run modal_boltzgen.py \
  --input-yaml gpc3_mb.yaml \
  --protocol protein-anything \
  --num-designs 20000 \
  --designs-per-gpu 800
