GPU=H100 TIMEOUT=1440 uvx --python 3.12 --with pandas modal run modal_boltzgen.py \
  --input-yaml pdl1_ab.yaml \
  --protocol nanobody-anything \
  --num-designs 20000 \
  --designs-per-gpu 200
