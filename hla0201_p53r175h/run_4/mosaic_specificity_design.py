import mosaic
import jax
from mosaic.optimizers import (
    simplex_APGM,
    gradient_MCMC,
)

import mosaic.losses.structure_prediction as sp
from mosaic.models.boltz1 import Boltz1

from mosaic.common import TOKENS
from mosaic.losses.transformations import SoftClip
from mosaic.notebook_utils import pdb_viewer
from jaxtyping import Float, Array
from mosaic.common import LossTerm
from mosaic.structure_prediction import TargetChain
from mosaic.models.af2 import AlphaFold2
from mosaic.proteinmpnn.mpnn import ProteinMPNN

boltz1 = Boltz1()

target_sequence = "PEPTIDESEQ"
