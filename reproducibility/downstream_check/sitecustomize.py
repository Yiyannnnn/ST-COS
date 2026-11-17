"""Initialize method-level random generators from STCOS_METHOD_SEED."""

import os
import random


seed_text = os.environ.get("STCOS_METHOD_SEED")
if seed_text:
    seed = int(seed_text)
    random.seed(seed)

    try:
        import numpy as np

        np.random.seed(seed)
    except ImportError:
        pass

    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ImportError:
        pass
