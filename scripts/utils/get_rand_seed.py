import numpy as np
from utils.consts import *

# Helper function to pull an
# integer to use as a random seed from the
# current random number generator
def get_rand_seed(random_number_generator):
    rand_seed_curr = random_number_generator.randint(
            0,
            INT32_MAX
            )
    return rand_seed_curr

