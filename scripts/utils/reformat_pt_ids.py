import numpy as np
import pandas as pd
from utils.consts import *

def reformat_pt_ids(col):
    reformatted_data = {}
    for idx in col.index:
        pt_id_raw = col.loc[idx]
        # Remove any spaces
        pieces = pt_id_raw.split(' ')
        pt_id = ''
        for piece in pieces:
            pt_id += piece
        # Make sure mhbb is all lowercase if it's there
        # Else make sure patient IDs are 6-digit
        # with leading zeros
        if pt_id[0] in ['m','M']:
            pt_id='mhbb'+pt_id[4:]
        else:
            pt_id = pt_id.rjust(NUMERIC_PT_ID_WIDTH,'0')
        reformatted_data[idx]=pt_id
    pt_ids_reformatted = pd.Series(data = reformatted_data)
    return pt_ids_reformatted

