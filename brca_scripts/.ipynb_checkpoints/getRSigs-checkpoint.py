import os
from rpy2.robjects import pandas2ri, conversion
import rpy2.robjects as robjects
import pandas as pd
import numpy as np


def getRSigs():

    # Activate the pandas2ri conversion
    pandas2ri.activate()
    
    # Path to the .rds file
    rds_path = os.path.join(os.getcwd(), 'gene_signatures', 'all_signature_list_symbol.Rds')
    
    # Load the readRDS function from R
    readRDS = robjects.r['readRDS']
    
    # Read the .rds file
    r_object = readRDS(rds_path)
    
    signatures_dict = {}
    singatures_groups_list = []
    
    for i in range(len(r_object)):
        signatures_dict[r_object.names[i]] = r_object[i].tolist()
        singatures_groups_list.append(r_object.names[i])
    
    return signatures_dict, singatures_groups_list