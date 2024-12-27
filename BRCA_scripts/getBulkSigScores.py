import os
from rpy2.robjects import pandas2ri, conversion
import rpy2.robjects as robjects
import pandas as pd
import numpy as np
import mygene

def getBulkSigs():

    # Activate the pandas2ri conversion
    pandas2ri.activate()
    
    # Path to the .rds file
    rds_path = os.path.join(os.getcwd(), 'gene_signatures', 'all_signature_scores.Rds')
    
    # Load the readRDS function from R
    readRDS = robjects.r['readRDS']
    
    # Read the .rds file
    r_object = readRDS(rds_path)
    
    return r_object