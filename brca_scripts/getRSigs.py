import os
from rpy2.robjects import pandas2ri, conversion
import rpy2.robjects as robjects
import pandas as pd
import numpy as np
import mygene

def getRSigs():
    # Activate the pandas2ri conversion
    pandas2ri.activate()
    
    # Path to the .rds file
    rds_path = '/tmp/work/Visium/BRCA_2024/gene_signatures/all_signature_list_ensembl.Rds'
    # rds_path = os.path.join(os.getcwd(), 'gene_signatures', 'all_signature_list_ensembl.Rds')
    
    # Load the readRDS function from R
    readRDS = robjects.r['readRDS']
    
    # Read the .rds file
    r_object = readRDS(rds_path)
    
    signatures_dict = {}
    singatures_groups_list = []
    
    for i in range(len(r_object)):
        signatures_dict[r_object.names[i]] = r_object[i].tolist()
        singatures_groups_list.append(r_object.names[i])

    mg = mygene.MyGeneInfo()
    
    symbols_dict_from_ensembl = {}
    
    for key in signatures_dict:
        ensembl_ids = signatures_dict[key]
        symbol_list = []
        results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
        symbols = [item.get('symbol', 'No Symbol') for item in results]  #will add 'No symbol' string if no symbol key is present
        symbols_dict_from_ensembl[key] = symbols
    
    return symbols_dict_from_ensembl, singatures_groups_list