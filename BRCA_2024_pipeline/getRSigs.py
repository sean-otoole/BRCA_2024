import os  # Import the os module to interact with the operating system
from rpy2.robjects import pandas2ri, conversion  # Import rpy2 for R-Python interface
import rpy2.robjects as robjects  # Import robjects for interacting with R from Python
import pandas as pd  # Import pandas for working with dataframes
import numpy as np  # Import numpy for numerical operations
import mygene  # Import the mygene package for querying gene information

def getRSigs():
    # Activate the pandas2ri conversion to convert between pandas and R dataframes
    pandas2ri.activate()

    # Get the current working directory
    current_directory = os.getcwd()
    
    # Define the path to the .rds file containing the gene signatures
    rds_path = current_directory + '/gene_signatures/all_signature_list_ensembl.Rds'
    print(rds_path)  # Print the path to ensure it's correct
    
    # Load the R function readRDS using robjects
    readRDS = robjects.r['readRDS']
    
    # Read the .rds file to get the gene signatures in R
    r_object = readRDS(rds_path)
    
    # Initialize empty dictionary and list to store signature data
    signatures_dict = {}
    singatures_groups_list = []
    
    # Iterate over the elements in the R object and store them in a dictionary
    for i in range(len(r_object)):
        signatures_dict[r_object.names[i]] = r_object[i].tolist()  # Convert R list to Python list
        singatures_groups_list.append(r_object.names[i])  # Append group names to the list

    # Initialize MyGeneInfo object to query gene symbols
    mg = mygene.MyGeneInfo()
    
    # Dictionary to store symbols corresponding to Ensembl IDs
    symbols_dict_from_ensembl = {}
    
    # Iterate over each key in the signatures dictionary
    for key in signatures_dict:
        ensembl_ids = signatures_dict[key]  # Get the Ensembl IDs for the current signature
        symbol_list = []  # Initialize an empty list to store the symbols
        # Query MyGeneInfo to get the gene symbols for the Ensembl IDs
        results = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species='human')
        # Extract the gene symbols from the results, defaulting to 'No Symbol' if none is found
        symbols = [item.get('symbol', 'No Symbol') for item in results]
        symbols_dict_from_ensembl[key] = symbols  # Store the symbols in the dictionary
    
    # Return the dictionary of symbols and the list of signature groups
    return symbols_dict_from_ensembl, singatures_groups_list
