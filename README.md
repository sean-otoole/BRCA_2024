# Spatial transcriptomics analysis of pre- and postmenopausal breast cancer samples:

For a visual summary of this project please take a look at the [report](https://sean-otoole.github.io/BRCA_2024/) associated with this repository.

<sub>**Please note** that this repository, even with the appropriate libraries and packages installed, will not operate independently. For confidentiality reasons the **original datasets** are not included.
<br>

## Project Organization
```

┌── BRCA_2024_pipeline/                            : contains the main pipeline
│   ├── BRCA_report.Rmd                            : Generates the report file
│   ├── acquirePseudoBulkRef.py                    : Acquires the breast cancer reference data set and constructs pseudo bulk profiles from the reference data
│   ├── applyCosineScores.py                       : Calculates cosine similarity scores across the tissue for a number of cell types from the reference data set
│   ├── applyImageMask.py                          : Applies an image mask to the visium data to restrict the analysis to tumors based on pathology recommendations
│   ├── geneSigCosSimHeatmap.py                    : Generates a heatmap correlating cosine similarity scores with gene signature values
│   ├── getRSigs.py                                : Acquires gene signatures utilizing the AUCell method
│   ├── helpFuncs.py                               : Contains various helper functions
│   ├── label_transfer_seurat_deprecrated.r        : No longer relevant script which had utilized a Seurat-based integration technique
│   ├── ma_et_al_ref_analysis.r                    : Processes and acquires the b cell subtype data set
│   ├── preProcessPipeline.py                      : Processes all of the visium samples
│   ├── pyAUCell.py                                : A python implementaiton of the AUcell method for calculating gene signatures
│   ├── sctransform.r                              : Script for calling the transform normalization technique which regresses out sources of technical noise.
│   ├── spatial_metric_description.r               : Graphic script which will display example images and calculate statistics given an already present metric
│   ├── summary_stats.r                            : Quickly displays some relevant stats
│   ├── visiumPythonPipeline.py                    : Pipeline for calling the script written in python.
│   ├── reports/
│   │   └── BRCA_report.html                       : The html file produced by the markdown script which can be viewed using the link mentioned above.
│   └── visiumRPipeline                            : A script for calling the scripts that were written in R
├── BRCA_nbs_previous/                             : contains previous notebook files for prototyping pipeline scripts
├── BRCA_scripts/                                  : contains early prototypes for scripts of the main pipeline                        
|── figures/
├── gene_signatures/                               : list of genes associated with various immune cell type-related gene signatures
├── index.html                               : report file generated with R markdown script produced with Python and R code from the main pipeline
└── README.md                                      : project description

```

<br>
