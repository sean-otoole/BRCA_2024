# Spatial transcriptomics analysis of pre- and postmenopausal breast cancer samples:

For a visual summary of this project please download and view the HTML file associated with this repository [BRCA_report.html](./BRCA_report.html).

<sub>**Please note** that this repository, even with the appropriate libraries and packages installed, will not operate independently. For confidentiality reasons the **original datasets** are not included.
<br>

## Project Organization
```

┌── BRCA_2024_pipeline/                            : contains the main pipeline
│   ├── BRCA_report.Rmd
│   ├── acquirePseudoBulkRef.py
│   ├── applyCosineScores.py
│   ├── applyImageMask.py
│   ├── geneSigCosSimHeatmap.py
│   ├── getRSigs.py
│   ├── helpFuncs.py
│   ├── label_transfer_seurat_deprecrated.r
│   ├── ma_et_al_ref_analysis.r
│   ├── preProcessPipeline.py
│   ├── pyAUCell.py
│   ├── sctransform.r
│   ├── spatial_metric_description.r
│   ├── summary_stats.r
│   ├── visiumPythonPipeline.py
│   ├── reports/
│   │   └── BRCA_report.html
│   └── visiumRPipeline
├── BRCA_nbs_previous/                             : contains previous notebook files for prototyping pipeline scripts
├── BRCA_scripts/                                  : contains early prototypcs for scripts of the main pipeline                        
|── figures/
├── gene_signatures/                               : list of genes associated with various immune cell type related gene signatures
├── BRCA_report.html                               : report file generated with R mardown script produced with python and R code from the main pipline
└── README.md                                      : project description

```

<br>
