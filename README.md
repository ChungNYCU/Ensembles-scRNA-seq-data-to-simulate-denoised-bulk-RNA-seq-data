# Ensembles-scRNA-seq-data-to-simulate-denoised-bulk-RNA-seq-data
Sample data can download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141834 <br>
Before running the code, you should install umap and scanpy. <br>
Umap: https://umap-learn.readthedocs.io/en/latest/ <br>
Scanpy: https://scanpy.readthedocs.io/en/stable/ <br>
Version: scanpy==1.7.0 anndata==0.7.5 umap==0.5.1 numpy==1.19.2 scipy==1.5.2 pandas==1.1.3 scikit-learn==0.23.2 statsmodels==0.12.0 <br>
```bash
├── 00_select_scRNA_and_bulkRNA_co_gene.py
├── 01_data_preprocess.py
├── 02_ensembled_data.py
├── README.md
└── dataset
    ├── GSE141834_bulkRNAseq_normalized_counts.txt
    └── GSE141834_scRNAseq_seuratV3_normalized.txt
```
