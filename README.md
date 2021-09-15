# Ensembles-scRNA-seq-data-to-simulate-denoised-bulk-RNA-seq-data

## Thesis Link: https://github.com/ChungNYCU/Ensembles-scRNA-seq-data-to-simulate-denoised-bulk-RNA-seq-data/blob/master/Weiche%20Chung%20Master%20Thesis.pdf

## Abstract
This research attempts to develop an algorithm to simulate traditional bulk RNA sequencing (bulk RNA-seq) data by single-cell RNA sequencing (scRNA-seq) data. scRNA-seq was first published in 2009. Nowadays, due to the decline in sequencing costs, a large amount of scRNA- seq data have been generated and provided to researchers. The emergence of scRNA-seq allows researchers to study new biology issues. For example, identification of cell types and the study of heterogeneous cells, such as cancer cells, are due to scRNA-seq. This allows us to observe expression in specific cells. Currently, there are several traditional bulk RNA-seq datasets and many methods suitable for traditional bulk RNA-seq but none for scRNA-seq, such as the Gene Co-Expression Analysis Method. If we compare the gene expression matrix of traditional bulk RNA-seq and scRNA-seq, the gene expression matrix of traditional RNA sequencing has less noise and fewer missing values. For scRNA-seq, the experimental procedures are more sophisticated, resulting in greater noise. At the same time, the depth of sequencing is shallower, producing more missing values. Because the two RNA sequencing methods are the same, this study tried to reduce the noise and missing values of the scRNA-seq data to mimic traditional bulk RNA-seq data using the algorithm way. <br>
We can generate ensembled data for each scRNA-seq sample, and every sample has the same attributes with bulk RNA-seq, include noise, missing value decrease and sample’s correlation, low expression gene value increase. <br>
In conclusion, our method can use in diverse situations, like comparing existing bulk RNA-seq data and scRNA-seq. From the perspective of analysis, we can generate a lot of purely bulk RNA-seq samples, such as turn cancer cells into pure cancer tissue. <br>

## Package version and Dataset 
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
