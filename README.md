# Batch corrected t-SNE
joint work with David B. Dunson, Kirk C. Wilhelmsen et al.

### Motivation
In a variety of biomedical applications, there is growing interest in removing systematic biases fromhigh-dimensional data. Some examples include systematic variations and different protocols that compromise reliable sample comparison.Without proper correction, low-dimensional representation of high-dimensional data might encode and reproduce the same systematic variations observed in the original data, and compromise the interpretation of the results. For example, in single-cell transcriptional profiling (scRNAseq) analysis the observed data is sparse, effected by technical variation and cell specific cell-specific variation in capture efficiency. Normalization, clustering of cell type specific data and imputation is used estimate cell-type specific gene expression levels. Currently available methods attempt to remove technical variation by supervised ignoring less significant principal components or independent components that are highly correlated with known “nuisance” covariates.

In the article associated with this repository, we propose two simple and scalable procedures to remove batch effects prior to dimensionality reduction. We demonstrate the utility of this approach using t-SNE a very popular method for dimensionality reduction and data visualisation.

### Results
Our methods are able to adjust different types of batch effects, which might be represented as categorical or continuous measurements. The procedures illustrated in this article are based on algebraic results and constrained optimisation, leading to efficient algorithms and fast computation. The algorithms are applied on single cell gene expression. Results provide interesting insights on cell structures while demonstrating effective in removing nuisance covariates.


# Contents of this repo
`SIMULATION` contains the code to reproduce results from the paper
