# Batch corrected t-SNE
Biomedical research often produces high-dimensional data confounded by batch effects such as systematic experimental variations, different protocols and subject identifiers.
Without proper correction, low-dimensional representation of high-dimensional data might encode and reproduce the same systematic variations observed in the original data, and compromise the interpretation of the results.
In this article, we propose a novel procedure to remove batch effects from low-dimensional embeddings obtained with t-SNE dimensionality reduction.
The proposed methods are based on linear algebra and constrained optimization, leading to efficient algorithms and fast computation in many high-dimensional settings.
Results on artificial single-cell transcription profiling data show that the proposed procedure successfully removes multiple batch effects from t-SNE embeddings, while retaining fundamental information on cell types.
When applied to single-cell gene expression data to investigate mouse medulloblastoma, the proposed method successfully removes batches related with mice identifiers and the date of the experiment, while preserving clusters of oligodendrocytes, astrocytes, and endothelial cells and microglia, which are expected to lie in the stroma within or adjacent to the tumors.

This repository is associated with the paper: Aliverti, Tilson,  Filer, Babcock, Colaneri, Ocasio, Gershon, Wilhelmsen, and Dunson, D (2019). ['Batch-correction of high dimensional data'](https://arxiv.org/abs/1911.06708)

# Contents 
- `SIMULATION` contains the code to reproduce results from the paper
- `./tsne_julia/bctsne.jl` Julia implementation of bcTsne. It can be loaded direcly as a module. See `./SIMULATION/run_sim.jl` for an example.

```
.
├── README.md
├── SIMULATION
│   ├── README.md
│   ├── run_sim.jl
│   ├── sim.pdf
│   └── SIMULATION_julia.R
└── tsne_julia
    └── bctsne.jl
```
