## Simulation study

Reproduce the simulations running the script `./SIMULATION_revision.R`, doing for example

```R
Rscript ./SIMULATION_julia.R
```
Note that a `julia` exectuable is also required (i.e. `julia`'s path should be in your `$PATH` variables and you are able to call `julia` from the command line).
Julia can be freely downloaded at [its official page](https://julialang.org).

The scripts writes the file `sim.pdf` which correspond to Figure 1 of the paper. 


Tested with

```
julia version 1.1.1
```
and

```
R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] scater_1.10.1               ggplot2_3.2.0
 [3] splatter_1.6.1              SingleCellExperiment_1.4.1
 [5] SummarizedExperiment_1.12.0 DelayedArray_0.8.0
 [7] BiocParallel_1.16.6         matrixStats_0.54.0
 [9] Biobase_2.42.0              GenomicRanges_1.34.0
[11] GenomeInfoDb_1.18.2         IRanges_2.16.0
[13] S4Vectors_0.20.1            BiocGenerics_0.30.0

loaded via a namespace (and not attached):
 [1] beeswarm_0.2.3           tidyselect_0.2.5         locfit_1.5-9.1
 [4] purrr_0.3.2              reshape2_1.4.3           HDF5Array_1.10.1
 [7] lattice_0.20-38          rhdf5_2.26.2             colorspace_1.4-1
[10] viridisLite_0.3.0        rlang_0.4.0              pillar_1.4.2
[13] glue_1.3.1               withr_2.1.2              GenomeInfoDbData_1.2.0
[16] plyr_1.8.4               stringr_1.4.0            zlibbioc_1.28.0
[19] munsell_0.5.0            gtable_0.3.0             labeling_0.3
[22] vipor_0.4.5              Rcpp_1.0.1               scales_1.0.0
[25] backports_1.1.4          checkmate_1.9.4          XVector_0.22.0
[28] gridExtra_2.3            digest_0.6.20            stringi_1.4.3
[31] dplyr_0.8.3              grid_3.6.1               tools_3.6.1
[34] bitops_1.0-6             magrittr_1.5             lazyeval_0.2.2
[37] RCurl_1.95-4.12          tibble_2.1.3             crayon_1.3.4
[40] pkgconfig_2.0.2          Matrix_1.2-17            DelayedMatrixStats_1.4.0
[43] ggbeeswarm_0.6.0         assertthat_0.2.1         viridis_0.5.1
[46] Rhdf5lib_1.4.3           R6_2.4.0                 compiler_3.6.1


```
