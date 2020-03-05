suppressPackageStartupMessages({
  library(bcTSNE)
  library(data.table)
  library(SingleCellExperiment)
  library(RSpectra)
  library(batchelor)
  library(kBET)
  library(Rtsne)
  library(lisi)
  library(harmony)
  library(dlfUtils)
  library(scater)
})

path <- "mouse_raw_data_allcells_Notscaled.tsv"
dat <- fread(path)

## Subset to veh only
dat <- dat[TREATMENT == "Vehicle"]
dat[ , TREATMENT := NULL]

mCols <- c("CELL", "SEX", "DATE", "NODEID", "MOUSE")
mDat <- dat[ , .SD, .SDcols = mCols]
mDat[ , CELLTYPE := NODEID]
mDat[grepl("Node", NODEID),     CELLTYPE := "Other"]
mDat[grepl("Neurons", NODEID),  CELLTYPE := "Other"]
mDat[grepl("Vascular", NODEID), CELLTYPE := "Other"]

cCols <- setdiff(names(dat), mCols)
cDat <- data.table(GENE = cCols)

tmp <- t(as.matrix(dat[ , .SD, .SDcols = cCols]))
sce <- SingleCellExperiment(assays = list(counts = tmp),
                            colData = mDat,
                            rowData = cDat)
sizeFactors(sce) <- librarySizeFactors(sce)
sce <- normalize(sce, return_log = FALSE)
sce <- normalize(sce)
assay(sce, "centered") <- t(scale(t(normcounts(sce)), scale = FALSE))

SVD <- svds(A = t(assay(sce, "centered")), k = 50)
reducedDim(sce) <- SVD$u %*% diag(SVD$d)
reducedDimNames(sce) <- "PCA"

saveRDS(sce, "real.sce")
rm(dat, tmp, mDat, cDat, SVD, mCols, cCols); gc();

res <- vector(mode = "list", length = 4)
names(res) <- c("tsne", "bcts", "hmny", "mnn")

## Calculate unadjusted t-SNE
set.seed(1234)
res$tsne <- Rtsne(reducedDim(sce), pca = FALSE)

## Calculate BC-t-SNE
set.seed(1234)
res$bcts <- bctsne(X = t(assay(sce, "centered")), 
                   Z = model.matrix( ~ -1 + factor(colData(sce)$MOUSE)), 
                   k = 50, 
                   perplexity = 30, 
                   maxIter = 1000)


## Calcualte Harmony
set.seed(1234)
sce <- RunHarmony(sce, group.by.vars = "MOUSE")
res$hmny = Rtsne(as.matrix(reducedDim(sce, "HARMONY")), pca = FALSE)

## Calculate MNN
set.seed(1234)
mnnSCE <- fastMNN(sce, batch = factor(colData(sce)$MOUSE))
res$mnn <- Rtsne(reducedDim(mnnSCE), pca = FALSE)


res <- lapply(res, "[[", "Y")
saveRDS(res, "real.res")

## Function to calculate separation metrics
calcMetrics <- function(Y, bchLst) {
  calcSil <- function(x) {
    s <- batch_sil(pca.data = list(x = Y), batch = x, nPCs = 2)
    1 - abs(s)
  }
  calcKBET <- function(x) {
    kBET(Y, batch = x, do.pca = FALSE, plot = FALSE)$average
  }
  calcPCA <- function(x) {
    pcRegression(pca.data = prcomp(Y), batch = x, n_top = 2)$pcReg
  }
  sil  <- sapply(bchLst, calcSil)
  kbet <- sapply(bchLst, calcKBET)
  lisi <- compute_lisi(Y, 
                       meta_data = as.data.frame(bchLst),
                       label_colnames = names(bchLst))
  lisi <- colMeans(lisi)
  sizes <- sapply(bchLst, function(x) length(unique(x)))
  lisi <- (lisi - 1)/(sizes - 1)
  pca  <- sapply(bchLst, calcPCA)
  res <- list(sil = sil, kbet = kbet, lisi = lisi, pca = pca)
  do.call(cbind, res)
}

bchLst <- as.list(colData(sce)[ , c("SEX", "DATE", "MOUSE", "CELLTYPE")])
bchLst <- lapply(bchLst, factor)

metrics <- lapply(res, calcMetrics, bchLst = bchLst)

saveRDS(metrics, "real.metrics")








