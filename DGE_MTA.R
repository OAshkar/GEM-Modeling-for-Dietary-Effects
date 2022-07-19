library(edgeR)
library(tidyverse)
library(org.Mm.eg.db)
library(data.table)

load("./counts/Count.Rdata")

y <- DGEList(counts)
anno <- select(org.Mm.eg.db,
  keys = row.names(y$counts),
  columns = c("ENTREZID", "GENETYPE"),
  keytype = "SYMBOL"
)
y$genes <- anno
y$samples <- y$samples %>%
  mutate(genotype = gsub(x = colnames(y$counts), pattern = "(ob_ob|WT)(.*)", replacement = "\\1")) %>%
  mutate(diet = gsub(x = colnames(y$counts), pattern = "(.*)(HFD|ND|HDF)(.*)", replacement = "\\2")) %>%
  mutate(rep = gsub(x = colnames(y$counts), pattern = "(.*)_(\\d)_(.*)", replacement = "\\2")) %>%
  mutate(tissue = gsub(x = colnames(y$counts), pattern = "(.*)(Ao|Hi|Sk|He|Li|Ep|Hy)(.*)", replacement = "\\2"))
y$samples <- y$samples %>%
  mutate(diet = gsub(x = y$samples$diet, pattern = "HDF", replacement = "HFD"))
factors <- paste(y$samples$diet, y$samples$genotype, sep = ".")
factors2 <- paste(y$samples$diet, y$samples$genotype, y$samples$tissue, sep = ".")
y$samples <- y$samples %>% mutate(
  factors = paste(y$samples$diet, y$samples$genotype, sep = "."),
  factors2 = paste(y$samples$diet, y$samples$genotype, y$samples$tissue, sep = ".")
)
# filter other than coding genes.
y <- y[which(y$genes$GENETYPE == "protein-coding"), ]
y <- calcNormFactors(y)


tissuenames <- unique(y$samples$tissue)
for (i in seq_along(tissuenames)) {
  print(tissuenames[i])
  tissueMat <- y[, y$samples$tissue == tissuenames[i]]
  tissueMat <- calcNormFactors(tissueMat)
  factors <- paste(tissueMat$samples$diet, tissueMat$samples$genotype, sep = ".")
  factors2 <- paste(tissueMat$samples$diet, tissueMat$samples$genotype, tissueMat$samples$tissue, sep = ".")
  #############################
  #
  design <- model.matrix(~ 0 + factors)
    # Contrasts
  cont_matrix <- makeContrasts(
    # WTHFD - WTND
    dietWT = (factorsHFD.WT - factorsND.WT),
    # ObObHFD - WTND
    dietgenotype = (factorsHFD.ob_ob - factorsND.WT),
    # ObObND - WTND
    GenotypeND = (factorsND.ob_ob - factorsND.WT),
    # interaction = (factorsHFD.ob_ob - factorsND.ob_ob) / 2 - (factorsHFD.WT - factorsND.WT) / 2, # FIXME this one is mathematically meaninless to me, but what does it mean to give a result
    levels = design
  )
  v <- limma::voom(tissueMat, design)
  fit <- lmFit(v, design)
  vfit <- contrasts.fit(fit, cont_matrix)
  efit <- eBayes(vfit)

  for (ii in seq_along(colnames(efit$contrasts))) {
    topTable(efit, p.value = 1, number = Inf, coef = ii) %>%
      dplyr::rename(gene = SYMBOL) %>%
      # dplyr::rename(log.fc = logFC) %>%
      dplyr::rename(padj = adj.P.Val) %>%
      dplyr::rename(pval = P.Value) %>%
      dplyr::select(gene, logFC, padj, pval, t) %>%
      mutate(rank = rank(logFC, ties.method = "random")) %>%
      write.csv(paste0("rMTAResults/dge/", tissuenames[i], "_", colnames(efit$contrasts)[ii], ".csv"))
  }
}
lapply(de.res, nrow)

