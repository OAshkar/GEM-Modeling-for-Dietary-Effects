library(edgeR)
library(ExploreModelMatrix)
library(glimma)
library(tidyverse)
library(data.table)
library(ggh4x)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(grid)

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

# DONE remove ob_ob_HFD_2_Sk # NOTE this caused unbalance and seems unnessary
## grep(pattern = "ob_ob_HFD.*Sk.*", colnames(y), value = T)
## y <- y[,-14]

#################################################
# Diagnostics

y <- calcNormFactors(y, method = "TMMwsp")
lcpm <- cpm(y, log = T)

# We estimate the variance for each row in the logcounts matrix
var_genes <- apply(lcpm, 1, var)
# Get the gene names for the top 500 most variable genes
select_var <- names(sort(var_genes, decreasing = TRUE))[1:200]
# Subset logcounts matrix
highly_variable_lcpm <- lcpm[select_var, ]
dim(highly_variable_lcpm)

dev.new(res = 150, width = 6, height = 9)
png("paper/figure/top200lcpm.png", units = "in", width = 15, height = 25, res = 300)
pheatmap::pheatmap(highly_variable_lcpm,
  name = "LogCount",
  # fontsize_row = 5,
  # fontsize_col = 8,
  # cellheight = 4.5,
  # cellwidth = 8,
  main = "Top 200 Variable Genes",
  show_rownames = F,
  show_colnames = F
)
dev.off()

## dev.off()

png("paper/figure/top100lcpm.png", units = "in", width = 20, height = 18, res = 300)
pheatmap::pheatmap(highly_variable_lcpm,
  fontsize_row = 5,
  fontsize_col = 8,
  cellheight = 4.5,
  cellwidth = 8,
  main = "Top 100 Variable Genes"
)
dev.off()



## png("paper/figure/volcanoTreatment.png", units = "in", width = 8, height = 7, res = 400)
## ggplot(x) +
##   geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = threshold)) +
##   geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
##   geom_vline(xintercept = -1, linetype = "dashed", color = "red") +
##   # geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val), label = ifelse(symbol %in% OkGenes$surv_dge_lasso, symbol, "")), max.overlaps = Inf) +
##   scale_x_continuous(breaks = -10:10) +
##   # ggtitle("Colorectal Cancer Gene Expression") +
##   xlab("log2 fold change") +
##   ylab("-log10 adjusted p-value") +
##   theme(
##     legend.position = "none",
##     plot.title = element_text(size = rel(1.5), hjust = 0.5),
##     axis.title = element_text(size = rel(1.25))
##   ) +
##   theme_classic()
## dev.off()

## png("paper/figure/vocalnoStrain.png", units = "in", width = 8, height = 7, res = 400)
## ggplot(ys) +
##   geom_point(aes(x = logFC, y = -log10(adj.P.Val), colour = threshold)) +
##   geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
##   geom_vline(xintercept = -1, linetype = "dashed", color = "red") +
##   # geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val), label = ifelse(symbol %in% OkGenes$surv_dge_lasso, symbol, "")), max.overlaps = Inf) +
##   scale_x_continuous(breaks = -10:10) +
##   # ggtitle("Colorectal Cancer Gene Expression") +
##   xlab("log2 fold change") +
##   ylab("-log10 adjusted p-value") +
##   theme(
##     legend.position = "none",
##     plot.title = element_text(size = rel(1.5), hjust = 0.5),
##     axis.title = element_text(size = rel(1.25))
##   ) +
##   theme_classic()
## dev.off()


##################################################################################################################################################################
# DGE stratified by tissues .
# volcano plots are generated for visualization
# Need DietI for GSEA, genotypeND and genotypeHFD are used

tissuenames <- unique(y$samples$tissue)
res1 <- data.frame(matrix(ncol = 4, nrow = 0))
de.res <- list()
de.res.HFD <- list()
de.res.diet.wt <- list()
de.res.diet.ob <- list()
AllSet <- list()

png("./paper/figures/Diagnostics.VolcanoPlot.png", units = "in",  width = 12, height = 15, res=200)
par(mfrow = c(7,2))
for (i in seq_along(tissuenames)) {
  print(tissuenames[i])
  tissueMat <- y[, y$samples$tissue == tissuenames[i]]
  tissueMat <- calcNormFactors(tissueMat , method = "TMMwsp")
  factors <- paste(tissueMat$samples$diet, tissueMat$samples$genotype, sep = ".")
  factors2 <- paste(tissueMat$samples$diet, tissueMat$samples$genotype, tissueMat$samples$tissue, sep = ".")
  #############################
  #
  design <- model.matrix(~ 0 + factors)
  cont_matrix <- makeContrasts(
    dietWT = (factorsHFD.WT - factorsND.WT),
    dietob_ob = (factorsHFD.ob_ob - factorsND.ob_ob),
    GenotypeND = (factorsND.ob_ob - factorsND.WT),
    GenotypeHFD = (factorsHFD.ob_ob - factorsHFD.WT),
    StrainI = (factorsND.WT + factorsHFD.WT) / 2 - (factorsND.ob_ob + factorsHFD.ob_ob) / 2, # Strain main effect
    DietI = (factorsND.WT + factorsND.ob_ob) / 2 - (factorsHFD.WT + factorsHFD.ob_ob) / 2, # Treatment main effect
    # interaction = (factorsHFD.ob_ob - factorsND.ob_ob) / 2 - (factorsHFD.WT - factorsND.WT) / 2, # FIXME this one is mathematically meaninless to me, but what does it mean to give a result
    levels = design
  )
  v <- limma::voom(tissueMat, design)
  fit <- lmFit(v, design)
  vfit <- contrasts.fit(fit, cont_matrix)
  efit <- eBayes(vfit)

  limma::volcanoplot(efit, coef = 3, names = efit$genes$SYMBOL, highlight = 100) + title(paste(tissuenames[i], "Genotype.ND"))
  limma::volcanoplot(efit, coef = 4, names = efit$genes$SYMBOL, highlight = 100) + title(paste(tissuenames[i], "Genotype.HFD"))

  # dat <- topTable(efit, coef = "GenotypeND",  n = 2000)
  # xlsx::write.xlsx(dat, "DiffExp_ObOb.xlsx", sheetName = tissuenames[i], row.names = F, append = T)
  #################################################################################
  # Summary Table
  # 7 rows. 7 cols (the coefficients)
  # Each cell is a split of up and
  results <- decideTests(efit, p.value = 0.5, lfc = log2(1.5))
  print(summary(results))
  results <- data.frame(summary(results))
  res2 <- results %>%
    dplyr::filter(Var1 != "NotSig", Var2 != "interaction") %>%
    group_by(Var2) %>%
    mutate(tissue = tissuenames[i])
  res1 <- rbind(res1, res2)
  colnames(res1) <- colnames(res2)
  diffTable <- topTable(efit, p.value = 1, number = Inf, coef = 3, lfc = 0) %>% # Genotype ND with adjp<0.1
    dplyr::rename(gene = SYMBOL) %>%
    # dplyr::rename(log.fc = logFC) %>%
    dplyr::rename(padj = adj.P.Val) %>%
    dplyr::rename(pval = P.Value) %>%
    dplyr::select(gene, logFC, padj, pval, t) %>%
    mutate(rank = rank(logFC, ties.method = "random"))
  de.res[[tissuenames[i]]] <- as.data.table(diffTable)
  ## Genotype HFD
  diffTable.HFD <- topTable(efit, p.value = 1, number = Inf, coef = 4, lfc = 0) %>% # Genotype ND with adjp<0.1
    dplyr::rename(gene = SYMBOL) %>%
    # dplyr::rename(log.fc = logFC) %>%
    dplyr::rename(padj = adj.P.Val) %>%
    dplyr::rename(pval = P.Value) %>%
    dplyr::select(gene, logFC, padj, pval, t) %>%
    mutate(rank = rank(logFC, ties.method = "random"))
  de.res.HFD[[tissuenames[i]]] <- as.data.table(diffTable.HFD)
  ######################################################################################
  # Diet.ND
  diffTable.diet.wt <-
      ## topTable(efit, p.value = 1, number = 1000, coef = 1, lfc = log2(1.5)) %>% # Genotype ND with adjp<0.1
      topTable(efit, p.value = 1, number = Inf, coef = 1) %>% # Genotype ND with adjp<0.1
      dplyr::rename(gene = SYMBOL) %>%
                                        # dplyr::rename(log.fc = logFC) %>%
      dplyr::rename(padj = adj.P.Val) %>%
      dplyr::rename(pval = P.Value) %>%
      dplyr::select(gene, logFC, padj, pval, t) %>%
      mutate(rank = rank(logFC, ties.method = "random"))
  de.res.diet.wt[[tissuenames[i]]] <- as.data.table(diffTable.diet.wt)
  ######################################################################################
  # Diet.ob
  diffTable.diet.ob <-
      ## topTable(efit, p.value = 1, number = 1000, coef = 2, lfc = log2(1.5)) %>% # Genotype ND with adjp<0.1
      topTable(efit, p.value = 1, number = Inf, coef = 2) %>% # Genotype ND with adjp<0.1
    dplyr::rename(gene = SYMBOL) %>%
    # dplyr::rename(log.fc = logFC) %>%
    dplyr::rename(padj = adj.P.Val) %>%
    dplyr::rename(pval = P.Value) %>%
    dplyr::select(gene, logFC, padj, pval, t) %>%
    mutate(rank = rank(logFC, ties.method = "random"))
  de.res.diet.ob[[tissuenames[i]]] <- diffTable.diet.ob
  ##########################################################################################3
  # Filter top 5% genes from each of the 4 contrasts by 1.5 FC
  top_n <- round(nrow(y)*0.05) # 1061
  AllSet[[paste0("dietWT.", tissuenames[i])]] <- topTable(efit, p.value = 1, number = top_n, coef = 1, lfc = log2(1.5)) %>% rownames()
  AllSet[[paste0("dietob_ob.", tissuenames[i])]] <- topTable(efit, p.value = 1, number = top_n, coef = 2, lfc = log2(1.5)) %>% rownames()
  AllSet[[paste0("GenotypeND.", tissuenames[i])]] <-   topTable(efit, p.value = 1, number = top_n, coef = 3, lfc = log2(1.5)) %>% rownames()
  AllSet[[paste0("GenotypeHFD.", tissuenames[i])]] <-  topTable(efit, p.value = 1, number = top_n, coef = 4, lfc = log2(1.5)) %>% rownames()
}
dev.off()

lapply(de.res, nrow)
save(de.res.diet.wt, de.res.diet.ob, file = "ResultsDge/de.res.new.RData")
saveRDS(AllSet, file = "./ResultsDge/AllSet.RDS")
#########################

## load("./ResultsDge/de.res.new.RData")

#### Summary Table
library(flextable)

res <- res1 %>%
  mutate(Var2 = paste0(Var2, ".", Var1)) %>%
  dplyr::select(-Var1) %>%
  group_by(tissue)

my_header <- data.frame(
  col_keys = c(
    "tissue", "dietWT.Down", "dietWT.Up",
    "dietob_ob.Down", "dietob_ob.Up",
    "GenotypeND.Down", "GenotypeND.Up", "GenotypeHFD.Down", "GenotypeHFD.Up",
    "StrainI.Down", "StrainI.Up", "DietI.Down", "DietI.Up"
  ),
  line2 = c(
    "Tissue", rep("Diet WT", 2), rep("Diet Ob/Ob", 2), rep("Genotype ND", 2), rep("Genotype HFD", 2),
    rep("Genotype Intr*", 2), rep("Diet Intr**", 2)
  ),
  line3 = c("Tissue", rep(c("Up", "Down"), 6))
)
ress <- tidyr::pivot_wider(as_tibble(res), names_from = Var2, values_from = Freq)

ft <- flextable::flextable(ress, col_keys = my_header$col_keys) %>%
  set_header_df(
    mapping = my_header,
    key = "col_keys"
  ) %>%
  # add_footer("* Interaction") %>%
  theme_booktabs() %>%
  merge_v(part = "header") %>%
  merge_h(part = "header") %>%
  align(align = "center", part = "all") %>%
  autofit() %>%
  empty_blanks()
ft
save_as_image(x = ft, path = "paper/figures/DGEsummaryTable_by_tissue.png", webshot = "webshot2")

tissue.labs <- tissues
names(tissue.labs) <- tissue

effect.lab <- c("Diet.WT", "Diet.Ob/Ob", "Diet Effect", "Genotype.ND", "Genotype.HFD", "Genotype Effect")
names(effect.lab) <- c("dietWT", "dietob_ob", "DietI", "GenotypeND", "GenotypeHFD", "StrainI")

png("./paper/figures/DGE.effectsizes.png", units = "in",  width = 10, height = 4, res=150)
tidyr::pivot_wider(res1, names_from = Var1, values_from = Freq) %>%
    mutate(Var2  = fct_relevel(Var2, "DietI", after = 2)) %>%
    ggplot() +
    geom_bar(aes(y = tissue, x = -Down),
             stat = "identity", position = "dodge", fill = "#F15412", width = 1) +
    geom_bar(aes(y = tissue, x = Up),
             stat = "identity", position = "dodge", fill = "#34B3F1", width = 1) +
    scale_x_continuous(labels = scales::label_number(suffix = " K", scale = 1e-3)) +
    ggh4x::facet_nested(tissue~Var2, scale = "free", switch = "y", labeller = labeller(tissue = tissue.labs, Var2 = effect.lab)) +
    labs(x = "Differential Gene Expression Contrasts Up/Down Counts", y = "Tissues") +
    theme_bw() +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_blank(),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.background = element_rect(fill="#EEEEEE"),
          strip.text.y.left = element_text(angle = 0))
dev.off()

