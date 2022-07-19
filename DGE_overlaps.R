library(tidyverse)
library(ComplexHeatmap)
library(dtplyr)
library(data.table)
library(patchwork)

AllSet.l <- readRDS("./ResultsDge/AllSet.RDS")
load("./ResultsDge/de.res.new.RData")

AllSet <- list_to_matrix(AllSet.l) %>% as.data.frame()

##############################################################################################################################################
# DEPRECATED Overlap All
m2 <- make_comb_mat(AllSet, mode = "intersect",
                    min_set_size = 15)
UpSet(m2[comb_size(m2) > 4]) #

comb_name(m2)
comb_degree(m2)
comb_size(m2)

##############################################################################################################################################
# Overlap 4 contrasts ~ strata(tissue)
tissue <- AllSet %>% colnames() %>% str_split("\\.", simplify = T) %>% .[,2] %>% unique()
tissues <- case_when(
    tissue == "Ao" ~ "Aorta",
    tissue == "Ep" ~ "Epididymal Fat",
    tissue == "He" ~ "Heart",
    tissue == "Li" ~ "Liver",
    tissue == "Hy" ~ "Hypothalamus",
    tissue == "Hi" ~ "Hippocampus",
    tissue == "Sk" ~ "Skeletal Muscle",
    )

tissue.fig <- list()
tis.inter.dietEff <- list()
for( i in seq_along(tissue)){
    m2 <- make_comb_mat(AllSet %>% dplyr::select(contains(tissue[i])), mode = "distinct")
    m2 <- m2[comb_degree(m2) > 1]
    tissue.fig[[i]] <- UpSet(m2,
                             right_annotation = rowAnnotation(MainEffect = c("Diet",
                                                                             "Diet",
                                                                             "Genotype","Genotype"),
                                                              show_legend = ifelse(i == 7, T , F),
                                                              col = list(MainEffect = c("Diet" = "#F15412",
                                                                                        "Genotype" = "#34B3F1")),
                                                              show_annotation_name = F
                                                              ),
                             column_title = tissues[i],
                             row_labels = c("Diet WT", "Diet Ob/Ob", "Genotype ND", "Genotype HFD"),
                             height = unit(40, "mm"),
                             width = unit(55, "mm"),
                             comb_col = c("black", "#F15412", "#34B3F1")[c(rep(1,5),2, rep(1,4), 3)],
                             lwd = 2,
                             column_title_gp = gpar(fontsize = 21),
                             column_names_gp = grid::gpar(fontsize = 15),
                             row_names_gp = grid::gpar(fontsize = 15)
                             ) %>% draw() %>% grid::grid.grabExpr()
    tis.inter.dietEff[[tissue[i]]] <- extract_comb(m2, "1100") # TODO GSEA
}

png(file = "./paper/figures/overlap_by_tissue.png", units =  "in", width=17, height = 10, res = 150)
## append(tissue.fig, space, after = 5 ) %>%
    wrap_plots(tissue.fig, ncol = 3, guides = "keep")
dev.off()

##############################################################################################################################################
# Overlap 7 tissue ~ strata(contrast)
contrast <- AllSet %>% colnames() %>% str_split("\\.", simplify = T) %>% .[,1] %>% unique()
contrast2 <- c("Diet.WT", "Diet.Ob/Ob", "Genotype.ND", "Genotype.HFD")
contrast.fig <- list()
inter.contrast <- list()

for( i in seq_along(contrast)){
    m3 <- make_comb_mat(AllSet %>% dplyr::select(contains(contrast[i])), mode = "distinct")
    m3 <- m3[comb_degree(m3) > 3]
    contrast.fig[[i]] <- UpSet(m3,
                               right_annotation = NULL,
                               column_title = contrast2[i],
                               row_labels = tissues,
                               height = unit(50, "mm"),
                               width = unit(150, "mm"),
                               comb_col = c(ifelse(i ==1, "black", "#F15412"), rep("black", ncol(m3)-1)),
                               lwd = 2,
                               column_title_gp = gpar(fontsize = 25),
                               column_names_gp = grid::gpar(fontsize = 15),
                               row_names_gp = grid::gpar(fontsize = 15),
                               ) %>%
        draw() %>%  grid::grid.grabExpr()
    if(i != 1) inter.contrast[[contrast[i]]] <- extract_comb(m3, names(comb_degree(m3))[comb_degree(m3) == 7])
}

png(file = "./paper/figures/overlap_by_contrast.png", width = 16, height = 8, res = 150, units = "in")
wrap_plots(contrast.fig, ncol = 2)
dev.off()

de.common <- unique(unlist(inter.contrast))


tidy.common.lcpm <- as.data.frame(lcpm)[row.names(lcpm) %in% de.common, ] %>%
    rownames_to_column() %>%
    tidyr::pivot_longer(
               cols = -rowname
           ) %>%
    mutate(strain = gsub(x = name, pattern = "(ob_ob|WT)(.*)", replacement = "\\1")) %>%
    mutate(diet = gsub(x = name, pattern = "(.*)(HFD|ND|HDF)(.*)", replacement = "\\2")) %>%
    mutate(diet = gsub(x = diet, pattern = "HDF", replacement = "HFD")) %>%
    mutate(genotype = paste0(strain, diet)) %>%
    mutate(rep = gsub(x = name, pattern = "(.*)_(\\d)_(.*)", replacement = "\\2")) %>%
    mutate(tissue = gsub(x = name, pattern = "(.*)(Ao|Hi|Sk|He|Li|Ep|Hy)(.*)", replacement = "\\2")) %>%
    mutate( tissue = case_when(
    tissue == "Ao" ~ "Aorta",
    tissue == "Ep" ~ "Epididymal Fat",
    tissue == "He" ~ "Heart",
    tissue == "Li" ~ "Liver",
    tissue == "Hy" ~ "Hypothalamus",
    tissue == "Hi" ~ "Hippocampus",
    tissue == "Sk" ~ "Skeletal Muscle",
    ))


#################################################################################
png(file = "./paper/figures/DGE.overlap.4genes.png", units = "in", width = 8, height = 7, res = 400)
ggplot(tidy.common.lcpm, aes(x = diet, y = value)) +
    geom_boxplot(aes(x = diet, y = value, fill = strain)) +
    scale_fill_hue(c = 45, l = 45) +
    ## scale_fill_manual(values=c("#F15412", "#34B3F1")) +
    ## geom_jitter(aes(x = diet, y = value, color = strain), binaxis = "y", size = 0.3) +
    facet_nested(rowname ~ tissue+diet, scales = "free") + theme_bw() +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(fill="#EEEEEE")
    ) +
    labs(y = "Log CPM Value")
dev.off()


#################################################################################
## Plot effect size
## de.res.diet.wt de.res.diet.ob

de.table.wt <- as_tibble(Reduce(rbind, lapply(de.res.diet.wt, function(x) x[, c("gene", "logFC", "pval")])))
de.table.wt$tissue <- rep(names(de.res.diet.wt), Reduce(c, lapply(de.res.diet.wt, nrow)))
de.table.wt$Contrast <- "Diet.WT"
## write.csv(
##   de.table.ND %>% filter(tissue == "Ao"),
##   "ResultsDge/diff_exprs.csv"
## )

## de.res2.HFD <- lapply(de.res.HFD, function(x) x$gene)
## de.common.HFD <- Reduce(intersect, de.res2.HFD) ## Genes DE in all sets
## de.common <- unique(c(de.common.HFD, de.common))
de.table.ob <- as.tibble(Reduce(rbind, lapply(de.res.diet.ob, function(x) x[, c("gene", "logFC", "pval")])))
de.table.ob$tissue <- rep(names(de.res.diet.ob), Reduce(c, lapply(de.res.diet.ob, nrow)))
de.table.ob$Contrast <- "Diet.Ob/Ob"

de.table <- rbind(de.table.wt, de.table.ob) %>%
    mutate( tissue = case_when(
    tissue == "Ao" ~ "Aorta",
    tissue == "Ep" ~ "Epididymal Fat",
    tissue == "He" ~ "Heart",
    tissue == "Li" ~ "Liver",
    tissue == "Hy" ~ "Hypothalamus",
    tissue == "Hi" ~ "Hippocampus",
    tissue == "Sk" ~ "Skeletal Muscle",
    ))

png(file = "./paper/figures/DGE.ES.4genes.png", units = "in", width = 10, height = 5, res = 150)
ggplot(
  de.table %>%
    filter(gene %in% de.common) %>%
    group_by(tissue),
  aes(x = tissue, y = logFC, fill = Contrast)
) +
  geom_bar(stat = "identity", width = 0.8,
           position = position_dodge2(width = 1.2, preserve = "single")) +
  facet_nested(gene ~ tissue, scales = "free") +
  theme_bw() +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    ## legend.position="none"
  ) +
  labs(y = "Expression Log Fold Change") +
  scale_fill_hue(c = 45, l = 45)
dev.off()
