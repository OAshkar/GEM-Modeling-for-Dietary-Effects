library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ComplexUpset)
library(patchwork)
library(grid)
library(tidyverse)
library(ggh4x)
library(furrr)

#load("samplesAllTidy.Rdata") #ACHR
load(file = "samplesAllTidy_ACHR2.Rdata")
transportRxns <- readRDS("./tmp/transportRxns.RDS")
dat <- dat %>% select(-any_of(transportRxns))

#################################################################################
# DEPRECATED Comparison of Flux ObND-WTND & ObHFD-WTHFD
# DONE Comparison of Flux Ob.HFD-Ob.ND & WTHFD-WTND
#################################################################################

df.wilcox <- function(mat0, mat1, padj.cutoff, r.cutoff, df.cutoff, rdf.cutoff) {
  if (!(all(rownames(mat0) == rownames(mat1))) | !(all(colnames(mat0) == colnames(mat1)))) {
    stop("Samples are not consistent")
  }
  dflux.test <- function(s0, s1) {
    # an inner helper function to run wilcoxon's rank-sum test
    ## tryCatch(
    {
      wilcox.res <- wilcox.test(s0, s1)
      # p value
      wilcox.p <- wilcox.res$p.value
      # effect size for wilcoxon's rank-sum test: rank biserial correlation
      wilcox.r <- unname(1 - 2 * wilcox.res$statistic / (sum(!is.na(s0)) * sum(!is.na(s1))))
      # another effect size measure: difference of median fluxes (I use median instead of median to facilitate some specific downstream analyses)
      ## m0 <- median(s0)
      ## m1 <- median(s1)
      m0 <- median(s0, na.rm = T)
      m1 <- median(s1, na.rm = T)
      data.frame(
        rxn = parent.frame(1)$.,
        lb0 = min(s0), ub0 = max(s0), median0 = m0,
        lb1 = min(s1), ub1 = max(s1), median1 = m1,
        diff.median = m1 - m0, rel.diff = (m1 - m0) / abs(m0), lfc.abs = log2(abs(m1) / abs(m0)), lfc = log2(m1 / m0),
        r = wilcox.r, pval = wilcox.p
      )
    }
  }
  # This part try to check if any empty vector is empty or equal length
  # The warning is fine to skip, as lapply keep track of the vector
  res <- Reduce(rbind, furrr::future_map(colnames(mat0), function(.) {
    if ((length(mat0[[.]]) != length(mat1[[.]])) | any(is.na(mat0[[.]])) | any(is.na(mat1[[.]]))) {
      ## m0 <- median(mat0[[.]])
      ## m1 <- median(mat1[[.]])
      m0 <- median(mat0[[.]], na.rm = T)
      m1 <- median(mat1[[.]], na.rm = T)
      data.frame(
        rxn = .,
        lb0 = min(mat0[[.]]), ub0 = max(mat0[[.]]), median0 = m0,
        lb1 = min(mat1[[.]]), ub1 = max(mat1[[.]]), median1 = m1,
        diff.median = m1 - m0, rel.diff = (m1 - m0) / abs(m0), lfc.abs = log2(abs(m1) / abs(m0)), lfc = log2(m1 / m0),
        r = NA, pval = NA
      )
    } else {
      dflux.test(mat0[[.]], mat1[[.]])
    }
  }))
  res %>%
      mutate(padj = p.adjust(pval, method="fdr")) %>%
      filter(padj < padj.cutoff)
  ## %>%
  ##     mutate(dir =ifelse(!(padj<padj.cutoff & abs(r)>r.cutoff & abs(diff.median)>df.cutoff & abs(rel.diff)>rdf.cutoff),
  ##                        0,
  ##                 ifelse(diff.median>0, 1, -1))
             ##)
}


# Split by tissue & subsplit by group (genotype&diet)
# NOTE Comparison of Flux Ob.HFD-Ob.ND & WTHFD-WTND
grid <- expand.grid(tissue = unique(dat$tissue),
            group = combn(unique(dat$group), 2 , FUN = paste0, collapse='-')) %>%
    separate(group, into = c("c1", "c0"), sep = "-") %>%
    filter((c1 == "WT HFD" & c0 == "WT ND") | (c1 == "Ob/Ob HFD" & c0 == "Ob/Ob ND"))


unique(paste(grid$c1, grid$c0)) # TEST Correct Comparsion

fluxs <- split(dat, paste(dat$tissue, dat$group))
rm(dat)

multisession(6)
options(future.globals.maxSize = 8000 * 1024^2)

fluxes <- furrr::future_map(fluxs, function(x) {
    y <- x %>%
        select_if(is.numeric) %>%
        #dplyr::select(-c("group", "tissue", "source", "...1")) %>%
        dplyr::select(-c("...1")) %>%
        as.data.frame()
    y[is.na(y)] <- 0
    y
})

df.res <- furrr::future_map(1:nrow(grid), function(x){
                        df.wilcox(fluxes[[paste(grid[x,1], grid[x,3])]], # TEST ND
                                  fluxes[[paste(grid[x,1], grid[x,2])]], # TEST HFD
                                  ## by = 2,
                                  ## nc = 3,
                                  padj.cutoff = 0.05,
                                  r.cutoff = Inf,
                                  df.cutoff = 0, rdf.cutoff = 0
                                  )
                    }
                    )
names(df.res) <- paste0(grid[,1], ':',grid[,2], 'vs' ,grid[,3])
save(df.res, grid , file = "./tmp/tmpFluxes_ACHR2_extended.Rdata")

#################################################################################
load("./tmp/tmpFluxes_ACHR2_extended.Rdata")

df.res.sig <- lapply(df.res, function(x) x %>%
                                           filter(padj < 0.001 & abs(lfc) > 0.7)
                       )

####################################################################################
# Effect Size
df.res.sig.df <- lapply(df.res.sig, function(x) summarise(x, Up = sum(diff.median > 1) ,
                                                          Down = sum(diff.median < 1) ))

## stack(df.res.names)
df.res.sig.df <- enframe(df.res.sig.df) %>% unnest()  %>%
    mutate(Tissue = str_split(name, ":", simplify = T)[,1]) %>%
    mutate(Comparison = str_split(name, ":", simplify = T)[,2])

png("./paper/figures/DF.effectsizes.png", units = "in",  width = 10, height = 4, res=150)
df.res.sig.df  %>%
    ggplot() +
    geom_bar(aes(y = Tissue, x = -Down),
             stat = "identity", position = "dodge", fill = "#F15412", width = 1) +
    geom_bar(aes(y = Tissue, x = Up),
             stat = "identity", position = "dodge", fill = "#34B3F1", width = 1) +
    scale_x_continuous(labels = scales::label_number(suffix = " K", scale = 1e-3)) +
    ggh4x::facet_nested(Tissue~Comparison, scale = "free", switch = "y") +
    labs(x = "Differential Gene Expression Contrasts Up/Down Counts", y = "Tissues") +
    theme_bw() +
    theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          axis.text.y = element_blank(),
          axis.title.x = element_text(color = "grey20", size = 12, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          strip.background = element_rect(fill="#EEEEEE"),
          strip.text.y.left = element_text(angle = 0))
dev.off()


####################################################################################
## Overlap analaysis
df.res.names <- lapply(df.res.sig, function(x) x %>% .[["rxn"]])

uniquecontrasts <- unique(paste0(grid$c1, 'vs', grid$c0))

###


names(df.res.names) <- names(df.res.names) %>% make.names()
wt <- grep(pattern = make.names(uniquecontrasts[1]), names(df.res.names), value = T)
ob <- grep(make.names(uniquecontrasts[2]), names(df.res.names), value = T)
all <- c(ob, wt)

png(file = "./paper/figures/df.wt.overlap.png", units =  "in", width = 8, height = 5.5, res = 150)
wt.df <- UpSetR::fromList(df.res.names[wt])
colnames(wt.df) <- colnames(wt.df) %>% stringr::str_split(pattern = "\\.", simplify = T) %>% .[,1]
upset(wt.df, colnames(wt.df), name='Geneset',
                 sort_intersections_by = "degree",
                 n_intersections = 8,
      min_degree = 3,
                 width_ratio=0.2,
                 height_ratio = 0.5,
                 queries = list(upset_query(intersect = colnames(wt.df),
                                            color = "#F15412", fill = "#F15412")),
                    themes=upset_modify_themes(
                        list(
                            'Intersection size'=theme(
                                axis.text=element_text(size=12, face='bold'),
                                axis.title=element_text(size=12)
                            ),
                            'intersections_matrix'=theme(
                                axis.text=element_text(size=12, face='bold'),
                                axis.title=element_text(size=15)
                            )
                        )
                    ))
dev.off()


png(file = "./paper/figures/df.ob.overlap.png", units =  "in", width = 8, height = 5.5, res = 150)
ob.df <- UpSetR::fromList(df.res.names[ob])
colnames(ob.df) <- colnames(ob.df) %>% stringr::str_split(pattern = "\\.", simplify = T) %>% .[,1]
upset(ob.df,
      colnames(ob.df), name='geneset',
      sort_intersections_by = "degree",
      n_intersections = 8,
      min_degree = 3,
      width_ratio=0.2,
      height_ratio = 0.5,
      queries = list(upset_query(intersect = colnames(ob.df),
                                 color = "#F15412", fill =  "#F15412")),
      themes=upset_modify_themes(
          list(
              'Intersection size'=theme(
                  axis.text=element_text(size=12, face='bold'),
                  axis.title=element_text(size=12)
              ),
              'intersections_matrix'=theme(
                  axis.text=element_text(size=12, face='bold'),
                  axis.title=element_text(size=15)
              )
          )
      ))
dev.off()

png(file = "./paper/figures/df.wilcox.all.overlap.png", units =  "in", width = 10.5, height = 8.5, res = 150)
all.df <- UpSetR::fromList(df.res.names)
colnames(all.df) <- gsub(x = colnames(all.df), pattern = "Ob.Ob" , replacement = "Ob/Ob") %>%
    str_extract_all("^.*?(WT|Ob/Ob)", simplify = T) %>% paste(., "Context")

upset(all.df, colnames(all.df),
      name='Geneset',
      sort_intersections_by = "degree",
      min_degree = 12,
      max_degree = 14,
      n_intersections = 8,
      width_ratio=0.2,
      ## intersections = "all",
      height_ratio = 0.5,
      ## max_combinations_datapoints_n = 100^100 ,
      queries = list(upset_query(intersect = colnames(all.df),
                          color = "#F15412", fill = "#F15412")),
      themes=upset_modify_themes(
          list(
              'Intersection size'=theme(
                  axis.text=element_text(size=12, face='bold'),
                  axis.title=element_text(size=12)
              ),
              'intersections_matrix'=theme(
                  axis.text=element_text(size=12, face='bold'),
                  axis.title=element_text(size=15)
              )
          )
      )
      )
dev.off()

## png("./paper/figures/overlap.wilcox.png", width = 11, height = 12, units = "in", res = 150)
## (wt.plot+ob.plot)/ all.plot + plot_annotation(tag_levels = 'A')
## dev.off()

Reduce(intersect, df.res.names) ## Genes DE in all sets
sapply(df.wilcox.res.sig, length) # 7 tissues

# Load Reaction Names
humangem.df <- readxl::read_xlsx("./models/Human-GEM.xlsx")

## wilRes.diet <- humangem.df %>%
##     filter(ID %in% Reduce(intersect, df.res.names[grepl(paste(uniquecontrasts[c(1,4)], collapse="|"), names(df.res.names))])) %>% ## Genes DE in all sets
##     dplyr::select(ID, EQUATION, SUBSYSTEM, 'GENE ASSOCIATION') %>%
##     mutate(Method = "Wilcox-Rank", Interaction = "Diet")
## wilRes.genotype <- humangem.df %>%
##     filter(ID %in% Reduce(intersect, df.res.names[grepl(paste(uniquecontrasts[c(2,3)], collapse="|"), names(df.res.names))])) %>% ## Genes DE in all sets
##     dplyr::select(ID, EQUATION, SUBSYSTEM, 'GENE ASSOCIATION') %>%
##     mutate(Method = "Wilcox-Rank", Interaction = "Genotype")

wilRes <- humangem.df %>%
    filter(ID %in% Reduce(intersect, df.res.names)) %>% ## Genes DE in all sets
    dplyr::select(ID, EQUATION, SUBSYSTEM, 'GENE ASSOCIATION') %>%
    separate_rows(4 , sep = "or") %>%
    mutate(GeneAssociation = trimws(`GENE ASSOCIATION`)) %>%
    dplyr::select(-4) %>%
    mutate(GeneAssociation = mapIds(org.Hs.eg.db,
                                 column = "SYMBOL",
                                 keys = wilRes$GeneAssociation,
                      keytype = "ENSEMBL")
                      ) %>%
        with_groups(ID, mutate, GeneAssociation = paste0(GeneAssociation, collapse = " or ")) %>% distinct()
save(wilRes, file = "final_wil_res.Rdata")
######################################################################
load(file = "final_wil_res.Rdata")

unique(wilRes$SUBSYSTEM)


# Rx with GENES only
wilRes <- wilRes %>%
    ## filter(!(SUBSYSTEM %in% c("Exchange/demand reactions", "Transport reactions"))) %>%
    dplyr::select(SUBSYSTEM,  EQUATION, GeneAssociation) %>% #EQUATION
    separate_rows(GeneAssociation , sep = "or") %>%
    mutate(GeneAssociation = trimws(GeneAssociation)) %>% unique() %>%
    arrange(SUBSYSTEM, GeneAssociation) %>%
    group_by(SUBSYSTEM, EQUATION) %>%
    summarise(across(everything(), ~toString(na.omit(unique(.))))) %>% unique() %>%
    flextable::as_grouped_data(groups = c("SUBSYSTEM", "EQUATION")) %>%
    mutate(GeneAssociation = gsub( x= GeneAssociation , "NULL", replacement = "") ) %>%
    dplyr::rename(`Associated Genes` = GeneAssociation) %>%
    dplyr::rename( `Metabolic Pathway` = SUBSYSTEM )

wilRes %>% gt::gt() %>%
    gt::gtsave("paper/figures/tableRxWilcoxTrimmed.tex")

wilRes %>%
    ## dplyr::rename( Contrast = Interaction )  %>%
    flextable::flextable() %>%
    flextable::width(j = 1:3, width = c(2,5,2), unit = "in") %>%
    flextable::fontsize(size = 13, j =1) %>%
    ## flextable::theme_booktabs() %>%
    ## print(., preview = "pdf")
    ## flextable::autofit() %>%
    flextable::save_as_image(path = "./paper/figures/tableRxWilcoxTrimmed.png", webshot = "webshot2")


# Full Table
library(openxlsx)
wb <- createWorkbook()
for (i in 1:length(names(df.res))) {
    message("Writing sheet", i)

    addWorksheet(wb, sheetName = make.names(names(df.res)[i]))

    left_join(df.res[[i]],
        humangem.df %>% dplyr::select(ID, EQUATION, SUBSYSTEM, "GENE ASSOCIATION"),
        by = c("rxn" = "ID")
    ) %>%
        arrange(SUBSYSTEM) %>%
        flextable::as_grouped_data(groups = c("SUBSYSTEM")) %>%
        unique() %>% writeDataTable(wb, sheet = i, .)
}
saveWorkbook(wb, "./paper/supp/dflux.xlsx", overwrite = TRUE) ## save to working directory

## df.wilcox.res <- future_map(names(fluxes), function(xtissue) {
##            out <- future_map(fluxes[[xtissue]], function(xcol){
##                              pairwise.wilcox.test(xcol,
##                                                   fluxs[[xtissue]]$group,
##                                                   p.adjust.method = "bonferroni")})
##            names(out) <- names(fluxes[[xtissue]])
##            x <- sapply(out, function(y) { # just to make sure nothing random happens
##                p <- y$p.value
##                n <- outer(rownames(p), colnames(p), paste, sep=' V ')
##                p <- as.vector(p)
##                names(p) <- n
##                p
##         })
##            x
## })
## names(df.wilcox.res) <- names(fluxes)

## ## save(df.wilcox.res, file = "./tmp/tmpFluxes_ACHR2.Rdata") # ACHR
## #################################################################################
## load(file = "./tmp/tmpFluxes_ACHR2.Rdata") # ACHR

sapply(df.wilcox.res, length) # 7 tissues

## View(head(df.wilcox.res$Ao[-c(4,7,8),]))
## rownames(df.wilcox.res$Ao)[-c(4, 7, 8)]
## df.wilcox.res.sig <- furrr::future_map(df.wilcox.res, function(x) x %>%
##                                   as.data.frame %>%
##                                   .[-c(4, 7, 8),] %>%
##                                   select_if(~all(. < 10^-10)) %>%
##                                   names(.)
##        )


#################################################################################
# If still need the mean values

## df.ND.res <- furrr::future_map(
##   fluxes,
##   ~ df.wilcox(.[["WT ND"]], .[["Ob/Ob ND"]],
##     ## by = 2,
##     ## nc = 3,
##     padj.cutoff = 0.05,
##     r.cutoff = Inf,
##     df.cutoff = 0, rdf.cutoff = 0
##   )
## )
## df.hfd.res <- furrr::future_map(
##   fluxes,
##   ~ df.wilcox(.[["WT HFD"]], .[["Ob/Ob HFD"]],
##     #by = 2,
##     padj.cutoff = 0.05,
##     r.cutoff = Inf,
##     df.cutoff = 0, rdf.cutoff = 0
##   )
## )

## ## Merge, select rxn, diff.median, pval, r and add comparison name
## n <- names(df.ND.res)
## df.ND.res2 <- map(
##   n,
##   ~ df.ND.res[[.]] %>%
##     ## dplyr::select(c(rxn, lb0, ub0, median0, median1, lb1, ub1, diff.median, r, pval, lfc.abs)) %>%
##     mutate(comparison = "ND") %>%
##     mutate(tissue = .x)
## )
## df.ND.res2 <- Reduce(rbind, df.ND.res2)
## df.hfd.res2 <- map(
##   n,
##   ~ df.hfd.res[[.]] %>%
##     ## dplyr::select(c(rxn, , lb0, ub0, median0, median1, lb1, ub1, diff.median, r, pval, lfc.abs)) %>%
##     mutate(comparison = "HFD") %>%
##     mutate(tissue = .x)
## )
## df.hfd.res2 <- Reduce(rbind, df.hfd.res2)
## df.res <- rbind(df.ND.res2, df.hfd.res2) %>%
##     arrange(pval)

## head(df.res)

## ## Filter by diff.median > 0 or top 250 reactions (more managable sets!)
## df.res.filtered <- df.res %>%
##   filter(padj < 0.001) %>%
##   filter(abs(r) > 0.7 | is.na(r))  # %>%
##   ## group_by(comparison, tissue) %>%
##   ## top_n(abs(rel.diff), n = 1000)  %>%
##   ## ungroup()
## dim(df.res.filtered)

## x <- split(df.res.filtered, df.res.filtered$comparison)
## x <- map(x, ~ split(., .$tissue))
## de.res.list.ND <- map(x[["ND"]], ~ .[["rxn"]])
## names(de.res.list.ND) <- paste0(names(de.res.list.ND), ".ND")
## de.res.list.HFD <- map(x[["HFD"]], ~ .[["rxn"]])
## names(de.res.list.HFD) <- paste0(names(de.res.list.HFD), ".HFD")
## de.res.list <- c(de.res.list.HFD, de.res.list.ND)


## ## dev.new(width = 8, height = 9, res = 200)

## png(file = "./paper/figures/upset.ND.png", width = 8, height = 9, res = 200, units = "in")
## UpSetR::upset(fromList(de.res.list.ND), order.by = "degree", nsets = 7, nintersects = 10)
## dev.off()

## ## combos <- Reduce(c, lapply(2:length(de.res.list.ND),
## ##             function(x) combn(1:length(de.res.list.ND),x,simplify=FALSE)))
## ## # 4-ways is between 91 and 57
## ## nd.unique <- lapply(combos, function(x) Reduce(intersect, de.res.list.ND[x]))[57:91] %>%
## ##     unlist %>% unique #%>% length

## png(file = "./paper/figures/upset.HFD.png", width = 8, height = 9, res = 200, units = "in")
## UpSetR::upset(fromList(de.res.list.HFD), order.by = "degree", nsets = 7, nintersects = 10)
## dev.off()
## ## combos <- Reduce(c, lapply(2:length(de.res.list.HFD),
## ##             function(x) combn(1:length(de.res.list.HFD),x,simplify=FALSE)))
## ## # 4-ways is between 91 and 57
## ## hfd.unique <- lapply(combos, function(x) Reduce(intersect, de.res.list.HFD[x]))[57:91] %>%
## ##     unlist %>% unique #%>% length
## ## intersect(hfd.unique, nd.unique)

## png(file = "./paper/figures/upset.reactions.both.ND.HFD.png", width = 8, height = 9, res = 200, units = "in")
## UpSetR::upset(fromList(de.res.list), order.by = "degree", nsets = 14, nintersects = 20)
## dev.off()

## de.common.ND <- Reduce(intersect, de.res.list.ND) ## Genes DE in all sets
## de.common.HFD <- Reduce(intersect, de.res.list.HFD) ## Genes DE in all sets

#save(df.ND.res, df.hfd.res, file = "tmpFluxes_ACHR.Rdata") # ACHR


## https://stackoverflow.com/questions/56361724/multiple-tests-with-pairwise-combinations-using-dplyr-tidyverse

## dat <- dat %>%
##     pivot_longer(
##     -c(`...1`, source, tissue, group),
##     names_to = "Reaction"
##   )
## combos <- unique(dat$group) %>%
##   combn(2, simplify = F) %>%
##   set_names(map_chr(., ~ paste(., collapse = "_")))

##  dat2 %>%
##      group_split(tissue) %>%
##   set_names(map_chr(., ~ unique(.$tissue))) %>% # add tissue name to list items
##   map_df(function(x) {
##     map_df(combos, function(y) {
##       filter(x, group %in% y) %>%
##         wilcox.test(count ~ group , data = .) %>%
##         broom::tidy()
##     }, .id = "contrast")
##   }, .id = "trial")


## test.fun <- function(dat, col) {

## ## https://stackoverflow.com/questions/21271449/how-to-apply-the-wilcox-test-to-a-whole-dataframe-in-r
##  c1 <- combn(unique(dat$group),2)
##  sigs <- list()
##  for(i in 1:ncol(c1)) {
##     sigs[[i]] <- wilcox.test(
##                    dat[dat$group == c1[1,i],col],
##                    dat[dat$group == c1[2,i],col]
##                  )
##     }
##     names(sigs) <- paste("Group",c1[1,],"by Group",c1[2,])

##  tests <- data.frame(Test=names(sigs),
##                     W=unlist(lapply(sigs,function(x) x$statistic)),
##                     p=unlist(lapply(sigs,function(x) x$p.value)),row.names=NULL)

##  return(tests)
## }


## tests <- lapply(colnames(dat)[-1],function(x) test.fun(dat,x))
## names(tests) <- colnames(dat)[-1]
## # tests <- do.call(rbind, tests) reprints as data.frame

## # This solution is not "slow" and outperforms the other answers significantly:
## system.time(
##   rep(
##    tests <- lapply(colnames(dat)[-1],function(x) test.fun(dat,x)),10000
##   )
## )
## group <- dat[,9100:9207]$group
## tissue <- dat[,9100:9207]$tissue
## x <- t(dat[, 9100:9207] %>% select_if(is.numeric))
## pwt <- scran::pairwiseWilcox(x, groups = group, block = tissue)

## x <-
##     wilcox.paired.multcomp(MAR03802 ~ group | tissue, data = dat)

