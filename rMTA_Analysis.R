library(tidyverse)
library(ComplexUpset)
library(ComplexHeatmap)
library(org.Mm.eg.db)

rmtaFiles <- list.files(path = "./rMTAResults", recursive = T, pattern = "rmtares.csv", full.names = T)
dat <- map_dfr(rmtaFiles, ~mutate(read_csv(.), source = . ))
# Sort by rTS or mTS
# Get top 500
# .[[detetedGenes]]
# Find most intersecting genes
dat <- dat %>% arrange(desc(rTS)) %>%
    mutate(class = gsub(x = source ,
                        pattern = "(.*rMTAResults/)(.*)(/rmtares..*)",
                        replacement = "\\2")) %>%
    mutate(genotype =  gsub(x= class,
                            pattern = "(ob_ob|WT)_(.*)_(.*)",
                            replacement = "\\1")) %>%
    mutate(diet =  gsub(x= class,
                            pattern = "(ob_ob|WT)_(.*)_(.*)",
                            replacement = "\\2")) %>%
    mutate(tissue =  gsub(x= class,
                            pattern = "(ob_ob|WT)_(.*)_(.*)",
                            replacement = "\\2"))

mtaList <- dat %>% with_groups(class, slice_max, order_by = rTS, n = 20) %>%
    dplyr::select(deletedGenes, class) %>% as.data.frame() %>%
    split(.$class) # group_split causes error
mtaList <- lapply(mtaList, function(x)x[["deletedGenes"]])

## names(mtaList) <- unique(dat$class)

mtaList

## mta.res <- list_to_matrix(mtaList) %>% as.data.frame()

## mta.comb <- make_comb_mat(mta.res, mode = "intersect", min_set_size = 10)
## mta.comb <- mta.comb[comb_degree(mta.comb) > 6]

pdf(file = "test.pdf")
UpSet(mta.comb)
dev.off()



UpSetR::fromList(mtaList) %>%  
  mutate_all(as.logical) -> x

colnames(x) %>%
    gsub(pattern =  "ob_ob", replacement  = "Ob/Ob") %>%
    gsub(pattern = "_", replacement = " ") -> colnames(x)

all_labs <- colnames(x)
Ob_labs <- all_labs[grep("Ob/Ob", all_labs)]
WT_labs <- all_labs[grep("WT", all_labs)]


png(file = "./paper/figures/rMTAUpset.png", width = 8, height = 11.5, res = 200, units = "in")
ComplexUpset::upset(x , colnames(x), name='Geneset',
                    sort_intersections_by = "degree", set_sizes = F,
                    intersections = "observed", # inclusive intersection,
                    width_ratio=0.2,
                    height_ratio = 0.5,
                    queries = list(
                        upset_query(intersect = Ob_labs,
                                    color = "purple", fill = "purple"),
                        upset_query(intersect = WT_labs,
                             color = "#F15412", fill = "#F15412"),
                        upset_query(intersect = all_labs,
                             color = "#34B3F1", fill = "#34B3F1")
                    ),
                    themes=upset_modify_themes(
                        list(
                            'Intersection size'=theme(
                                axis.text=element_text(size=20, face='bold'),
                                axis.title=element_text(size=20)
                            ),
                            'intersections_matrix'=theme(
                                axis.text=element_text(size=12, face='bold'),
                                axis.title=element_text(size=15)
                            )
                        )
                    ))
dev.off()

##########################
all_labs <- names(mtaList)
Ob_labs <- all_labs[grep("ob_ob", all_labs)]
WT_labs <- all_labs[grep("WT", all_labs)]

rmtaTopG.all <- AnnotationDbi::select(x = org.Mm.eg.db, keys  = Reduce(intersect, mtaList),
                      columns = c("GENENAME") , keytype =  "SYMBOL")

rmtaTopG.ob <- AnnotationDbi::select(x = org.Mm.eg.db, keys  = Reduce(intersect, mtaList[Ob_labs]),
                      columns = c("GENENAME") , keytype =  "SYMBOL")

rmtaTopG.wt <- AnnotationDbi::select(x = org.Mm.eg.db, keys  = Reduce(intersect, mtaList[WT_labs]),
                      columns = c("GENENAME") , keytype =  "SYMBOL")

intersect(rmtaTopG.all$SYMBOL, rmtaTopG.ob$SYMBOL)

rmtaTopG.all %>% gt::gt() %>%
    gt::gtsave("paper/figures/rMTATopG.tex")

flextable::flextable(rmtaTopG.all) %>%  flextable::theme_booktabs() %>%
    flextable::width(j = 1:2, width = c(1,3), unit = "in") %>%
    flextable::fontsize(size = 14, part = "header") %>%
    flextable::fontsize(size = 12, part = "header") %>%
    ## flextable::autofit() %>%
    flextable::save_as_image(x = .,
                  path = "paper/figures/rMTATopG.png",
                  webshot = "webshot2")
