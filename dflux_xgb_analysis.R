library(tidyverse)
library(tidymodels)
library(vip)

final_fit_xgb <- readRDS("tmp/final_fit_xgb.RDS")

#final_xgb %>%
#  fit(data = train_trans)


## 2.9 Evaluate Model
final_fit_xgb %>%
    collect_metrics() # NOTE No over-fit

final_res_pred <- final_fit_xgb %>%
    collect_predictions()
final_res_pred %>% View()

View(head(final_res_pred))
View(final_res_pred)

# Confusion Matrix
final_fit_xgb %>%
    collect_predictions() %>% conf_mat(group, .pred_class) %>%
    autoplot(type = "heatmap") + scale_fill_viridis_c()


### 2.9.2 AUC
auc <- final_res_pred %>%
  roc_auc(truth = group,
          4:7)

### 2.9.3 Accuracy
accu <- final_res_pred %>%
  metrics(truth = group,
          estimate = .pred_class)

### 2.9.4 Recall
recall <- final_res_pred %>%
  recall(truth = group,
          estimate = .pred_class)

### 2.9.5 Precision
precision <- final_res_pred %>%
  precision(truth = group,
          estimate = .pred_class)


### ROC curve
png(file = "./paper/figures/xgb.ROC.png", res = 100)
roc_curve(final_res_pred,
          truth = group,
            4:7) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity)) +
    geom_path(aes(color=.level, lty = .level)) +
    annotate(geom = "text", label = paste("ROC AUC =", auc$.estimate), x = 0.25, y  = 0.75)  +
    annotate(geom = "text", label = paste("Accuracy =", round(accu$.estimate[1], 6)), x = 0.25, y  = 0.65)  +
    annotate(geom = "text", label = paste("Kappa =", round(accu$.estimate[2],6)), x = 0.75, y  = 0.35)  +
    annotate(geom = "text", label = paste("Recall =", round(recall$.estimate,6)), x = 0.75, y  = 0.25)  +
    annotate(geom = "text", label = paste("Precision =", round(precision$.estimate,6)), x = 0.75, y  = 0.15)  +
    geom_abline(lty = 3) +
    coord_equal() +
    theme_bw() +
    labs(color='Class', lty = 'Class')
                                        #+ autoplot()
dev.off()

#######
#######
v <- final_fit_xgb %>%
  extract_fit_parsnip() %>%
    vip::vip(geom = "point")
png(file = "./paper/figures/xgb.top10.png")
v
dev.off()


# Load Reaction Names
humangem.df <- readxl::read_xlsx("./models/Human-GEM.xlsx")
# Tidy table with VIP, BIGG name (manual), formula

transportRxns <- humangem.df %>%
    filter(SUBSYSTEM %in% c("Transport reactions", "Artificail reactions")) %>%
    .$ID
saveRDS(transportRxns, file = "./tmp/transportRxns.RDS")


xgbRes <- humangem.df %>%
    filter(ID %in% v$data[[1]]) %>%
    dplyr::select(ID, EQUATION, SUBSYSTEM, 'GENE ASSOCIATION') %>%
    mutate(Method = "xgboost") %>%
    dplyr::slice(match(v$data[[1]], ID)) %>%
    mutate(Rank = 1:10) %>%
    separate_rows(4 , sep = "and") %>%
    mutate(GeneAssociation = trimws(`GENE ASSOCIATION`)) %>%
    dplyr::select(-4)

xgbRes <-
    xgbRes %>% mutate(GeneAssociation = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                 column = "SYMBOL",
                                 keys = xgbRes$GeneAssociation,
                      keytype = "ENSEMBL")
                      ) %>%
        with_groups(ID,
                    mutate, GeneAssociation = paste0(GeneAssociation, collapse = " and ")) %>%
    distinct()

xgbRes[which(xgbRes$GeneAssociation == "NULL"), "GeneAssociation"] <- ""

xgbRes %>%
  select(SUBSYSTEM, EQUATION ,GeneAssociation, Rank) %>%
    arrange(SUBSYSTEM, Rank, GeneAssociation) %>%
    group_by(SUBSYSTEM) %>%
    flextable::as_grouped_data("SUBSYSTEM") %>%
    rename(`Associated Genes` = GeneAssociation) %>%
    rename( `Metabolic Pathway` = SUBSYSTEM )  %>%
    gt::gt() %>%
    gt::fmt_missing(columns = everything(),
                rows = everything(),
                missing_text = "") %>%
    gt::gtsave("paper/figures/tableXgbVIP.tex")


## flextable::flextable() %>%
##     flextable::width(j = 1:3, width = c(2,6,2)) %>%
##     flextable::theme_booktabs() %>%
##     flextable::fontsize(size = 14, part = "header") %>%
##     flextable::fontsize(size = 12, part = "header") %>%
##     flextable::autofit() %>%
##     flextable::save_as_image(path = "./paper/figures/tableXgbVIP.png", webshot = "webshot2")



v2 <- final_fit_xgb %>%
    extract_fit_parsnip() %>%
    vip::vip(geom = "point", num_features = 200)
v2 <- data.frame(rx = v2$data[[1]], importance = v2$data$Importance)

v2 %>%
    left_join(humangem.df, by = c("rx" = "ID")) %>%
     select(rx, EQUATION, SUBSYSTEM, importance) %>%
    arrange(desc(importance)) %>%
    rename( `Reaction ID` = rx )  %>%
    xlsx::write.xlsx(file = "./paper/supp/xgb top 200 reactions.xlsx",
                     sheetName = "xgb results",
                     showNA = F, row.names = F)
