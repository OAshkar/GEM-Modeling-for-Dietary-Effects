library(tidymodels)
set.seed(123)

vb_split <- readRDS("./tmp/vb_split.RDS")
xgb_wf <- readRDS("./tmp/xgb_wf.RDS")
xgb_res <- readRDS("./tmp/xgb_res.RDS")

## 2.6 Explore Results
### 2.6.1 Metrics
collect_metrics(xgb_res) %>% gt::gt() %>%
    gt::gtsave("./Diagnostics/xgb_hypertune_metrics.tex")
collect_predictions(xgb_res)  %>% gt::gt() %>%
    gt::gtsave("./Diagnostics/xgb_hypertune_predications.tex") #%>% View()  # TODO Compare with final model to check overfitting

### 2.6.2 ROC
xgb_res_long <-
    xgb_res %>%
    collect_metrics() %>%
    filter(.metric == "roc_auc") %>%
    select(mean, mtry:sample_size) %>%
    pivot_longer(mtry:sample_size,
                 values_to = "value",
                 names_to = "parameter"
                 )

png(file = "./Diagnostics/AUC_hypertuning.png")
xgb_res_long %>%
  ggplot(aes(value, mean, color = parameter)) +
  geom_point(alpha = 0.8, show.legend = FALSE) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = NULL, y = "AUC")
dev.off()

## 2.7 Find best model
### 2.7.1 Show best parameters
show_best(xgb_res, "roc_auc") # %>% View()

### 2.7.2 Select best parameters
best_auc <- select_best(xgb_res, "roc_auc")
best_auc %>% gt::gt() %>% gt::gtsave("./Diagnostics/xgb_hypertune_best_auc.tex")

## 2.8 Finalize workflow
final_xgb <- finalize_workflow(
  xgb_wf,
  best_auc
)
final_xgb

rm(xgb_res, xgb_wf); gc()

### 2.8.1 Fit final model
 # fit the model on train and test on test by the last WF.
final_fit_xgb <- last_fit(final_xgb,
                          split = vb_split,
                          metrics = metric_set(
                              recall, precision, f_meas,
                              accuracy, kap, roc_auc, sens, spec)
                        )
saveRDS(final_fit_xgb, file = "tmp/final_fit_xgb.RDS")
