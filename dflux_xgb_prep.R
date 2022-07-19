# Run the recipe, workflow and hypertuning results
library(tidymodels)
library(doParallel)
set.seed(123)

load("samplesAllTidy_ACHR2.Rdata")
transportRxns <- readRDS(file = "./tmp/transportRxns.RDS")

not_all_na <- function(x) any(!is.na(x))
dat <-
  dat %>%
  ## mutate(factor = factor(paste0(tissue, group))) %>%
  ## select(-c(tissue, group, source, `...1`)) %>%
  select(-c( source, `...1`)) %>%
  select(-any_of(transportRxns)) %>%
  select(where(not_all_na)) # %>%
  ## .[,1500:ncol(.)] %>%  # FIXME to be removed
  ## with_groups(c(tissue, group), sample_n, size = 100) # FIXME to be removed
 # mutate_if(is.numeric, ~ replace(., is.na(.), 0)) # Rstats needs missing to be NA, but fluxe of zero is also a possibility

# 1. Data Manipulation
## 1.1 Data splitting
vb_split <- initial_split(dat, prop = 0.8, strata = tissue)
vb_split


# extract training and testing sets
vb_train <- training(vb_split)
vb_test <- testing(vb_split)

rm(dat); gc()

## 1.2 Setting Recipe
vb_recipe <-
  recipe(group~., data = vb_train) %>%
  #step_impute_knn(all_numeric()) %>%
  #step_impute_median(all_numeric()) %>%
  step_nzv(all_numeric()) %>%
  step_normalize(all_numeric()) %>%
      ###
  ## step_unknown(all_nominal_predictors()) %>%
  step_dummy(all_nominal_predictors())
vb_recipe

## vb_prep <- prep(vb_recipe)
## summary(vb_prep)

## train_trans <- juice(vb_prep)
## test_trans <- bake(vb_prep, new_data = vb_test)

xgb_spec <- parsnip::boost_tree(
  trees = 200,
  tree_depth = tune(),
  min_n = tune(),
  sample_size = tune(),
  mtry = tune(),         ## randomness
  loss_reduction = tune(),                     ## first three: model complexity
  learn_rate = tune(),                         ## step size
) %>% set_engine("xgboost", num_class = 4,
                 objective = "multi:softprob",
                 tree_method = 'gpu_hist',
                 verbose=1) %>%
    set_mode("classification")



## 2.2 set Grid Search
xgb_grid <- grid_latin_hypercube(
  min_n(), #2nd
  tree_depth(), #3d
  loss_reduction(),
  sample_size = sample_prop(),
  finalize(mtry(), vb_train),
  learn_rate(),
  size = 10
)
xgb_grid



## 2.3 Model Workflow
xgb_wf <- workflow() %>%
    ## add_formula(group~.) %>%
  add_recipe(vb_recipe) %>%
  add_model(xgb_spec)
xgb_wf


## 2.4 Cross validation
vb_folds <- vfold_cv(vb_train,
                     v = 5,
                     strata = tissue)
vb_folds

#################################
#################################
# Hypertuning
doParallel::registerDoParallel(2)


# TODO bake first better, no need for recipe in wf then!

xgb_res <- tune_grid(
    xgb_wf,
    resamples = vb_folds,
    grid = xgb_grid,
    control = control_grid(save_pred = TRUE, # To obtain ROC curve
                         verbose = TRUE)
)
#xgb_res

saveRDS(xgb_res, file = "tmp/xgb_res.RDS")
saveRDS(xgb_wf, file = "tmp/xgb_wf.RDS")
saveRDS(vb_split, file = "tmp/vb_split.RDS")
