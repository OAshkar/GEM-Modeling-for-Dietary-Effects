library(tidyverse)
library(ggh4x)
library(SummarizedExperiment)

dat <- read.csv("simulations/original_EE.csv")

dat <- dplyr::select(dat, -X)
rownames(dat) <- as.character(seq(0.1, 0.8, 0.1))
dat <- rownames_to_column(dat, "medium")
x <- colnames(dat)
dat <- tidyr::pivot_longer(dat, -medium, names_to = "sample") %>%
  mutate(genotype = gsub(
    sample,
    pattern = "(.*)_(WT|ob_ob)_(HFD|ND)_(Ao|Ep|He|Sk|Li|Hy|Hi).*",
    replacement = "\\2",
    perl = t
  )) %>%
  mutate(diet = gsub(
    sample,
    pattern = "(.*)_(WT|ob_ob)_(HFD|ND)_(Ao|Ep|He|Sk|Li|Hy|Hi).*",
    replacement = "\\3",
    perl = t
  )) %>%
  mutate(tissue = gsub(
    sample,
    pattern = "(.*)_(WT|ob_ob)_(HFD|ND)_(Ao|Ep|He|Sk|Li|Hy|Hi).*",
    replacement = "\\4",
    perl = t
  )) %>%
  mutate(dietgenotype = paste(diet, genotype)) %>%
  mutate(medium = as.numeric(medium))

ggplot(dat, aes(x = medium, y = value)) +
  facet_wrap2(vars(tissue), nrow = 3, ncol = 3) +
  geom_point(aes(color = dietgenotype, shape = dietgenotype, size = 2)) +
  geom_smooth(aes(color = dietgenotype), se = F)


ggplot(dat, aes(x = medium, y = value)) +
  facet_grid(tissue ~ genotype + diet) +
  geom_point(aes(color = dietgenotype, shape = dietgenotype)) +
  geom_smooth(aes(color = dietgenotype), se = F)

dev.new(height = 6, width = 7, res = 150)
seq_n <- seq(0.1, 0.9, by = 0.1)
seq_labels <- paste0(seq_n * 100, "%")
ggplot(dat, aes(x = medium, y = value)) +
  facet_grid(tissue ~ diet) +
  geom_point(aes(color = genotype, shape = dietgenotype)) +
  geom_smooth(aes(color = genotype), se = F) +
  xlab("Percentage of Carbohydrates") +
  ylab("ATPc Flux Rate (mmol/hr/gmDW)") +
  scale_x_continuous(breaks = seq_n, labels = seq_labels) +
  theme(
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10, angle = 45),
  )

png(file = "paper/figure/ATPflux0.1Carbs.png", unit = "in", res = 150, height = 4, width = 7)
ggplot(dat %>% filter(medium == 0.1), aes(x = sample, y = value, color = diet)) +
  geom_bar(stat = "identity") +
  facet_grid(tissue ~ genotype, scales = "free_x") +
  xlab("") +
  ylab("ATPc Flux Rate (mmol/hr/gmDW)") +
  # scale_x_continuous(breaks = seq_n, labels = seq_labels) +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10, angle = 45),
    axis.text.x = element_blank()
  ) +
  geom_hline(yintercept = 48, linetype = "dashed", color = "red", size = 0.2)
dev.off()

# Get the median and diff of fluxes

png(file = "paper/figure/diffATPflux0.1Carbs.png", unit = "in", res = 150, height = 4, width = 7)
dat %>%
  filter(medium == 0.1) %>%
  group_by(tissue, genotype, diet) %>%
  summarise(median = median(value)) %>%
  ungroup() %>%
  group_by(tissue, diet) %>%
  summarize(diff = diff(median)) %>%
  ungroup() %>%
  ggplot(aes(x = tissue, y = diff, color = diet)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(
    strip.text.x = element_text(size = 10)
  ) +
  theme(
    strip.text.y = element_text(size = 10, angle = 45)
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_grid(tissue ~ diet, scales = "free_x") +
  ylab("Difference of ATPc Flux Rate Between Same Genotype (mmol/hr/gmDW)")
dev.off()


### Flux Sampling
temp <- list.files(path = "samplingRes", pattern = "*.csv", full.names = T)
myfiles <- purrr::map_df(temp, read_csv)
myfiles <- sjmisc::rotate_df(myfiles)


colnames(myfiles) <- paste0(rep(basename(temp), each = 40), "_", 1:40)
samples <- data.frame(samples = colnames(myfiles)) %>%
  mutate(Genotype = ifelse(grepl(pattern = ".*WT.*", x = samples), "WT", "Ob")) %>%
  mutate(Diet = ifelse(grepl(pattern = ".*ND.*", x = samples), "ND", "HFD")) %>%
  mutate(tissue = case_when(
    grepl(pattern = ".*Ao.*", x = samples) ~ "Ao",
    grepl(pattern = ".*Ep.*", x = samples) ~ "Ep",
    grepl(pattern = ".*He.*", x = samples) ~ "He",
    grepl(pattern = ".*Hy.*", x = samples) ~ "Hy",
    grepl(pattern = ".*Hi.*", x = samples) ~ "Hi",
    grepl(pattern = ".*Li.*", x = samples) ~ "Li",
    grepl(pattern = ".*Sk.*", x = samples) ~ "Sk"
  ))

se <- SummarizedExperiment(assay = myfiles, colData = samples)
library(edgeR)
plotDensities(log(assay(se)))

dge <- DGEList(counts = myfiles, samples = samples)


x <- calcNormFactors(x, method = "TMM")