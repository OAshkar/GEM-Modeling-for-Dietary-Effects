library(tidyverse)
library(reticulate)
setwd("~/Documents/Projects/ObAnalysis/")
load("samplesAllTidy_ACHR2.Rdata")
tissue <- dat$tissue
tissue <- case_when(
    tissue == "Ao" ~ "Aorta",
    tissue == "Ep" ~ "Epididymal Fat",
    tissue == "He" ~ "Heart",
    tissue == "Li" ~ "Liver",
    tissue == "Hy" ~ "Hypothalamus",
    tissue == "Hi" ~ "Hippocampus",
    tissue == "Sk" ~ "Skeletal Muscle",
    )


sublabs <- dat$group

mat <- dat %>%
    select(-"...1") %>% #write.csv(file = "tocluster.csv")
    select_if(is.numeric) %>%
    mutate_all(~replace(., is.na(.),0)) %>%
    t %>%
    M3C::featurefilter(percentile = 10, method = 'MAD') %>% .[["filtered_data"]] %>%
    scale #%>%
    #log1p
colMeans(mat) # zero
apply(mat, 2, sd) #1

all <- paste(labels, sublabs)
gc()

## pca <- M3C::pca(mat,labels = labels, dotsize = 3)
umap <- umap(mat,labels = labels,  dotsize = 3,
             n_neighbors = 200,
             min_dist = 0.5,
             metrci = "dice",
             seed = 333,
             min_dist = 0.1,
             method="umap-learn"
)
#save(x,y,filtered,labels, sublabs, file = "sampling_Clustering.Rdata")



# By Tissue
dat <- split(dat, tissue)

umap_plots <- list()
pca_plots <- list()
for(i in seq_along(names(dat))){
    print(i)
    labels <- dat[[i]]$group
    x <- dat[[i]] %>%
        select(-"...1") %>%
        select_if(is.numeric) %>%
        t %>%
        M3C::featurefilter(percentile = 20, method = 'MAD') %>% .[["filtered_data"]] %>%
        as.data.frame %>% mutate_all(function(y) bestNormalize::orderNorm(y) %>% .$x)
    pca_plots[[i]] <- M3C::pca(x,
             controlscale=T,
             scale=3,
             labels = labels,
             ## perplex = 3,
             ## legendtitle = paste(names(dat)[i], "Groups"),
             legendtitle = "Group",
             ## n_neighbors = 6,
             axistextsize = 10, legendtextsize = 11, dotsize = 3,
             textlabelsize = 3 #, seed = 333
             )
    umap_plots[[i]] <- umap(x,
             controlscale=T,
             scale=3,
             labels = labels,
             ## perplex = 3,
             ## legendtitle = paste(names(dat)[i], "Groups"),
             legendtitle = "Group",
             ggtitle = names(dat)[i],
             n_neighbors = 60,
             axistextsize = 10,
             legendtextsize = 11, dotsize = 3,
             textlabelsize = 3,
             seed = 1000,
             min_dist = 0.99,
             metric = "dice",
             method="umap-learn"
             )
}

umap_plots <- append(umap_plots, list(x=plot_spacer()), 6)
png(file = "paper/figures/samplingUmapTissuesACHR2.png", res = 150, width = 15, height = 10, units = "in")
patchwork::wrap_plots(umap_plots) +
    guide_area() +
    plot_layout(guides = "collect")
dev.off()

png(file = "paper/figures/samplingPCATissuesACHR2.png", res = 150, width = 15, height = 10, units = "in")
patchwork::wrap_plots(pca_plots)
dev.off()


patchwork::wrap_plots(tsne_plots)

#################################################################################
# Diagnostics Block
dat <- split(dat, dat$tissue)
labels <- dat[[i]]$group
x <- dat[[i]] %>%
    select(-"...1") %>%
    select_if(is.numeric) %>%
    t() %>%
    M3C::featurefilter(percentile = 10, method = "MAD") %>%
    .[["filtered_data"]] %>%
    scale()

neighbors <- seq(1,50,3) * 15
mindists <- c(0.001, 0.01, 0.05, 0.1, 0.5, 0.99)
dists <- c("euclidean", "manhattan", "canberra", "cosine", "hamming", "dice")

dotsize = 2
grid <-expand.grid(min_dist = mindists, meteric = dists, n_neighbors = neighbors)
# man, cranb, hamming causes warnings. First two for connectivity, last for inverse_transform not ava.
scaled <- list()
normed <- list()
groupes <- list()
for(i in seq_along(dat)){
    groupes[[i]] <- dat[[i]]$group
    tmp <- dat[[i]] %>%
        select(-"...1") %>%
        select_if(is.numeric) %>%
        t %>% M3C::featurefilter(percentile = 10, method = 'MAD') %>% .[["filtered_data"]]

    normed[[i]] <-   tmp  %>%
        as.data.frame %>% mutate_all(function(y) bestNormalize::orderNorm(y) %>% .$x)
    # regular scale
    scaled[[i]] <- tmp %>%
        scale
}
rm(dat)

pdf(file = "Diagnostics/umap_diagnostics196.pdf", width=20, height = 15)
for(i in 1:nrow(grid)){
    umaps <- list()
    umaps2 <- list()
    print(grid[i,])
    for(ii in seq_along(groupes)){
        labels <-  groupes[[ii]]
        x <- normed[[ii]]
        umaps[[ii]] <- umap(x,
                                controlscale=T,
                                scale=3,
                                labels = labels,
                                ## perplex = 3,
                                legendtitle = paste(grid[i,], "/" ,collapse = ''),
                                axistextsize = 10, legendtextsize = 11, dotsize = 3,
                                textlabelsize = 3,
                                seed = 1000,
                                n_neighbors =  grid$n_neighbors[i],
                                min_dist = grid$min_dist[i],
                                metric = grid$meteric[i],
                                method="umap-learn"
             )
        x <- scaled[[ii]]
        umaps2[[ii]] <- umap(x,
                        controlscale=T,
                        scale=3,
                        labels = labels,
                        ## perplex = 3,
                        legendtitle = paste("scale", grid[i,], "/" ,collapse = ''),
                        axistextsize = 10, legendtextsize = 11, dotsize = 3,
                        textlabelsize = 3,
                        seed = 1000,
                        n_neighbors =  grid$n_neighbors[i],
                        min_dist = grid$min_dist[i],
                        metric = grid$meteric[i],
                        method="umap-learn"
        )
}
print(patchwork::wrap_plots(umaps))
print(patchwork::wrap_plots(umaps2))
}
dev.off()



# PCA not bad on ACHR, but umap is very bad with ACHR. I will
# OPTG is horrible even with single samples


# Umap custom Function
umap <- function (mydata, labels = FALSE, printres = FALSE, seed = FALSE,
    axistextsize = 18, legendtextsize = 18, dotsize = 5, textlabelsize = 4,
    legendtitle = "Group", controlscale = FALSE, scale = 1, low = "grey",
    high = "red", colvec = c("skyblue", "gold", "violet", "darkorchid",
        "slateblue", "forestgreen", "violetred", "orange", "midnightblue",
        "grey31", "black"), printheight = 15, printwidth = 22,
    text = FALSE, ggtitle = NULL, ...)
{
    if (controlscale == TRUE && class(labels) %in% c("character",
        "factor") && scale %in% c(1, 2)) {
        stop("when categorical labels, use scale=3")
    }
    if (controlscale == TRUE && class(labels) %in% c("numeric") &&
        scale %in% c(3)) {
        stop("when continuous labels, use scale=1 or scale=2")
    }
    if (controlscale == FALSE && scale %in% c(2, 3)) {
        warning("if your trying to control the scale, please set controlscale=TRUE")
    }
    if (sum(is.na(labels)) > 0 && class(labels) %in% c("character",
        "factor")) {
        warning("there is NA values in the labels vector, setting to unknown")
        labels <- as.character(labels)
        labels[is.na(labels)] <- "Unknown"
    }
    if (sum(is.na(text)) > 0 && class(text) %in% c("character",
        "factor")) {
        warning("there is NA values in the text vector, setting to unknown")
        text <- as.character(text)
        text[is.na(text)] <- "Unknown"
    }
    message("***UMAP wrapper function***")
    message("running...")
    if (seed != FALSE) {
        set.seed(seed)
    }
    if (labels[1] == FALSE && text[1] == FALSE) {
        umap <- umap::umap(t(as.matrix(mydata)), ...)
        scores <- data.frame(umap$layout)
        p <- ggplot(data = scores, aes(x = X1, y = X2)) + geom_point(colour = "skyblue",
            size = dotsize) + theme_bw() + theme(legend.position = "none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = "black"),
            axis.text.x = element_text(size = axistextsize, colour = "black"),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
            scale_colour_manual(values = colvec)+ ggtitle(ggtitle)  + theme(plot.title = element_text(hjust = 0.5)) +
                    theme(plot.title = element_text(size = 15, face = "bold"))
        if (printres == TRUE) {
            message("printing UMAP to current directory...")
            png("UMAP.png", height = printheight, width = printwidth,
                units = "cm", res = 900, type = "cairo")
            print(p)
            dev.off()
        }
    }
    else if (labels[1] != FALSE && text[1] == FALSE) {
        umap <- umap::umap(t(as.matrix(mydata)), ...)
        scores <- data.frame(umap$layout)
        if (controlscale == TRUE) {
            if (scale == 1) {
                p <- ggplot(data = scores, aes(x = X1, y = X2)) +
                  geom_point(aes(colour = labels), size = dotsize) +
                  theme_bw() + theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                    colour = "black"), axis.text.x = element_text(size = axistextsize,
                    colour = "black"), axis.title.x = element_text(size = axistextsize),
                  axis.title.y = element_text(size = axistextsize),
                  legend.title = element_text(size = legendtextsize),
                  legend.text = element_text(size = legendtextsize)) +
                  labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral")+ ggtitle(ggtitle) + theme(plot.title = element_text(hjust = 0.5)) +
                    theme(plot.title = element_text(size = 15, face = "bold"))
            }
            else if (scale == 2) {
                p <- ggplot(data = scores, aes(x = X1, y = X2)) +
                  geom_point(aes(colour = labels), size = dotsize, ) +
                  theme_bw() + theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                    colour = "black"), axis.text.x = element_text(size = axistextsize,
                    colour = "black"), axis.title.x = element_text(size = axistextsize),
                  axis.title.y = element_text(size = axistextsize),
                  legend.title = element_text(size = legendtextsize),
                  legend.text = element_text(size = legendtextsize)) +
                  labs(colour = legendtitle) + scale_colour_gradient(low = low,
                  high = high)+ ggtitle(ggtitle) + theme(plot.title = element_text(hjust = 0.5)) +
                    theme(plot.title = element_text(size = 15, face = "bold"))
            }
            else if (scale == 3) {
                p <- ggplot(data = scores, aes(x = X1, y = X2)) +
                  geom_point(aes(colour = labels), size = dotsize, ) +
                  theme_bw() + theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                    colour = "black"), axis.text.x = element_text(size = axistextsize,
                    colour = "black"), axis.title.x = element_text(size = axistextsize),
                  axis.title.y = element_text(size = axistextsize),
                  legend.title = element_text(size = legendtextsize),
                  legend.text = element_text(size = legendtextsize)) +
                  labs(colour = legendtitle) + scale_colour_manual(values = colvec)+ ggtitle(ggtitle)  + theme(plot.title = element_text(hjust = 0.5)) +
                    theme(plot.title = element_text(size = 15, face = "bold"))
            }
        }
        else {
            p <- ggplot(data = scores, aes(x = X1, y = X2)) +
                geom_point(aes(colour = labels), size = dotsize) +
                theme_bw() + theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                  colour = "black"), axis.text.x = element_text(size = axistextsize,
                  colour = "black"), axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
                labs(colour = legendtitle)+ ggtitle(ggtitle)  + theme(plot.title = element_text(hjust = 0.5)) +
                    theme(plot.title = element_text(size = 15, face = "bold"))
        }
        if (printres == TRUE) {
            message("printing UMAP to current directory...")
            png("UMAPlabeled.png", height = printheight, width = printwidth,
                units = "cm", res = 900, type = "cairo")
            print(p)
            dev.off()
        }
    }
    else if (labels[1] != FALSE && text[1] != FALSE) {
        umap <- umap::umap(t(as.matrix(mydata)), ...)
        scores <- data.frame(umap$layout)
        scores$label <- text
        if (controlscale == TRUE) {
            if (scale == 1) {
                p <- ggplot(data = scores, aes(x = X1, y = X2,
                  label = label)) + geom_point(aes(colour = labels),
                  size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                    colour = "black"), axis.text.x = element_text(size = axistextsize,
                    colour = "black"), axis.title.x = element_text(size = axistextsize),
                  axis.title.y = element_text(size = axistextsize),
                  legend.title = element_text(size = legendtextsize),
                  legend.text = element_text(size = legendtextsize)) +
                  labs(colour = legendtitle) + scale_colour_distiller(palette = "Spectral") +
                  geom_text(vjust = "inward", hjust = "inward",
                    size = textlabelsize)+ ggtitle(ggtitle)  + theme(plot.title = element_text(hjust = 0.5)) +
                    theme(plot.title = element_text(size = 15, face = "bold"))
            }
            else if (scale == 2) {
                p <- ggplot(data = scores, aes(x = X1, y = X2,
                  label = label)) + geom_point(aes(colour = labels),
                  size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                    colour = "black"), axis.text.x = element_text(size = axistextsize,
                    colour = "black"), axis.title.x = element_text(size = axistextsize),
                  axis.title.y = element_text(size = axistextsize),
                  legend.title = element_text(size = legendtextsize),
                  legend.text = element_text(size = legendtextsize)) +
                  labs(colour = legendtitle) + scale_colour_gradient(low = low,
                  high = high) + geom_text(vjust = "inward",
                  hjust = "inward", size = textlabelsize)+ gtitle(ggtitle) + theme(plot.title = element_text(hjust = 0.5)) +
                    theme(plot.title = element_text(size = 15, face = "bold"))
            }
            else if (scale == 3) {
                p <- ggplot(data = scores, aes(x = X1, y = X2,
                  label = label)) + geom_point(aes(colour = labels),
                  size = dotsize) + theme_bw() + theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                    colour = "black"), axis.text.x = element_text(size = axistextsize,
                    colour = "black"), axis.title.x = element_text(size = axistextsize),
                  axis.title.y = element_text(size = axistextsize),
                  legend.title = element_text(size = legendtextsize),
                  legend.text = element_text(size = legendtextsize)) +
                  labs(colour = legendtitle) + scale_colour_manual(values = colvec) +
                  geom_text(vjust = "inward", hjust = "inward",
                    size = textlabelsize)+ ggtitle(ggtitle)+ theme(plot.title = element_text(hjust = 0.5)) +
                    theme(plot.title = element_text(size = 15, face = "bold"))
            }
        }
        else {
            p <- ggplot(data = scores, aes(x = X1, y = X2, label = label)) +
                geom_point(aes(colour = labels), size = dotsize ) +
                theme_bw() + theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.text.y = element_text(size = axistextsize,
                  colour = "black"), axis.text.x = element_text(size = axistextsize,
                  colour = "black"), axis.title.x = element_text(size = axistextsize),
                axis.title.y = element_text(size = axistextsize),
                legend.title = element_text(size = legendtextsize),
                legend.text = element_text(size = legendtextsize)) +
                labs(colour = legendtitle) + geom_text(vjust = "inward",
                hjust = "inward", size = textlabelsize)+ ggtitle(ggtitle)
        }
        if (printres == TRUE) {
            message("printing UMAP to current directory...")
            png("UMAPlabeled.png", height = printheight, width = printwidth,
                units = "cm", res = 900, type = "cairo")
            print(p)
            dev.off()
        }
    }
    else if (labels[1] == FALSE && text[1] != FALSE) {
        umap <- umap::umap(t(as.matrix(mydata)), ...)
        scores <- data.frame(umap$layout)
        scores$label <- text
        p <- ggplot(data = scores, aes(x = X1, y = X2, label = label)) +
            geom_point(aes(colour = factor(rep(1, ncol(mydata)))),
                size = dotsize) + theme_bw() + theme(legend.position = "none",
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.y = element_text(size = axistextsize, colour = "black"),
            axis.text.x = element_text(size = axistextsize, colour = "black"),
            axis.title.x = element_text(size = axistextsize),
            axis.title.y = element_text(size = axistextsize)) +
            scale_colour_manual(values = colvec) + geom_text(vjust = "inward",
            hjust = "inward", size = textlabelsize) + ggtitle(ggtitle) + theme(plot.title = element_text(hjust = 0.5)) +
            theme(plot.title = element_text(size = 15, face = "bold"))
        if (printres == TRUE) {
            message("printing UMAP to current directory...")
            png("UMAP.png", height = printheight, width = printwidth,
                units = "cm", res = 900, type = "cairo")
            print(p)
            dev.off()
        }
    }
    p + ggtitle(ggtitle) + theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size = 15, face = "bold"))
    message("done.")
    return(p)
}
