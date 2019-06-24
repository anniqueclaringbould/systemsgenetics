##########################################
#### Check SLE results DEPICT2.0.21/22 ###
##########################################

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(rtracklayer)
library(ggplot2)
library(reshape)
library(tidyr)
library(plotly)
library(gplots)

setwd('/Users/anniqueclaringbould/Documents/Projects/UMCG/Pathway-specific_PGS/DEPICT2/input')

#read in data
sle <- fread('./v2.0.21/SLE/SLE_genePvalues.txt')
reactome <- fread('reactome_predictions.txt')

sle_files <- list.files(path = 'v2.0.22/SLE/', pattern = '_zscore.txt', full.names = TRUE)
sle_z <- lapply(sle_files, fread)

sle2_files <- list.files(path = 'v2.0.22/SLE_same_again/', pattern = '_zscore.txt', full.names = TRUE)
sle2_z <- lapply(sle2_files, fread)

#adjust 
reactome <- reactome[, c('-', 'R-HSA-877300')]
dm <- merge(sle, reactome, by = '-')

dms <- list()

for (i in 1:6){
  #save all enrichment terms
  terms <- sle_z[[i]][,1]
  
  #remove enrichment terms from each dataset
  sle_z[[i]] <- sle_z[[i]][,2:ncol(sle_z[[i]])]
  sle2_z[[i]] <- sle2_z[[i]][,2:ncol(sle2_z[[i]])]
  
  #change column names for each dataset
  colnames(sle_z[[i]]) <- paste0('SLE1', colnames(sle_z[[i]]))
  colnames(sle2_z[[i]]) <- paste0('SLE2', colnames(sle2_z[[i]]))
  
  #bind together all datasets
  dms[[i]] <- do.call("cbind", list(terms,
                                    sle_z[[i]],
                                    sle2_z[[i]]))
}

dm_gc <- as.data.frame(dms[[1]])
dm_gc2 <- dm_gc[, c(2:ncol(dm_gc))]
cors_gc <- cor(dm_gc2, method = 'spearman')

dm_gf <- as.data.frame(dms[[2]])
dm_gf2 <- dm_gf[, c(2:ncol(dm_gf))]
cors_gf <- cor(dm_gf2, method = 'spearman')

dm_gp <- as.data.frame(dms[[3]])
dm_gp2 <- dm_gp[, c(2:ncol(dm_gp))]
cors_gp <- cor(dm_gp2, method = 'spearman')

dm_h <- as.data.frame(dms[[4]])
dm_h2 <- dm_h[, c(2:ncol(dm_h))]
cors_h <- cor(dm_h2, method = 'spearman')

dm_k <- as.data.frame(dms[[5]])
dm_k2 <- dm_k[, c(2:ncol(dm_k))]
cors_k <- cor(dm_k2, method = 'spearman')

dm_r <- as.data.frame(dms[[6]])
dm_r2 <- dm_r[, c(2:ncol(dm_r))]
cors_r <- cor(dm_r2, method = 'spearman')

#plots
hist(-log10(sle$SLE))
hist(sle$SLE)

ggplot(dm, aes(x = SLE, y = `R-HSA-877300`)) +
  geom_point(size = 0.4) +
  xlab('SLE gene p-values') +
  ylab('Interferon gamma signaling Z-score') +
  theme_bw()
ggsave("../plots/sle_reactome_interferon.png")

ggplot(dm, aes(x = -log10(SLE), y = `R-HSA-877300`)) +
  geom_point(size = 0.1) +
  xlab('SLE gene -log10(p-values)') +
  ylab('Interferon gamma signaling Z-score') +
  theme_bw()
ggsave("../plots/sle_reactome_interferon_logp.png")

ggplot(dm_r, aes(x=SLE1SLE, y=SLE2SLE)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-scores SLE, same settings twice \n(DEPICT2.0.22)')
ggsave("../plots/reactome_sle_same_again_v22.png")

ggplot(dm_h, aes(x=SLE1SLE, y=SLE2SLE)) +
  geom_point() +
  theme_bw() +
  ggtitle('HPO Z-scores SLE, same settings twice \n(DEPICT2.0.22)')
ggsave("../plots/hpo_sle_same_again_v22.png")
