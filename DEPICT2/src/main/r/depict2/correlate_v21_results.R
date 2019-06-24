#################################################
#### Correlate GWAS Z-scores for DEPICT2.0.21 ###
#################################################

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

setwd('/Users/anniqueclaringbould/Documents/Projects/UMCG/Pathway-specific_PGS/DEPICT2/input/v2.0.21')

##### READ IN DATA #####
sim_files <- list.files(path = './simulated', pattern = '_zscore.txt', full.names = TRUE)
sim <- lapply(sim_files, fread)
sim2_files <- list.files(path = './simulated2', pattern = '_zscore.txt', full.names = TRUE)
sim2 <- lapply(sim2_files, fread)
sim_no_cl_files <- list.files(path = './simulated_no_cl', pattern = '_zscore.txt', full.names = TRUE)
sim_no_cl <- lapply(sim_no_cl_files, fread)
sim2_no_cl_files <- list.files(path = './simulated2_no_cl', pattern = '_zscore.txt', full.names = TRUE)
sim2_no_cl <- lapply(sim2_no_cl_files, fread)
height_files <- list.files(path = './Height', pattern = '_zscore.txt', full.names = TRUE)
height <- lapply(height_files, fread)
height_no_cl_files <- list.files(path = './Height_no_cl', pattern = '_zscore.txt', full.names = TRUE)
height_no_cl <- lapply(height_no_cl_files, fread)
#ibd_files <- list.files(path = './IBD', pattern = '_zscore.txt', full.names = TRUE)
#ibd <- lapply(ibd_files, fread)
#ibd_no_cl_files <- list.files(path = './IBD_no_cl', pattern = '_zscore.txt', full.names = TRUE)
#ibd_no_cl <- lapply(ibd_no_cl_files, fread)
sle_files <- list.files(path = './SLE', pattern = '_zscore.txt', full.names = TRUE)
sle <- lapply(sle_files, fread)
sle_no_cl_files <- list.files(path = './SLE_no_cl', pattern = '_zscore.txt', full.names = TRUE)
sle_no_cl <- lapply(sle_no_cl_files, fread)

dms <- list()

for (i in 1:6){
  #save all enrichment terms
  terms <- sim[[i]][,1]
  
  #remove enrichment terms from each dataset
  sim[[i]] <- sim[[i]][,2:ncol(sim[[i]])]
  sim2[[i]] <- sim2[[i]][,2:ncol(sim2[[i]])]
  sim_no_cl[[i]] <- sim_no_cl[[i]][,2:ncol(sim_no_cl[[i]])]
  #sim2_no_cl[[i]] <- sim2_no_cl[[i]][,2:ncol(sim2_no_cl[[i]])]
  height[[i]] <- height[[i]][,2:ncol(height[[i]])]
  height_no_cl[[i]] <- height_no_cl[[i]][,2:ncol(height_no_cl[[i]])]
  sle[[i]] <- sle[[i]][,2:ncol(sle[[i]])]
  sle_no_cl[[i]] <- sle_no_cl[[i]][,2:ncol(sle_no_cl[[i]])]
  #ibd[[i]] <- ibd[[i]][,2:ncol(ibd[[i]])]
  #ibd_no_cl[[i]] <- ibd_no_cl[[i]][,2:ncol(ibd_no_cl[[i]])]
  
  #change column names for each dataset
  colnames(sim[[i]]) <- paste0('sim1_', colnames(sim[[i]]))
  colnames(sim2[[i]]) <- paste0('sim2_', colnames(sim2[[i]]))
  colnames(sim_no_cl[[i]]) <- paste0('sim1_no_cl_', colnames(sim_no_cl[[i]]))
  #colnames(sim2_no_cl[[i]]) <- paste0('sim2_no_cl_', colnames(sim2_no_cl[[i]]))
  colnames(height_no_cl[[i]]) <- paste0('no_cl_', colnames(height_no_cl[[i]]))
  colnames(sle_no_cl[[i]]) <- paste0('no_cl_', colnames(sle_no_cl[[i]]))
  #colnames(ibd_no_cl[[i]]) <- paste0('no_cl_', colnames(ibd_no_cl[[i]]))
  
  #bind together all datasets
  dms[[i]] <- do.call("cbind", list(terms,
                        sim[[i]],
                        sim2[[i]],
                        sim_no_cl[[i]],
                        #sim2_no_cl[[i]],
                        height[[i]],
                        height_no_cl[[i]],
                        sle[[i]],
                        sle_no_cl[[i]]))
                        #ibd[[i]],
                        #ibd_no_cl[[i]]))

  }

##### CORRELATIONS #####
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

##### PLOT #####
Colors=c("blue","white","red")
Colors=colorRampPalette(Colors)(100)

pdf('../../plots/heatmaps_v2.0.21.pdf', height = 12, width = 12)
heatmap.2(cors_r, cellnote = round(cors_r,2), main = 'Reactome', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'both', margins =c(11,11),Colv = NA, Rowv = NA,
          col = Colors)
heatmap.2(cors_gc, cellnote = round(cors_gc,2), main = 'GO C', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(11,11),Colv = NA, Rowv = NA,
          col = Colors)
heatmap.2(cors_gf, cellnote = round(cors_gf,2), main = 'GO F', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(11,11),Colv = NA, Rowv = NA,
          col = Colors)
heatmap.2(cors_gp, cellnote = round(cors_gp,2), main = 'GO P', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(11,11),Colv = NA, Rowv = NA,
          col = Colors)
heatmap.2(cors_h, cellnote = round(cors_h,2), main = 'HPO', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(11,11),Colv = NA, Rowv = NA,
          col = Colors)
heatmap.2(cors_k, cellnote = round(cors_k,2), main = 'KEGG', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(11,11),Colv = NA, Rowv = NA,
          col = Colors)
dev.off()

ggplot(dm_r, aes(x=sim2_simulated_z_1, y=sim2_simulated_z_2)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-scores simulated GWAS \n(DEPICT2.0.21)')
ggsave("../../plots/reactome_sim_gwas_z1_z2_v21.png")

ggplot(dm_r, aes(x=sim2_simulated_z_1, y=sim2_expected_z)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-score simulated GWAS expected vs.\nsimulated Z1 (DEPICT2.0.21)')
ggsave("../../plots/reactome_sim_gwas_z1_expz_v21.png")

ggplot(dm_r, aes(x=sim2_expected_z, y=HEIGHT)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-score simulated GWAS vs height \n(DEPICT2.0.21)')
ggsave("../../plots/reactome_sim_gwas_expz_height_v21.png")

ggplot(dm_r, aes(x=SLE, y=HEIGHT)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-score SLE vs height \n(DEPICT2.0.21)')
ggsave("../../plots/reactome_sle_height_v21.png")

p <- ggplot(dm_r, aes(x=IBD, y=HEIGHT, text = `-`)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-score IBD vs height \n(DEPICT2.0.21)')
ggsave("../../plots/reactome_ibd_height_v21.png")

ggplot(dm_r, aes(x=no_cl_HEIGHT, y=HEIGHT)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-score height with and without -cl \n(DEPICT2.0.21)')
ggsave("../../plots/reactome_height_height_no_cl_v21.png")

p <- ggplot(dm_h, aes(x=no_cl_HEIGHT, y=HEIGHT, text = `-`)) +
  geom_point() +
  theme_bw() +
  ggtitle('HPO Z-score height with and without -cl \n(DEPICT2.0.21)')
ggplotly(p)
ggsave("../../plots/hpo_height_height_no_cl_v21.png")
