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

setwd('/Users/anniqueclaringbould/Documents/Projects/UMCG/Pathway-specific_PGS/DEPICT2/input/v2.0.22')

##### READ IN DATA #####
height_files <- list.files(path = './Height_v22', pattern = '_zscore.txt', full.names = TRUE)
height <- lapply(height_files, fread)
height_fngp_files <- list.files(path = './Height_v22_fngp', pattern = '_zscore.txt', full.names = TRUE)
height_fngp <- lapply(height_fngp_files, fread)
height_fngp_2mb_files <- list.files(path = './Height_v22_fngp_2Mb', pattern = '_zscore.txt', full.names = TRUE)
height_fngp_2mb <- lapply(height_fngp_2mb_files, fread)
height_fngp_fnpp_files <- list.files(path = './Height_v22_fngp_fnpp', pattern = '_zscore.txt', full.names = TRUE)
height_fngp_fnpp <- lapply(height_fngp_fnpp_files, fread)
height_fnpp_files <- list.files(path = './Height_v22_fnpp', pattern = '_zscore.txt', full.names = TRUE)
height_fnpp <- lapply(height_fnpp_files, fread)

dms <- list()

for (i in 1:7){
  #save all enrichment terms
  terms <- height[[i]][,1]
  
  #remove enrichment terms from each dataset
  height[[i]] <- height[[i]][,2:ncol(height[[i]])]
  height_fngp[[i]] <- height_fngp[[i]][,2:ncol(height_fngp[[i]])]
  height_fngp_2mb[[i]] <- height_fngp_2mb[[i]][,2:ncol(height_fngp_2mb[[i]])]
  height_fngp_fnpp[[i]] <- height_fngp_fnpp[[i]][,2:ncol(height_fngp_fnpp[[i]])]
  height_fnpp[[i]] <- height_fnpp[[i]][,2:ncol(height_fnpp[[i]])]
  
  #change column names for each dataset
  colnames(height[[i]]) <- paste0('height_', colnames(height[[i]]))
  colnames(height_fngp[[i]]) <- paste0('height_fngp_', colnames(height_fngp[[i]]))
  colnames(height_fngp_2mb[[i]]) <- paste0('height_fngp_2mb_', colnames(height_fngp_2mb[[i]]))
  colnames(height_fngp_fnpp[[i]]) <- paste0('height_fngp_fnpp_', colnames(height_fngp_fnpp[[i]]))
  colnames(height_fnpp[[i]]) <- paste0('height_fnpp_', colnames(height_fnpp[[i]]))
  
  #bind together all datasets
  dms[[i]] <- do.call("cbind", list(terms,
                        height[[i]],
                        height_fngp[[i]],
                        height_fngp_2mb[[i]],
                        height_fngp_fnpp[[i]],
                        height_fnpp[[i]]))
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

dm_g <- as.data.frame(dms[[7]])
dm_g2 <- dm_r[, c(2:ncol(dm_g))]
cors_g <- cor(dm_g2, method = 'spearman')

##### PLOT #####

Colors=c("blue","white","red")
Colors=colorRampPalette(Colors)(200)
Breaks = c(seq(-1,1,by=0.01))

pdf('../../plots/heatmaps_height_v2.0.22.pdf', height = 8, width = 8)
heatmap.2(cors_r, cellnote = round(cors_r,2), main = 'Reactome', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'both', margins =c(16,16),Colv = NA, Rowv = NA,
          col = Colors, breaks = Breaks)
heatmap.2(cors_gc, cellnote = round(cors_gc,2), main = 'GO C', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(16,16),Colv = NA, Rowv = NA,
          col = Colors, breaks = Breaks)
heatmap.2(cors_gf, cellnote = round(cors_gf,2), main = 'GO F', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(16,16),Colv = NA, Rowv = NA,
          col = Colors, breaks = Breaks)
heatmap.2(cors_gp, cellnote = round(cors_gp,2), main = 'GO P', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(16,16),Colv = NA, Rowv = NA,
          col = Colors, breaks = Breaks)
heatmap.2(cors_h, cellnote = round(cors_h,2), main = 'HPO', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(16,16),Colv = NA, Rowv = NA,
          col = Colors, breaks = Breaks)
heatmap.2(cors_k, cellnote = round(cors_k,2), main = 'KEGG', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(16,16),Colv = NA, Rowv = NA,
          col = Colors, breaks = Breaks)
heatmap.2(cors_g, cellnote = round(cors_g,2), main = 'GTEx', 
          notecol="black", density.info='none', trace = 'none',
          dendogram = 'none', margins =c(16,16),Colv = NA, Rowv = NA,
          col = Colors, breaks = Breaks)
dev.off()

ggplot(dm_r, aes(x=height_HEIGHT, y=height_fngp_HEIGHT)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-scores height \n(DEPICT2.0.22)')
ggsave("../../plots/reactome_height_height_fngp.png")

ggplot(dm_r, aes(x=height_HEIGHT, y=height_fngp_fnpp_HEIGHT)) +
  geom_point() +
  theme_bw() +
  ggtitle('Reactome Z-scores height \n(DEPICT2.0.22)')
ggsave("../../plots/reactome_height_height_fngp_fnpp.png")
