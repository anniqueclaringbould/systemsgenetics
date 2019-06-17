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
sim_r <- fread('simulated/simulated_Reactome_Enrichment_zscore.txt')
sim_gc <- fread('simulated/simulated_GO_C_Enrichment_zscore.txt')
sim_gf <- fread('simulated/simulated_GO_F_Enrichment_zscore.txt')
sim_gp <- fread('simulated/simulated_GO_P_Enrichment_zscore.txt')
sim_h <- fread('simulated/simulated_HPO_Enrichment_zscore.txt')
sim_k <- fread('simulated/simulated_KEGG_Enrichment_zscore.txt')

sim2_r <- fread('simulated2/simulated_Reactome_Enrichment_zscore.txt')
sim2_gc <- fread('simulated2/simulated_GO_C_Enrichment_zscore.txt')
sim2_gf <- fread('simulated2/simulated_GO_F_Enrichment_zscore.txt')
sim2_gp <- fread('simulated2/simulated_GO_P_Enrichment_zscore.txt')
sim2_h <- fread('simulated2/simulated_HPO_Enrichment_zscore.txt')
sim2_k <- fread('simulated2/simulated_KEGG_Enrichment_zscore.txt')

height_r <- fread('Height_Reactome_Enrichment_zscore.txt')
height_gc <- fread('Height_GO_C_Enrichment_zscore.txt')
height_gf <- fread('Height_GO_F_Enrichment_zscore.txt')
height_gp <- fread('Height_GO_P_Enrichment_zscore.txt')
height_h <- fread('Height_HPO_Enrichment_zscore.txt')
height_k <- fread('Height_KEGG_Enrichment_zscore.txt')

height2_r <- fread('Height_no_cl/Height_Reactome_Enrichment_zscore.txt')
height2_gc <- fread('Height_no_cl/Height_GO_C_Enrichment_zscore.txt')
height2_gf <- fread('Height_no_cl/Height_GO_F_Enrichment_zscore.txt')
height2_gp <- fread('Height_no_cl/Height_GO_P_Enrichment_zscore.txt')
height2_h <- fread('Height_no_cl/Height_HPO_Enrichment_zscore.txt')
height2_k <- fread('Height_no_cl/Height_KEGG_Enrichment_zscore.txt')

sle_r <- fread('SLE_Reactome_Enrichment_zscore.txt')
sle_gc <- fread('SLE_GO_C_Enrichment_zscore.txt')
sle_gf <- fread('SLE_GO_F_Enrichment_zscore.txt')
sle_gp <- fread('SLE_GO_P_Enrichment_zscore.txt')
sle_h <- fread('SLE_HPO_Enrichment_zscore.txt')
sle_k <- fread('SLE_KEGG_Enrichment_zscore.txt')

ibd_r <- fread('IBD_Reactome_Enrichment_zscore.txt')
ibd_gc <- fread('IBD_GO_C_Enrichment_zscore.txt')
ibd_gf <- fread('IBD_GO_F_Enrichment_zscore.txt')
ibd_gp <- fread('IBD_GO_P_Enrichment_zscore.txt')
ibd_h <- fread('IBD_HPO_Enrichment_zscore.txt')
ibd_k <- fread('IBD_KEGG_Enrichment_zscore.txt')

##### ADJUST
colnames(sim_r) <- paste0('sim1_', colnames(sim_r))
colnames(sim2_r) <- paste0('sim2_', colnames(sim2_r))
colnames(height2_r) <- paste0('no_cl_', colnames(height2_r))

dm_r <- sle_r %>%
  left_join(height_r, by = '-') %>%
  left_join(height2_r, by = c('-' = 'no_cl_-')) %>%
  left_join(ibd_r, by = '-') %>%
  #left_join(sim_r, by = c('-' = 'sim_-')) %>%
  left_join(sim2_r, by = c('-' = 'sim2_-'))

colnames(sim_gc) <- paste0('sim1_', colnames(sim_gc))
colnames(sim2_gc) <- paste0('sim2_', colnames(sim2_gc))
colnames(height2_gc) <- paste0('no_cl_', colnames(height2_gc))

dm_gc <- sle_gc %>%
  left_join(height_gc, by = '-') %>%
  left_join(height2_gc, by = c('-' = 'no_cl_-')) %>%
  left_join(ibd_gc, by = '-') %>%
  #left_join(sim_gc, by = c('-' = 'sim_-')) %>%
  left_join(sim2_gc, by = c('-' = 'sim2_-'))

colnames(sim_gf) <- paste0('sim1_', colnames(sim_gf))
colnames(sim2_gf) <- paste0('sim2_', colnames(sim2_gf))
colnames(height2_gf) <- paste0('no_cl_', colnames(height2_gf))

dm_gf <- sle_gf %>%
  left_join(height_gf, by = '-') %>%
  left_join(height2_gf, by = c('-' = 'no_cl_-')) %>%
  left_join(ibd_gf, by = '-') %>%
  #left_join(sim_gf, by = c('-' = 'sim_-')) %>%
  left_join(sim2_gf, by = c('-' = 'sim2_-'))

colnames(sim_gp) <- paste0('sim1_', colnames(sim_gp))
colnames(sim2_gp) <- paste0('sim2_', colnames(sim2_gp))
colnames(height2_gp) <- paste0('no_cl_', colnames(height2_gp))

dm_gp <- sle_gp %>%
  left_join(height_gp, by = '-') %>%
  left_join(height2_gp, by = c('-' = 'no_cl_-')) %>%
  left_join(ibd_gp, by = '-') %>%
  #left_join(sim_gp, by = c('-' = 'sim_-')) %>%
  left_join(sim2_gp, by = c('-' = 'sim2_-'))

colnames(sim_h) <- paste0('sim1_', colnames(sim_h))
colnames(sim2_h) <- paste0('sim2_', colnames(sim2_h))
colnames(height2_h) <- paste0('no_cl_', colnames(height2_h))

dm_h <- sle_h %>%
  left_join(height_h, by = '-') %>%
  left_join(height2_h, by = c('-' = 'no_cl_-')) %>%
  left_join(ibd_h, by = '-') %>%
  #left_join(sim_h, by = c('-' = 'sim_-')) %>%
  left_join(sim2_h, by = c('-' = 'sim2_-'))

colnames(sim_k) <- paste0('sim1_', colnames(sim_k))
colnames(sim2_k) <- paste0('sim2_', colnames(sim2_k))
colnames(height2_k) <- paste0('no_cl_', colnames(height2_k))

dm_k <- sle_k %>%
  left_join(height_k, by = '-') %>%
  left_join(height2_k, by = c('-' = 'no_cl_-')) %>%
  left_join(ibd_k, by = '-') %>%
  #left_join(sim_k, by = c('-' = 'sim_-')) %>%
  left_join(sim2_k, by = c('-' = 'sim2_-'))

##### CORRELATIONS #####
dm_r2 <- dm_r[, c(2:ncol(dm_r))]
cors_r <- cor(dm_r2, method = 'spearman')

dm_gc2 <- dm_gc[, c(2:ncol(dm_gc))]
cors_gc <- cor(dm_gc2, method = 'spearman')

dm_gf2 <- dm_gf[, c(2:ncol(dm_gf))]
cors_gf <- cor(dm_gf2, method = 'spearman')

dm_gp2 <- dm_gp[, c(2:ncol(dm_gp))]
cors_gp <- cor(dm_gp2, method = 'spearman')

dm_h2 <- dm_h[, c(2:ncol(dm_h))]
cors_h <- cor(dm_h2, method = 'spearman')

dm_k2 <- dm_k[, c(2:ncol(dm_k))]
cors_k <- cor(dm_k2, method = 'spearman')


##### PLOT #####
Colors=c("blue","white","red")
Colors=colorRampPalette(Colors)(100)

pdf('../../plots/heatmaps_v2.0.21.pdf')
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
  geom_point(size = 0.5) +
  theme_bw() +
  ggtitle('HPO Z-score height with and without -cl \n(DEPICT2.0.21)')
ggplotly(p)
ggsave("../../plots/reactome_height_height_no_cl_v21.png")
