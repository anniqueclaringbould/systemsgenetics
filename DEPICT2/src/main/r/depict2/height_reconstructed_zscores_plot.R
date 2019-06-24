################################################
########## Plot height gene Z-scores ###########
################################################

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(rtracklayer)
library(ggplot2)
library(reshape)
library(tidyr)
library(ggrepel)

setwd('/Users/anniqueclaringbould/Documents/Projects/UMCG/Pathway-specific_PGS/DEPICT2/input/')

##### READ IN DATA #####
gene_score <- fread('Height_eigen_Enrichment_normalizedGwasGeneScores_ExHla.txt')
reconst_gene_score <- fread('geneScoresExHla2.txt')
gene_names <- fread('/Users/anniqueclaringbould/Documents/Projects/UMCG/eQTLGen/Annotation_files/ProbeAnnotation_STARv2.3.0e_Ensembl71.txt')
gene_mapping <- fread('ensgV83_skewness025.txt') #from /groups/umcg-wijmenga/scr02/depict2
snp_mapping <- fread('SNPMappings.txt', header = F) #from /groups/umcg-wijmenga/scr02/depict2/HarmoniserAllChrs/SNPMappings.txt 
height_gwas <- fread('Height_snps_pvals.txt') #from /groups/umcg-wijmenga/scr02/depict2/Height
maf <- fread('1000GENOMES-phase_1_EUR_MAF.vcf')

##### ADJUST #####
#merge reconstructed and original height Z-scores
colnames(gene_score) <- c('Gene', 'HeightGeneZScore')
colnames(reconst_gene_score) <- c('Gene', 'ReconstructedGeneZScore')
dm <- merge(reconst_gene_score, gene_score, by = 'Gene')

#add gene names (new dataframe because I do not have names for each ENSG number and don't want to lose anything)
dm2 <- gene_names %>%
  mutate(Gene = Probe) %>%
  select(Gene, Symbol) %>%
  inner_join(dm, by = 'Gene')

#combine SNP mapping information with P-values from height GWAS 
maf$ID <- paste0(maf$CHROM,"_", maf$POS)

snps <- snp_mapping %>%
  mutate(chr = V1,
         pos = V2,
         rs = V3,
         ID = paste0(chr,"_",pos)) %>%
  inner_join(height_gwas, by = 'rs') %>%
  mutate(GWAS_pval = HEIGHT) %>%
  select(chr, pos, rs, GWAS_pval, ID)
rm(snp_mapping)
rm(height_gwas)

maf <- maf[maf$ID %in% snps$ID, ]

snps <- snps %>%
  inner_join(maf, by = 'ID') %>%
  mutate(AF = as.numeric(AF)) %>%
  mutate(MAF = ifelse(AF > 0.5, 1-AF, AF)) %>%
  select(chr, pos, rs, GWAS_pval, MAF)

#make window of 50kb around each tested gene
gene_mapping <- gene_mapping %>%
  mutate(start_window = `Gene Start (bp)` - 50000,
         end_window = `Gene End (bp)` + 50000,
         Gene = `Ensembl Gene ID`,
         chr = as.numeric(`Chromosome Name`)) %>%
  select(Gene, chr, start_window, end_window)

#select the top genes with high reconstructed Z-score
top_genes <- dm %>%
  filter(ReconstructedGeneZScore >3) %>%
  left_join(gene_mapping, by = 'Gene')

#add all SNPs within the gene window of the top genes
top_SNP_per_gene <- setDT(snps)[top_genes,
                   on = .(chr, pos > start_window, pos < end_window)]

#select the most significant SNP in each window
top_SNP_per_gene <- top_SNP_per_gene %>%
  filter(!is.na(Gene)) %>%
  group_by(Gene) %>%
  filter(GWAS_pval == min(GWAS_pval)) %>%
  select(Gene, rs, ReconstructedGeneZScore, HeightGeneZScore, MAF)

cor <- cor(top_SNP_per_gene$MAF, top_SNP_per_gene$HeightGeneZScore, method = 'spearman')
cor.test(top_SNP_per_gene$MAF, top_SNP_per_gene$HeightGeneZScore, method = 'spearman')

##### PLOT #####
#plot reconstructed vs. original DEPICT Z-score
ggplot(dm, aes(y = ReconstructedGeneZScore, x = HeightGeneZScore)) +
  geom_point(size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=3, linetype="dashed", color = "red") +
  geom_text_repel(aes(label=ifelse
                      (ReconstructedGeneZScore >3 & HeightGeneZScore < 1, 
                        as.character(Gene), 
                        '')),
                  size = 3)
ggsave('../plots/height_reconstructed_zsore.png', height = 10, width = 10)

#same incl. gene names
ggplot(dm2, aes(y = ReconstructedGeneZScore, x = HeightGeneZScore)) +
  geom_point(size = 0.3) +
  theme_bw() +
  geom_hline(yintercept=3, linetype="dashed", color = "red") +
  geom_text_repel(aes(label=ifelse
                      (ReconstructedGeneZScore >3 & HeightGeneZScore < 1, 
                        as.character(Symbol), 
                        '')),
                  size = 3)
ggsave('../plots/height_reconstructed_zsore_gene_names.png', height = 10, width = 10)

#For each gene with ReconstructedGeneZScore > 3, take 50 KB around gene
#Take top SNP for height summary stats
#Plot MAF for this SNP with HeightGeneZScore:
ggplot(top_SNP_per_gene, aes(x = MAF, y = HeightGeneZScore)) +
  geom_point() +
  theme_bw() +
  ylab('Z-score for height, as calculated using DEPICT2') +
  xlab('Allele frequency for the most significant SNP \nin the 50kb window around this gene') +
  geom_text_repel(aes(label=as.character(Gene)), size = 3) +
  annotate("text", x = 0.15, y = 2, label = paste0("Correlation =\n", round(cor,2)))
ggsave('../plots/height_zsore_MAF_close_SNP.png', height = 8, width = 8)
 