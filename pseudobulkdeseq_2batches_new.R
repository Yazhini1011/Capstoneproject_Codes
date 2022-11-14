# script to perform pseudo-bulk DGA
setwd("C:\\Users\\yazhinir\\OneDrive - A STAR\\Desktop\\new_integration")

library(ExperimentHub)
library(Seurat)
library(DESeq2)
library(tidyverse)

RPCA_LISI@meta.data
RPCA_LISI@meta.data <- separate(RPCA_LISI@meta.data, col = 'Sample', into = c('samplenumber','ethnicity'), sep = '_')
RPCA_LISI@meta.data <- separate(RPCA_LISI@meta.data, col = 'ethnicity', into = c('race','healthy'), sep = '~')
healthy_only <- subset(RPCA_LISI, subset = healthy=="H")
saveRDS(healthy_only,'healthy_only.rds')
healthy_only<-readRDS("healthy_only.rds")
DimPlot(healthy_only,reduction = "umap", group.by = "MCT")
DefaultAssay(healthy_only)<-'integrated'
healthy
# get data

chinese_vs_malay <- subset (healthy_only, subset=race==c('chinese','malay' ))
malay_vs_indian<-subset (healthy_only, subset=race==c('malay','indian' ))
indian_vs_chinese<-subset(healthy_only,subset=race==c('indian','chinese'))

# QC and filtering
# explore QC


# get mito percent
seu.obj$mitoPercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')
View(seu.obj@meta.data)

# filter
seu.filtered <- subset(seu.obj, subset = nFeature_originalexp > 200 & nFeature_originalexp < 2500 &
                         nCount_originalexp > 800 & 
                         mitoPercent < 5 &
                         multiplets == 'singlet')
chinese_vs_malay <- subset (healthy_only, subset=race==c('chinese','malay' ))
malay_vs_indian<-subset (healthy_only, subset=race==c('malay','indian' ))

seu.obj
seu.filtered
DefaultAssay(chinese_vs_malay)<-'integrated'
# run Seurat's standard workflow steps
chinese_vs_malay<- NormalizeData(chinese_vs_malay)
chinese_vs_malay <- FindVariableFeatures(chinese_vs_malay)
chinese_vs_malay<- ScaleData(chinese_vs_malay)
chinese_vs_malay <- RunPCA(chinese_vs_malay)
ElbowPlot(chinese_vs_malay)
chinese_vs_malay <- RunUMAP(chinese_vs_malay, dims = 1:20)

# visualize 
healthy_onlyintegrated_final_plot <- DimPlot(chinese_vs_malay, reduction = 'umap', group.by = 'MCT', label = TRUE)
ethnicity <- DimPlot(chinese_vs_malay, reduction = 'umap', group.by = 'race')

healthy_onlyintegrated_final_plot|ethnicity

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(healthy_only@meta.data)
chinese_vs_malay$Sample<-NULL
chinese_vs_malay$race <- paste((chinese_vs_malay$race), "-", sep = "")

chinese_vs_malay$Sample <- paste0(chinese_vs_malay$race, chinese_vs_malay$samplenumber,chinese_vs_malay$batch)
malay_vs_indian
DefaultAssay(chinese_vs_malay)

cts <- AggregateExpression(chinese_vs_malay, 
                           group.by = c("MCT", "Sample"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')



# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_Monocytes <- cts.split.modified$Monocytes


# 2. generate sample level metadata
colData <- data.frame(Sample = colnames(counts_Monocytes))
colData$Sample
colData$condition <- coldata.df$ethni
coldata.df<-as.data.frame(colData)
coldata.df<- separate(coldata.df, col='Sample', into = c('ethni' ,'idb'), sep = '-')

rownames(colData)<-colData$Sample
colData


# get more information from metadata




# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_Tcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=5
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)



# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res05 <- results(dds, name = "condition_malay_vs_chinese",alpha = 0.05)
csv_malay_vs_chinese_Monocytes<-write.csv(res05)

summary(res05)

resLFC <- lfcShrink(dds, coef="condition_treated_vs_untreated", type="apeglm")


####################################################################malayvsindian################################

DefaultAssay(malay_vs_indian)<-'integrated'
# run Seurat's standard workflow steps
malay_vs_indian<- NormalizeData(malay_vs_indian)
malay_vs_indian <- FindVariableFeatures(malay_vs_indian)
malay_vs_indian<- ScaleData(malay_vs_indian)
malay_vs_indian <- RunPCA(malay_vs_indian)
ElbowPlot(malay_vs_indian)
malay_vs_indian <- RunUMAP(malay_vs_indian, dims = 1:20)

# visualize 
healthy_onlyintegrated_final_plot <- DimPlot(chinese_vs_malay, reduction = 'umap', group.by = 'MCT', label = TRUE)
ethnicity <- DimPlot(chinese_vs_malay, reduction = 'umap', group.by = 'race')

healthy_onlyintegrated_final_plot|ethnicity

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(healthy_only@meta.data)
malay_vs_indian$Sample<-NULL
malay_vs_indian$race <- paste((malay_vs_indian$race), "-", sep = "")

malay_vs_indian$Sample <- paste0(malay_vs_indian$race, malay_vs_indian$samplenumber,malay_vs_indian$batch)
malay_vs_indian
DefaultAssay(malay_vs_indian)

cts <- AggregateExpression(malay_vs_indian, 
                           group.by = c("MCT", "Sample"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')
make.names(unique = F,names = row.names(cts.split))


# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_T <- cts.split.modified$T


# 2. generate sample level metadata
colData <- data.frame(Sample = colnames(counts_Tcell))
colData$Sample
coldata.df<-as.data.frame(colData)
coldata.df<- separate(coldata.df, col='Sample', into = c('ethni' ,'idb'), sep = '-')
colData$condition <- coldata.df$ethni
coldata.df<-as.data.frame(colData)
rownames(colData)<-colData$Sample



# get more information from metadata




# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_Tcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=5
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds)
dds <- DESeq(dds,test="LRT",reduced = ~ 1,useT=TRUE,minmu=1e-6,  minReplicatesForReplace=Inf)



# Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res05 <- results(dds, name = "condition_malay_vs_chinese",alpha = 0.05)
csv_malay_vs_indian_B<-write.csv(res05)
plotMA(res05, ylim=c(-5,5))
summary(res05)

colData

rld <- rlog(dds, blind=TRUE)

# Plot PCA

DESeq2::plotPCA(rld, intgroup = "condition")
plotDispEsts(dds)

# Output results of Wald test for contrast for stim vs ctrl
levels(colData$condition)[2]
levels(cluster_metadata$group_id)[1]

contrast <- c("condition", levels(colData$condition)[2], levels(colData$condition)[1])

# resultsNames(dds)
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 contrast =  contrast,
                 res=res)
contrast <- c("condition", colData, levels(cluster_metadata$group_id)[1])
#########################################################indianvschinese#############################################
chinese_vs_indian <- subset (healthy_only, subset=race==c('chinese','indian' ))


DefaultAssay(chinese_vs_indian)<-'integrated'

# run Seurat's standard workflow steps
chinese_vs_indian<- NormalizeData(chinese_vs_indian)
chinese_vs_indian <- FindVariableFeatures(chinese_vs_indian)
chinese_vs_indian<- ScaleData(chinese_vs_indian)
chinese_vs_indian <- RunPCA(chinese_vs_indian)
ElbowPlot(chinese_vs_indian)
chinese_vs_indian <- RunUMAP(chinese_vs_indian, dims = 1:20)
(chinese_vs_indian$race)
# visualize 
healthy_onlyintegrated_final_plot <- DimPlot(chinese_vs_indian, reduction = 'umap', group.by = 'MCT', label = TRUE)
ethnicity <- DimPlot(chinese_vs_indian, reduction = 'umap', group.by = 'race')

healthy_onlyintegrated_final_plot|ethnicity

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(healthy_only@meta.data)
chinese_vs_indian$Sample<-NULL
chinese_vs_indian$race <- paste((chinese_vs_indian$race), "-", sep = "")

chinese_vs_indian$Sample <- paste0(chinese_vs_indian$race, chinese_vs_indian$samplenumber,chinese_vs_indian$batch)
malay_vs_indian
DefaultAssay(chinese_vs_indian)

cts <- AggregateExpression(chinese_vs_indian, 
                           group.by = c("MCT", "Sample"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)

cts <- cts$RNA

# transpose
cts.t <- t(cts)


# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame
cts.split <- split.data.frame(cts.t,
                              f = factor(splitRows))

# fix colnames and transpose

cts.split.modified <- lapply(cts.split, function(x){
  rownames(x) <- gsub('.*_(.*)', '\\1', rownames(x))
  t(x)
  
})

#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')



# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_Bcell <- cts.split.modified$B
counts_Monocytes<-cts.split.modified$Monocytes
counts_Bcell <- cts.split.modified$B
counts_NKcell <- cts.split.modified$NK
# 2. generate sample level metadata
colData <- data.frame(Sample = colnames(counts_Bcell))
colData$Sample
coldata.df<-as.data.frame(colData)
coldata.df<- separate(coldata.df, col='Sample', into = c('ethni' ,'idb'), sep = '--')
colData$condition <- coldata.df$ethni
coldata.df<-as.data.frame(colData)
rownames(colData)<-colData$Sample
colData



# get more information from metadata




# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_Bcell,
                              colData = colData,
                              design = ~ condition)

# filter
keep <- rowSums(counts(dds)) >=5
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds,test="LRT",reduced = ~ 1,useT=TRUE,minmu=1e-6,  minReplicatesForReplace=Inf)


rld <- rlog(dds, blind=TRUE)

plotDispEsts(dds)
plotMA(dds)

# Plot PCA
DESeq2::plotPCA()
DESeq2::plotPCA(rld, intgroup = "condition",returnData=F)
p+geom_label(aes(label = colData$Sample))


dds_replace<-replaceOutliers(
  dds,
  trim = 0.2,
  'chinese--219'
  )
  
                                                                                                                                                                        # Check the coefficients for the comparison
resultsNames(dds)

# Generate results object
res05 <- results(dds, name = "condition_indian_vs_chinese",alpha = 0.05)
csv_malay_vs_indian_B<-write.csv(res05)

summary(res05)

#############################################################forplotting#############################################
levels(factor(colData$condition))[2]
levels(colData$condition)[1]
contrast <- c("condition", levels(factor(colData$condition))[2], levels(factor(colData$condition))[1])
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 contrast =  contrast,
                 res=res)
summary(res)
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble()

padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

## ggplot of top genes
normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=5)


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)


# Check significant genes output
sig_res


colData$condition


padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)

# Check significant genes output
sig_res

gathered_top20_sig <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
gathered_top20_sig <- inner_join(colData[, c("Sample", "condition" )], gathered_top20_sig, by = c("Sample" = "samplename"))

## plot using ggplot2
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = condition), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top  Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))



rld <- rlog(dds, blind=TRUE)

plotDispEsts(dds)

# Plot PCA
DESeq2::plotPCA()
DESeq2::plotPCA(rld, intgroup = "Sample",returnData=F)
write.csv(sig_res,
          paste0("results", clusters[1], "_", levels(cluster_metadata$sample)[2], "_vs_", levels(cluster_metadata$sample)[1], "_sig_genes.csv"),
          quote = FALSE, 
          row.names = FALSE)


res_table_thres <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.5)

## Volcano plot
ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot of B MEMORY cells in INDIAN relative to Malay") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous_(limits = c(0,20)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

# Extract normalized counts for only the significant genes
sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% sig_res$gene)

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

