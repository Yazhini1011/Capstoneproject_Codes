setwd("C:\\Users\\yazhinir\\Downloads")
library(Seurat)
library(SeuratData)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(dittoSeq)
library(DESeq2)
BiocManager::install("dittoSeq")
install.packages("rlang")
remotes::install_github("jmcphers/rsrecovr")
rsrecovr::recovr()
Bcells_integrated<-readRDS("Bcellintegrated_withmeta.rds")
DefaultAssay(Bcells_integrated)<-'integrated'
# Run the standard workflow for visualization and clustering
Bcells_integrated  <- ScaleData(Bcells_integrated, verbose = FALSE)
Bcells_integrated  <- RunPCA(Bcells_integrated , npcs = 30, verbose = FALSE)
ElbowPlot(Bcells_intergrated, ndims =30)
Bcells_intergrated  <- RunUMAP(Bcells_intergrated , reduction = "pca", dims = 1:17)
Bcells_intergrated  <- FindNeighbors(Bcells_intergrated, reduction = "pca", dims = 1:17)
Bcells_intergrated  <- FindClusters(Bcells_intergrated, resolution = 1)
 #to plot by different metadatcolums 
p1<-DimPlot(Bcells_intergrated, reduction = "umap", group.by = "seurat_clusters")
p2<-DimPlot(Bcells_intergrated, reduction = "umap", group.by = "lib")
p3<-DimPlot(Bcells_intergrated, reduction = "umap", group.by = "ethnicity")
p1|p2|p3
DefaultAssay(Bcells_intergrated)<-'RNA'
##################featureplot###############################
FeaturePlot(Bcells_intergrated,features=c("CD19","MS4A1","IGHD","IL4R","TCL1A","PLPP5","CXCR4"))#NAIVE
FeaturePlot(Bcells_intergrated,features=c("MX1","SAMD9L”,”IFI44L”,”XAF1”,”IFITM1","MS4A1"))#NAIVE B1
FeaturePlot(Bcells_intergrated,features=c("CD27","TNFRSF13B","CRIP1","MS4A1"))#MEMORY
FeaturePlot(Bcells_intergrated,features=c("CD27","CD19","ITGAX","FCRL5","FCRL3"))#ATYPICAL
FeaturePlot(Bcells_intergrated,features=c("CD19","CD38","CD43","CD27","CD1C","LINC01857"))
FeaturePlot(Bcells_intergrated,features=c("CD19","CD38","CD43","CD27","CD1C","LINC01857"))
FeaturePlot(Bcells_integrated,features=c("CD19","IGHD","IL4R","TCL1A","PLPP5","CXCR4","MX1","SAMD9L","IFI44L","XAF1”,”IFITM1","MS4A1","CD27","TNFRSF13B","CRIP1","ITGAX","FCRL5","FCRL3","IGHM","IGHD"))

Bcells_intergrated <- RenameIdents(Bcells_intergrated, `0` = "Naive B", `1` = "Naive B", `2` = "Memory B",
                                `3` = "Naive B", `4` = "Memory B", `5` = "Naive B", `6` = "Memory B", `7` = "Memory B", `8` = "Memory B", `9` = "Naive B",
                                `10` = "Atypical B", `11` = " AV Naive B", `12` = "Memory B", `13` = "Memory B", `14` = "Naive B",`15` = "Memory B",`16` = "Naive B")
DimPlot(Bcells_intergrated, label = T)
Bcells_intergrated$bsubtypes<-Idents(Bcells_intergrated)
table(Bcells_intergrated$bsubtypes)


Idents(Bcells_intergrated) <- factor(Idents(Bcells_intergrated), levels = c("AV Naive B","Memory B","Naive B","Atypical B"))
markers.to.plot <- c("CD19","IGHD","IL4R","TCL1A","PLPP5","CXCR4","MX1","SAMD9L”,”IFI44L”,”XAF1”,”IFITM1","MS4A1","CD27","TNFRSF13B","CRIP1","ITGAX","FCRL5","FCRL3","IGHM","IGHD")
DotPlot(Bcells_intergrated, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()


t <- DotPlot( Bcells_intergrated, features=markers.to.plot, dot.scale=10, cols="RdBu", split.by="bsubtypes")

t<-VlnPlot(Bcells_integrated,features=markers.to.plot)



############################################crate data frame for barplot#################################
df<-data.frame(type=Bcells_intergrated$bsubtypes,ethnicity=Bcells_intergrated$ethnicity)
dittoBarPlot(Bcells_intergrated,"bsubtypes", group.by = "ethnicity",
                                             scale = "percent",legend.title='ethnicity',data.out=T)

saveRDS(Bcells_intergrated,"Bcells_intergrated.rds")
Bcells_intergrated<-ReadRDS("Bcells_intergrated.rds")
dittoBoxPlot(Bcells_integrated,)
plt<-dittoBarPlot(healthy_onlyy, "MCT", group.by = "race",data.out = T,
             scale = "count")
plt<-dittoBarPlot(Bcells_integrated, "bsubtypes", group.by = "ethnicity",data.out = T)#,
#scale = "count")
write.csv(plt$data,'proportionofMCT.csv')
healthy_onlyy <- subset (healthy_only, subset=race==c('chinese','malay','indian'))

#one sample test for prportions
anova(table(Bcells_intergrated$ethnicity))
table(Bcells_intergrated$ethnicity) #7822.667
length(Bcells_intergrated$ethnicity)
anova(Bcells_intergrated$ethnicity)
################################################################################################################
# get data
Bcells_intergrated<-readRDS("Bcells_intergrated.rds")
Chinese_vs_Malay <- subset (Bcells_intergrated, subset=ethnicity==c('Chinese','Malay'))
Malay_vs_Indian<-subset (Bcells_intergrated, subset=ethnicity==c('Malay','Indian' ))
Indian_vs_Chinese<-subset(Bcells_intergrated,subset=ethnicity==c('Indian','Chinese'))

# QC and filtering
# explore QC



DefaultAssay(Chinese_vs_Malay)<-'integrated'
# run Seurat's standard workflow steps
Chinese_vs_Malay<- NormalizeData(Chinese_vs_Malay)
Chinese_vs_Malay <- FindVariableFeatures(Chinese_vs_Malay)
Chinese_vs_Malay<- ScaleData(Chinese_vs_Malay)
Chinese_vs_Malay <- RunPCA(Chinese_vs_Malay)
ElbowPlot(Chinese_vs_Malay)
Chinese_vs_Malay <- RunUMAP(Chinese_vs_Malay, dims = 1:20)

# visualize 
Chinese_vs_Malay_final_plot <- DimPlot(Chinese_vs_Malay, reduction = 'umap', group.by = 'bsubtypes', label = TRUE)
Chinese_vs_Malay_final_ethnicity<- DimPlot(chinese_vs_malay, reduction = 'umap', group.by = 'ethnicity')

Chinese_vs_Malay_final_plot|Chinese_vs_Malay_final_ethnicity

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(healthy_only@meta.data)
Chinese_vs_Malay$ethnicity <- paste("",(Chinese_vs_Malay$ethnicity), sep = "~")
Chinese_vs_Malay$ethnicity <- paste((Chinese_vs_Malay$ethnicity),"", sep = "-")
Chinese_vs_Malay$Sample <- paste0(Chinese_vs_Malay$ethnicity, Chinese_vs_Malay$DCP_ID,Chinese_vs_Malay$lib)
Chinese_vs_Malay$Sample <- paste0(Chinese_vs_Malay$ethnicity, Chinese_vs_Malay$DCP_ID)
DefaultAssay(chinese_vs_malay)

cts <- AggregateExpression(Chinese_vs_Malay, 
                           group.by = c("bsubtypes", "Sample"),
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
  rownames(x) <- gsub('.*_~(.*)', '\\1', (rownames(x)))
  t(x)
  
  
})

#cts.split.modified <- lapply(cts.split, function(x){
 # rownames(x) <- substring(x, regexpr("_", x) + 1, nchar(x))
 # t(x)
#substring(x, regexpr("_", x) + 1, nchar(x))

#})
#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')

 
# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_MemoryB <- cts.split.modified$`Memory B`


# 2. generate sample level metadata
colData <- data.frame(Sample = colnames(counts_MemoryB))
coldata.df<-as.data.frame(colData)
coldata.df<- separate(coldata.df, col='Sample', into = c('ethni' ,'idb'), sep = '-')
colData$condition <-coldata.df$ethni




# get more information from metadata




# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_MemoryB,
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
res05 <- results(dds, name = "condition_Malay_vs_Chinese",alpha = 0.05)


summary(res05)

resLFC <- lfcShrink(dds, coef="condition_Malay_vs_Chinese", type="apeglm")

# filter
keep <- rowSums(counts(dds)) >=5
dds <- dds[keep,]

# run DESeq2
dds <- DESeq(dds,test="LRT",reduced = ~ 1,useT=TRUE,minmu=1e-6,  minReplicatesForReplace=Inf)
res05 <- results(dds, name = "condition_Malay_vs_Chinese",alpha = 0.05)
summary(res05)

rld <- rlog(dds, blind=TRUE)
rld <- vst(dds, blind=TRUE)

plotDispEsts(dds)
plotMA(dds)

# Plot PCA
DESeq2::plotPCA()
DESeq2::plotPCA(rld, intgroup = "condition",returnData=F)
p+geom_label(aes(label = colData$Sample))

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

plotCounts(dds, gene=which.min(res$padj), intgroup="condition",)
plotcoun


d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

###################################################################################MALAYVSINDIAN#########################################
Malay_vs_Indian<-subset (Bcells_intergrated, subset=ethnicity==c('Malay','Indian' ))
Indian_vs_Chinese<-subset(Bcells_intergrated,subset=ethnicity==c('Indian','Chinese'))

# QC and filtering
# explore QC



DefaultAssay(Chinese_vs_Malay)<-'integrated'
# run Seurat's standard workflow steps
Malay_vs_Indian<- NormalizeData(Malay_vs_Indian)
Malay_vs_Indian <- FindVariableFeatures(Malay_vs_Indian)
Malay_vs_Indian<- ScaleData(Malay_vs_Indian)
Malay_vs_Indian <- RunPCA(Malay_vs_Indian)
ElbowPlot(Malay_vs_Indian)
Malay_vs_Indian <- RunUMAP(Malay_vs_Indian, dims = 1:20)

# visualize 
Chinese_vs_Malay_final_plot <- DimPlot(Chinese_vs_Malay, reduction = 'umap', group.by = 'bsubtypes', label = TRUE)
Chinese_vs_Malay_final_ethnicity<- DimPlot(Chinese_vs_Malay, reduction = 'umap', group.by = 'ethnicity')

Chinese_vs_Malay_final_plot|Chinese_vs_Malay_final_ethnicity

# pseudo-bulk workflow -----------------
# Acquiring necessary metrics for aggregation across cells in a sample
# 1. counts matrix - sample level
# counts aggregate to sample level

View(healthy_only@meta.data)
Indian_vs_Chinese$ethnicity <- paste("",(Indian_vs_Chinese$ethnicity), sep = "~")
Indian_vs_Chinese$ethnicity <- paste((Indian_vs_Chinese$ethnicity),"", sep = "-")
Indian_vs_Chinese$Sample <- paste0(Indian_vs_Chinese$ethnicity, Indian_vs_Chinese$DCP_ID,Indian_vs_Chinese$lib)
Indian_vs_Chinese$Sample <- paste0(Indian_vs_Chinese$ethnicity, Indian_vs_Chinese$DCP_ID)
DefaultAssay(Indian_vs_Chinese)

cts <- AggregateExpression(Indian_vs_Chinese, 
                           group.by = c("bsubtypes", "Sample"),
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
  rownames(x) <- gsub('.*_~(.*)', '\\1', (rownames(x)))
  t(x)
  
  
})

#cts.split.modified <- lapply(cts.split, function(x){
# rownames(x) <- substring(x, regexpr("_", x) + 1, nchar(x))
# t(x)
#substring(x, regexpr("_", x) + 1, nchar(x))

#})
#gsub('.*_(.*)', '\\1', 'B cells_ctrl101')


# Let's run DE analysis with B cells
# 1. Get counts matrix
counts_MemoryB <- cts.split.modified$`Memory B`

write.csv(counts_MemoryB,"counts_memoryBchinesevsMalay.csv")


# 2. generate sample level metadata
colData <- data.frame(Sample = colnames(counts_AtypicalB))
coldata.df<-as.data.frame(colData)
coldata.df<- separate(coldata.df, col='Sample', into = c('ethni' ,'idb'), sep = '-')
colData$condition <-coldata.df$ethni

plotc









# perform DESeq2 --------
# Create DESeq2 object   
dds <- DESeqDataSetFromMatrix(countData = counts_AtypicalB,
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
res05 <- results(dds, name = "condition_Malay_vs_Indian",alpha = 0.05)


summary(res05)

resLFC <- lfcShrink(dds, coef="condition_Malay_vs_Indian", type="apeglm")

# filter
keep <- rowSums(counts(dds)) >=5
dds <- dds[keep,]

# run DESeq2 FOR SINGLE CELL AS SUGGESTED BY ARTICLES 
dds <- DESeq(dds,test="LRT",reduced = ~ 1,useT=TRUE,minmu=1e-6,  minReplicatesForReplace=Inf)
res05 <- results(dds, name = "condition_Indian_vs_Chinese",alpha = 0.05)
summary(res05)

rld <- rlog(dds, blind=TRUE)
rld <- vst(dds, blind=TRUE)

plotDispEsts(dds)
plotMA(dds)
?plotMA
# Plot PCA
DESeq2::plotPCA()
p<-DESeq2::plotPCA(rld, intgroup = "condition",returnData=F)
p+geom_label(aes(label = colData$Sample))

levels(factor(colData$condition))[2]
levels(colData$condition)[1]
contrast <- c("condition", levels(factor(colData$condition))[2], levels(factor(colData$condition))[1])
res <- results(dds, 
               contrast = contrast,
               alpha = 0.05)

res <- lfcShrink(dds, 
                 contrast =  contrast,
                 res=res)
summary(res05)
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
AtypicalB<-write.csv(sig_res,'AtypicalBchinese_indian.csv')

colData$condition


padj_cutoff <- 0.05

# Subset the significant results
sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff) %>%
  dplyr::arrange(padj)qa

# Check significant genes output
sig_res

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotcoun

plotCounts(dds, gene,intgroup="condition")

d <- plotCounts(dds, gene='MZT2A', intgroup="condition", 
                returnData=T)
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(2,4,6))

#SCATTERPLOT
## ggplot of top genes
normalized_counts <- counts(dds, 
                            normalized = TRUE)

## Order results by padj values
top20_sig_genes <- sig_res %>%
  dplyr::arrange(padj) %>%
  dplyr::pull(gene) %>%
  head(n=20)


top20_sig_norm <- data.frame(normalized_counts) %>%
  rownames_to_column(var = "gene") %>%
  dplyr::filter(gene %in% top20_sig_genes)

gathered_top20_sign <- top20_sig_norm %>%
  gather(colnames(top20_sig_norm)[2:length(colnames(top20_sig_norm))], key = "samplename", value = "normalized_counts")
gathered_top20_sig <- inner_join(colData[, c("Sample", "condition" )], gathered_top20_sig, by = c("Sample" = "samplename"))
gathered_top20_sign$condition<- separate(gathered_top20_sign, col='samplename', into = c('condition' ,'Sample_id'), sep = '.')
gathered_top20_sign<-NULL
gathered_top20_sig$Sample<-colData$Sample
gathered_top20_sig$condition<-colData$condition
gathered_top20_sig$gene<-top20_sig_norm$gene
gathered_top20_sig$normaized_counts<-vatop20_sig_norm$

## plot using ggplot2
ggplot(gathered_top20_sig) +
  geom_point(aes(x = gene, 
                 y = normalized_counts, 
                 color = condition), 
             position=position_jitter(w=0.1,h=0)) +
  scale_y_log10() +
  xlab("genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(plot.title = element_text(hjust = 0.5))








library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))





res_table_thres <- res_tbl %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58)

## Volcano plot
ggplot(res_table_thres) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  ggtitle("Volcano plot of naive B cells in CHINESE relative to Indian") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  scale_y_continuous(limits = c(-1,3)) +
  scale_x_continuous(limits = c(-20,20),breaks = c(-2,2,4)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))                