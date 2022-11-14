install.packages("remotes")
library(remotes)
install_github("prabhakarlab/RCAv2")

library(RCAv2)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(DoubletFinder)
library(Seurat)
library(SingleR)
install.packages('BiocManager')
BiocManager::install('multtest')
install.packages('metap')

setwd("C:/Users/yazhinir/Downloads/library 1")


#Generate a RCA object from scRNA-seq data
PBMC_r <- createRCAObjectFrom10X("C:\\Users\\yazhinir\\Downloads\\library 1\\matrix")

# plot NODG distribution
PBMC_data <- PBMC_r$raw.data
nGeneVec <- Matrix::colSums(PBMC_data>0)

# Compute nUMI vector
nUMIVec <- Matrix::colSums(PBMC_data)

# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(PBMC_data), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(PBMC_data[mito.genes, , drop = FALSE])/Matrix::colSums(PBMC_data)

all_ngene <- data.frame(nGene = nGeneVec, 
                        nUMI = nUMIVec,
                        pMito = pMitoVec,
                        group=rep("SG_HEL_B001_L002", length(nGeneVec)))

pdf("1_distribution_of_NODG.pdf", width = 9, height = 3)
p1 <- ggplot(all_ngene, aes(x=nGene, y=pMito) ) + geom_point(size=.1)+
  geom_density_2d()
p2 <- ggplot(all_ngene, aes(x=nGene, y=nUMI) ) + geom_point(size=.1)+
  geom_density_2d()
p3 <- ggplot(all_ngene, aes(x=nUMI, y=pMito) ) + geom_point(size=.1)+
  geom_density_2d()
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

pdf("1_distribution_of_NODG_vs_pMito_zoom_in.pdf", width = 5, height = 5)
ggplot(all_ngene, aes(x=nGene, y=pMito) ) + geom_point(size=.1, color="grey")+
  geom_density_2d()+theme_bw(11)+xlim(0,6000)+ylim(0,0.25)
dev.off()

# select cells with at least 300 detected genes
PBMC_r <- dataFilter(PBMC_r,
                     nGene.thresholds = c(300,Inf),
                     nUMI.thresholds = c(0, Inf),
                     percent.mito.thresholds = c(0, 1),
                     min.cell.exp = 0.01*ncol(PBMC_r$raw.data), 
                     plot = FALSE)


all_ngene <- all_ngene[colnames(PBMC_r$raw.data),]

pdf("2_distribution_of_NODG_after_300NODG.pdf", width = 9, height = 3)
p1 <- ggplot(all_ngene, aes(x=nGene, y=pMito) ) + geom_point(size=.1)+
  geom_density_2d()
p2 <- ggplot(all_ngene, aes(x=nGene, y=nUMI) ) + geom_point(size=.1)+
  geom_density_2d()
p3 <- ggplot(all_ngene, aes(x=nUMI, y=pMito) ) + geom_point(size=.1)+
  geom_density_2d()
grid.arrange(p1, p2, p3, nrow = 1)
dev.off()

pdf("2_distribution_of_NODG_after_300NODG_zoom_in.pdf", width = 5, height = 5)
ggplot(all_ngene, aes(x=nGene, y=pMito) ) + geom_point(size=.1, color="green")+
  geom_density_2d()+theme_dark(11)+xlim(0,5000)+ylim(0,0.1)
dev.off()

# remove gene "^MT-|^ERCC|^RPS|^RPL|^HSP"
gene.row <- rownames(PBMC_r$raw.data)
gene.row.good <- gene.row[grep(pattern = "^MT-|^ERCC|^RPS|^RPL|^HSP", 
                               x = gene.row, invert = TRUE)]
# keep the raw data so that we could know what pMT is for RCA heatmap
raw.data.old <- PBMC_r$raw.data
PBMC_r$raw.data <- PBMC_r$raw.data[gene.row.good,]

# normalize data
PBMC_r <- dataLogNormalise(PBMC_r)

############ project to all immune panels ################
# get projection results against global panel
PBMC_r <- dataProject(PBMC_r, method = "GlobalPanel",
                      corMeth = "pearson", scale = TRUE)

global.proj <- as.data.frame(PBMC_r$projection.data)


global.proj.immune <- read.table("rownames_of_global_projection_immune_cells.txt", 
                                 header = FALSE)
global.proj <- global.proj[global.proj.immune$V1,]

# get projection results against other two panels
PBMC_r <- dataProjectMultiPanel(PBMC_r,method = list("NovershternPanel", 
                                                     "MonacoPanel"),
                                scale = TRUE,corMeth = "pearson")
two.proj <- as.data.frame(PBMC_r$projection.data)

# combine these panels
proj.all <- rbind(global.proj,two.proj)
proj.all <- as.matrix(proj.all)
proj.all <- as(proj.all, "dgCMatrix")

# Assign projection result to RCA object
PBMC_r$projection.data <- proj.all

#Estimate the most probable cell type label for each cell
PBMC_r <- estimateCellTypeFromProjection(PBMC_r,confidence = NULL)

###Skip the doublet removal (demuxlet and doublet finder steps) if using the practice dataset from 10X/Seurat
### doublet from demuxlet

demuxlet <- read.table("B021L001.scRNA.barcodeSummary.tsv",
                       sep = "\t", header = TRUE)

demuxlet <- demuxlet[which(demuxlet$Barcode %in% colnames(PBMC_r$raw.data)),] 

demuxlet$type<-unlist(lapply(demuxlet$SampleIdentity,
                             function(x)unlist(strsplit(x,split = ":"))[1]))

demuxlet.singlet <- demuxlet[which(demuxlet$type=="SNG"),"Barcode"]

# estimate doublet ratio
number_of_sample <- 16
# change this value for each batch
estimated.doublet.ratio <- (sum(demuxlet$type=="AMB"))+(sum(demuxlet$type=="MIX"))/length(demuxlet$type)*number_of_sample/(number_of_sample-1)
estimated.doublet.ratio <- ((sum(demuxlet$type=="AMB"))+(sum(demuxlet$type=="MIX")))/length(demuxlet$type)*number_of_sample/(number_of_sample-1)

################### Doublet Finder ##################################
# generate seurat object for doubletfinder
PBMC_s <- CreateSeuratObject(counts = PBMC_r$raw.data, 
                             project = "AIDA", 
                             min.cells = 0, 
                             min.features = 0)

# add ground truth by demuxlet
GT.df <- data.frame(barcode = colnames(PBMC_s))
rownames(GT.df) <- GT.df$barcode
GT.df$GT <- "Doublet"
GT.df[(GT.df$barcode %in% demuxlet.singlet),"GT"] <- "Singlet"
# add ground-truth doublet-singlet information
PBMC_s$GT <- GT.df[colnames(PBMC_s), "GT"]


# add RCA projection annotation
RCA.cell.df <- data.frame(barcode = colnames(PBMC_r$projection.data),
                          cell = unlist(PBMC_r$cell.Type.Estimate))

rownames(RCA.cell.df) <- RCA.cell.df$barcode
PBMC_s$RCA.proj.cell <- RCA.cell.df[colnames(PBMC_s),"cell"]

# normalize the data
PBMC_s <- NormalizeData(PBMC_s, normalization.method = "LogNormalize", scale.factor = 10000)

# select feature genes
PBMC_s <- FindVariableFeatures(PBMC_s, selection.method = "vst", nfeatures = 2000)

# scale the data
all.genes <- rownames(PBMC_s)
PBMC_s <- ScaleData(PBMC_s, features = all.genes)

# run PCA
PBMC_s <- RunPCA(PBMC_s, features = VariableFeatures(object = PBMC_s))
PBMC_s
ElbowPlot(PBMC_s, ndims = 30)

# UMAP
PBMC_s <- RunUMAP(PBMC_s, dims = 1:15)
PBMC_s
# clustering
PBMC_s <- FindNeighbors(PBMC_s, dims = 1:15)
PBMC_s <- FindClusters(PBMC_s, resolution = 1)

pdf("3.UMAP_of_seurat_for_DoubletFinder_without_annotation_new.pdf")
DimPlot(PBMC_s, reduction = "umap", label = TRUE)+NoLegend()
dev.off()

######### RCA annotation
RCA.results <- FetchData(PBMC_s, vars = c("seurat_clusters","RCA.proj.cell"))

RCA.results$count <- 1
RCA.results <- aggregate(count ~ ., RCA.results, FUN = sum)
RCA.results <- RCA.results %>% group_by(seurat_clusters) %>% top_n(n = 5, wt = count)
RCA.results <- RCA.results[order(RCA.results$seurat_clusters),]
write.csv(RCA.results, "RCA_results.csv")

# annotate each cluster
umap.df.cluster <- FetchData(PBMC_s, vars = c("UMAP_1","UMAP_2","seurat_clusters"))
# singler.results <- FetchData(eils, vars = c("seurat_clusters","singler"))
# singler.results$count <- 1
# singler.results <- aggregate(count ~ ., singler.results, FUN = sum)
# singler.final <- singler.results %>% group_by(seurat_clusters) %>% top_n(n = 5, wt = count)
# singler.final <- singler.final[order(singler.final$seurat_clusters),]

umap.df.cluster$cell <- NA
umap.df.cluster[which(umap.df.cluster$seurat_clusters==0),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==1),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==2),"cell"] <- "CD4 T"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==3),"cell"] <- "CD4 T"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==4),"cell"] <- "CD14 Mono"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==5),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==6),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==7),"cell"] <- "CD4 T"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==8),"cell"] <- "B"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==9),"cell"] <- "CD8 T"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==10),"cell"] <-"CD4 T"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==11),"cell"] <-"B"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==12),"cell"] <-"CD8 T"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==13),"cell"] <-"CD16 Mono"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==14),"cell"] <-"NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==15),"cell"] <-"CD8 T"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==16),"cell"] <-"CD8 T"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==17),"cell"] <-"B"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==18),"cell"] <- "CD14 Mono"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==19),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==20),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==21),"cell"] <- "CD14 Mono"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==22),"cell"] <- "Platelet"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==23),"cell"] <- "mDC"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==24),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==25),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==26),"cell"] <- "CD16 Mono"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==27),"cell"] <- "CD16 Mono"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==28),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==29),"cell"] <- "NK"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==30),"cell"] <- "CD14 Mono"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==31),"cell"] <- "pDC"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==32),"cell"] <- "stem cells"
umap.df.cluster[which(umap.df.cluster$seurat_clusters==33),"cell"] <- "pDC"




PBMC_s$RCA.annotation <- umap.df.cluster[colnames(PBMC_s),"cell"]

pdf("4.UMAP_of_seurat_for_DoubletFinder_with_annotation_from_singleR.pdf")
DimPlot(PBMC_s, reduction = "umap", group.by = "RCA.annotation", label = TRUE, label.size = 3)+
  NoLegend()
dev.off()


pdf("4.UMAP_of_seurat_for_DoubletFinder_with_marker_genes.pdf", width = 4*3, height = 4*3)

FeaturePlot(PBMC_s, features = c("CD14", "FCGR3A", "CD8A", "CD4", "CD3D", "NCAM1", "CD38", "CD16", "MS4A1", "PF4","CD27","TNFRSF17","CD1C","LILRA4"))





FeaturePlot(PBMC_s, features = "CD205")

FeaturePlot(PBMC_s, features = c("BTG1", "BTG4" ,"MZB1"))
dev.off()
######################to find marker genes###########mycode#############
cluster0_conserved_markers <- FindConservedMarkers(PBMC_s,
                                                   ident.1 = "cell",
                                                   grouping.var = "RCA.annotation",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
cluster0_conserved_markers <- FindConservedMarkers(PBMC_s,
                                                   ident.1 = "CD16",
                                                   grouping.var = "RCA.annotation",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)


## pK Identification (no ground-truth)
memory.limit(size=100000)
sweep.res <- paramSweep_v3(PBMC_s, PCs = 1:15, sct = FALSE)
saveRDS(sweep.res, "DoubletFinder.sweep.res.rds")


gt.calls <- PBMC_s@meta.data[rownames(sweep.res[[1]]), "GT"]
sweep.stats <- summarizeSweep(sweep.res, GT = TRUE, GT.calls = gt.calls)

pdf("12.bcmvn_doubletfinder.pdf")
bcmvn <- find.pK(sweep.stats)
dev.off()
#bcmvn in a console ----head(bcmvn)----choose the pk value where the bc matrix is the highest ----its usually 0.001 or 0.005-----and then cn change the value of the line (300,303)
 
## Homotypic Doublet Proportion Estimate
annotations <- PBMC_s$RCA.annotation
homotypic.prop <- modelHomotypic(annotations)# usually 0.2 but can be higher depending on the data 
nExp_poi <- round(estimated.doublet.ratio*length(colnames(PBMC_s)))# total number of expected doublets , the doublet finder finds out.
#####nexp hould be generlly about 8000
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))# lower than the total number of expected doublets 
PBMC_s$seurat_clusters
annotations <- PBMC_s$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(estimated.doublet.ratio*length(colnames(PBMC_s)))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## Run DoubletFinder with varying classification stringencies
PBMC_s <- doubletFinder_v3(PBMC_s, PCs = 1:15, pN = 0.25, 
                           pK = 0.005, nExp = nExp_poi, 
                           reuse.pANN = FALSE, sct = FALSE)
PBMC_s$pANN_0.25_0.005_13086

#########should be run in console , check the colnames(PBMC_S@meta.data)
PBMC_s <- doubletFinder_v3(PBMC_s, PCs = 1:15, pN = 0.25, 
                           pK = 0.005, nExp = nExp_poi.adj, 
                           reuse.pANN = "pANN_0.25_0.005_13086", sct = FALSE)
PBMC_s$DF.classifications_0.25_0.005_12208
PBMC_s$pANN_0.25_0.005_1503483
doublet.DF <- FetchData(PBMC_s, vars = c("UMAP_1","UMAP_2","DF.classifications_0.25_0.005_12208"))
#MINE

#choose the smaller value for doublet.df (line 309)

pdf("5.UMAP_of_seurat_for_DoubletFinder_singlet_doublet_from_DoubletFinder_only.pdf")
ggplot()+geom_point(data = doublet.DF[which(doublet.DF$DF.classifications_0.25_0.005_12208=="Singlet"),],
                    mapping = aes(x = UMAP_1, UMAP_2), size=.5,stroke=0, color="blue")+
  geom_point(data = doublet.DF[which(doublet.DF$DF.classifications_0.25_0.005_12208=="Doublet"),],
             mapping = aes(x = UMAP_1, UMAP_2), size=.5,stroke=0, color="red")

dev.off()

# final singlet after removing doublet finder
final.singlet.barcode <- rownames(doublet.DF[which(doublet.DF$DF.classifications_0.25_0.005_12208=="Singlet"),])
#########################MINE#############################################

doublet.DF$
  PBMC_s$

pdf("6.density_plot_of_pANN_between_all_and_singlet_only.pdf") #############use the higher DF.classification #############
plot(density(PBMC_s$pANN_0.25_0.005_13086), ylim=c(0,30), main="")
lines(density(PBMC_s$pANN_0.25_0.005_13086[final.singlet.barcode]), col="red")
dev.off()  



PBMC_s$pANN_0.25_0.005_1347448
final.singlet.barcode <- final.singlet.barcode[final.singlet.barcode %in% demuxlet.singlet]


doublet.df <- FetchData(PBMC_s, vars = c("UMAP_1","UMAP_2"))
doublet.df$barcode <- rownames(doublet.df)
doublet.df$final.doublet <- "Doublet"
doublet.df[which(doublet.df$barcode %in% final.singlet.barcode),"final.doublet"] <- "Singlet"

PBMC_s$final.doublet_type <- doublet.df[colnames(PBMC_s),"final.doublet"]

pdf("5.UMAP_of_seurat_for_final_singlets_and_doublets.pdf")
ggplot()+geom_point(data = doublet.df[which(doublet.df$final.doublet=="Singlet"),],
                    mapping = aes(x = UMAP_1, UMAP_2), size=.5,stroke=0, color="blue")+
  geom_point(data = doublet.df[which(doublet.df$final.doublet=="Doublet"),],
             mapping = aes(x = UMAP_1, UMAP_2), size=.5, stroke=0,color="red")
dev.off()

## hightlighting doublets identified by different tools
doublet.df <- FetchData(PBMC_s, vars = c("UMAP_1","UMAP_2","DF.classifications_0.25_0.005_12208"))#smaller value
doublet.df$barcode <- rownames(doublet.df)
doublet.df$demuxlet <- "Doublet"
doublet.df[which(doublet.df$barcode %in% demuxlet.singlet),"demuxlet"] <- "Singlet"


pdf("5.UMAP_of_seurat_for_highlighting_doublets_from_different_tools.pdf")#smaller value
ggplot()+geom_point(data = doublet.df[which(doublet.df$DF.classifications_0.25_0.005_12208=="Singlet" &
                                              doublet.df$demuxlet=="Singlet"),],
                    mapping = aes(x = UMAP_1, UMAP_2), size=.5,stroke=0, color="grey")+
  geom_point(data = doublet.df[which(doublet.df$DF.classifications_0.25_0.005_122088=="Doublet" &
                                       doublet.df$demuxlet=="Singlet"),],
             mapping = aes(x = UMAP_1, UMAP_2), size=.5,stroke=0, color="red")+
  geom_point(data = doublet.df[which(doublet.df$DF.classifications_0.25_0.005_12208=="Singlet" &
                                       doublet.df$demuxlet=="Doublet"),],
             mapping = aes(x = UMAP_1, UMAP_2), size=.5,stroke=0, color="orange")+
  geom_point(data = doublet.df[which(doublet.df$DF.classifications_0.25_0.005_12208=="Doublet" &
                                       doublet.df$demuxlet=="Doublet"),],
             mapping = aes(x = UMAP_1, UMAP_2), size=.5,stroke=0, color="blue")
dev.off()
##to add in numbers into the plot
dim(doublet.df[which(doublet.df$DF.classifications_0.25_0.005_12208=="Singlet" &
                       doublet.df$demuxlet=="Singlet"),])
dim(doublet.df[which(doublet.df$DF.classifications_0.25_0.005_12208=="Doublet" &
                       doublet.df$demuxlet=="Singlet"),])
dim(doublet.df[which(doublet.df$DF.classifications_0.25_0.005_12208=="Singlet" &
                       doublet.df$demuxlet=="Doublet"),])
dim(doublet.df[which(doublet.df$DF.classifications_0.25_0.005_12208=="Doublet" &
                       doublet.df$demuxlet=="Doublet"),])
HERE...
# plot QC for doublets from different tool
# doubletfinder
all_ngene.df <- all_ngene[doublet.df[which(doublet.df$DF.classifications_0.25_0.005_12208=="Doublet"),"barcode"],]

pdf("7.distribution_of_NODG_for_doublet_identified_by_doubletfinder_zoom_in.pdf", width = 6, height = 5)
ggplot(all_ngene.df, aes(x=nGene, y=pMito) ) + geom_point(size=.1)+
  geom_density_2d()+theme_bw(11)+xlim(0,6000)+ylim(0,0.1)
dev.off()

# demuxlet
all_ngene.free <- all_ngene[doublet.df[which(doublet.df$demuxlet=="Doublet"),"barcode"],]
pdf("7.distribution_of_NODG_for_doublet_identified_by_demuxlet_zoom_in.pdf", width = 6, height = 5)
ggplot(all_ngene.free, aes(x=nGene, y=pMito) ) + geom_point(size=.1)+
  geom_density_2d()+theme_bw(11)+xlim(0,6000)+ylim(0,0.1)
dev.off()


# plot QC metrics for doublet and singlet

all_ngene.new <- all_ngene#smaller
all_ngene.new$type <- "Doublet"
all_ngene.new[final.singlet.barcode,"type"] <- "Singlet"

pdf("7.distribution_of_NODG_for_singlet_and_doublet_zoom_in.pdf", width = 6, height = 5)
ggplot(all_ngene.new, aes(x=nGene, y=pMito, color=type) ) + geom_point(size=.2, stroke=0)+
  geom_density_2d()+theme_bw(11)+xlim(0,6000)+ylim(0,0.1)
dev.off()



saveRDS(PBMC_s, "Seurat_obj_for_DoubletFinder.rds")
PBMC_s <- readRDS("Seurat_obj_for_DoubletFinder.rds")


###start from here if not removing doublets
#### final singlet in RCA
raw.data.old <- raw.data.old[,final.singlet.barcode] # old raw data with Mt for RCA heatmap
PBMC_r.new <- createRCAObject(rawData = raw.data.old)

saveRDS(raw.data.old, "raw.data.old.rds")
raw.data.old <- readRDS("raw.data.old.rds")

saveRDS(PBMC_r.new, "RCA_PBMC_after_doublet_removal_and_before_QC.rds")
PBMC_r <- readRDS("RCA_PBMC_after_doublet_removal_and_before_QC.rds")

