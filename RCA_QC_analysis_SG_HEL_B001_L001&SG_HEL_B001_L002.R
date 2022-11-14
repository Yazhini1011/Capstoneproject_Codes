#set working directory
setwd("C:/Users/tanlm/Desktop/SG_HEL_B001")
library(Seurat)
library(RCAv2)
library(gridExtra)
library(ggplot2)
library(dplyr)
library(RColorBrewer) 
library(gsubfn)


lib_b1_l1 <- readRDS("C:/Users/tanlm/Desktop/SG_HEL_B001_L001_5GEX/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
lib_b1_l2 <- readRDS("C:/Users/tanlm/Desktop/SG_HEL_B001_L002_5GEX/RCA_PBMC_after_doublet_removal_and_before_QC.rds")
lib_b1_l1_raw_data <- lib_b1_l1$raw.data
lib_b1_l2_raw_data <- lib_b1_l2$raw.data
colnames(lib_b1_l1_raw_data)
colnames(lib_b1_l1_raw_data) <- paste(colnames(lib_b1_l1_raw_data), "_11", sep = "")
colnames(lib_b1_l2_raw_data) <- paste(colnames(lib_b1_l2_raw_data), "_12", sep = "")

common_gene <- intersect(rownames(lib_b1_l1_raw_data), rownames(lib_b1_l2_raw_data))
lib_b1_l1_raw_data<- lib_b1_l1_raw_data[common_gene,]
lib_b1_l2_raw_data <- lib_b1_l2_raw_data[common_gene,]

rca.all <- cbind(lib_b1_l1_raw_data,lib_b1_l2_raw_data)

############## using seurat to detect if variation exits between two libs ############
s.all <- CreateSeuratObject(counts = rca.all, min.cells = 0, min.features = 0)
# normalize the data
s.all <- NormalizeData(s.all, normalization.method = "LogNormalize", scale.factor = 10000)

# select feature genes
s.all <- FindVariableFeatures(s.all, selection.method = "vst", nfeatures = 2000)

# scale the data
all.genes <- rownames(s.all)
s.all <- ScaleData(s.all, features = all.genes)

# run PCA
s.all <- RunPCA(s.all, features = VariableFeatures(object = s.all))
ElbowPlot(s.all, ndims = 30)


# UMAP
s.all <- RunUMAP(s.all, dims = 1:20)


umap.df.s <- FetchData(s.all, vars = c("UMAP_1","UMAP_2"))
umap.df.s$barcode <- rownames(umap.df.s)
umap.df.s$lib <- "lib1"
umap.df.s[grep("_12",umap.df.s$barcode,value = TRUE),"lib"] <- "lib2"
umap.df.s$lib <- factor(umap.df.s$lib, levels = c("lib1","lib2"))

pdf("0.Seurat_UMAP_with_library_origin.pdf",width=9,height=7)
ggplot(data = umap.df.s, mapping = aes(x = UMAP_1, y = UMAP_2, color=lib)) + 
  geom_point(size = .5, stroke=0) + 
  scale_color_manual(values=c("blue",
                              "darkorange")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()


#If the UMAP show that the 2 libraries are well-mixed, proceed to do merging as shown below. If not, you must do batch effect correction.

lib_b1 <- SeuratObject:::RowMergeSparseMatrices(lib_b1_l1_raw_data, lib_b1_l2_raw_data)

PBMC_new <- createRCAObject(rawData = lib_b1)
# normalize data
PBMC_new <- dataLogNormalise(PBMC_new)

############ project to all immune panels ################
# get projection results against global panel
PBMC_new <- dataProject(PBMC_new, method = "GlobalPanel",
                      corMeth = "pearson", scale = TRUE)

global.proj <- as.data.frame(PBMC_new$projection.data)

global.proj.immune <- read.table("C:/Users/tanlm/Desktop/SG_HEL_B001/rownames_of_global_projection_immune_cells.txt", 
                                 header = FALSE)
global.proj <- global.proj[global.proj.immune$V1,]

# get projection results against other two panels
PBMC_new <- dataProjectMultiPanel(PBMC_new,method = list("NovershternPanel", 
                                                     "MonacoPanel"),
                                scale = TRUE,corMeth = "pearson")
two.proj <- as.data.frame(PBMC_new$projection.data)
# combine these panels
proj.all <- rbind(global.proj,two.proj)
proj.all <- as.matrix(proj.all)
proj.all <- as(proj.all, "dgCMatrix")
# Assign projection result to RCA object
PBMC_new$projection.data <- proj.all
#Estimate the most probable cell type label for each cell
PBMC_new <- estimateCellTypeFromProjection(PBMC_new,confidence = NULL)

# cluster cells
PBMC_new <- dataSClust(PBMC_new, res = 1)

# projection using UMAP
PBMC_new <- computeUMAP(PBMC_new)
plotRCAUMAP(PBMC_new,filename = "2_UMAP_PBMCs.pdf")

#estimate celltype projection 
PBMC_new <- estimateCellTypeFromProjection(PBMC_new,confidence = NULL)

color.cluster <- data.frame(cell=unlist(PBMC_new$cell.Type.Estimate),
                            color=PBMC_new$clustering.out$dynamicColorsList[[1]])
color.cluster$count <- 1
color.cluster <- aggregate(count ~ ., color.cluster, FUN = sum)
cluster.final <- color.cluster %>% group_by(color) %>% top_n(n = 5, wt = count)
write.csv(cluster.final, "cluster.final.csv")

#annotate 
umap.df <- PBMC_new$umap.coordinates
umap.df$barcode <- rownames(umap.df)
umap.df$cluster <- PBMC_new$clustering.out$dynamicColorsList[[1]]
umap.df$cell <- NA
umap.df[which(umap.df$cluster=="black"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="green"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "mDC"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "pDC"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "Platelet"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="red"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "T cells"


pdf("3.RCA_UMAP_of_cell_types_before_QC_filtering.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cell)) + 
  geom_point(size = .5) + 
  scale_color_manual(values=c("blue","grey","brown","darkgreen",
                              "darkorange","red", "pink", "turquoise", "yellow")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

# marker genes
marker_genes <- c("CD3D","CD3E","CD8A","CD4","NKG7","NCAM1","FCGR3A","GZMA", "GZMB",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CXCR2","CD27","CD38","TNFRSF17","ITGAM","LILRA4","FCER1A","IL3RA","CD1C",
                  "CD74","HLA-DQB1","HLA-DRA")
# CD1C mDC
# "CXCR2" neutrophil
# "CD27","CD38","TNFRSF17" plasma b
# LILRA4 pDC
# "FCER1A","CD1C" mDC

unique.cluster <- unique(umap.df$cluster)

# measure expression levels and percentages of marker genes across all clusters
marker.df <- data.frame()

for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(PBMC_new$data)) {
    marker.i <- PBMC_new$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cluster==cluster.j),"barcode"]
      marker.i.j <- marker.i[,barcode.j]
      exp_per.i.j <- sum(marker.i.j>0)*100/length(marker.i.j)
      exp_mean.i.j <- mean(marker.i.j)
      exp_per.i <- c(exp_per.i,exp_per.i.j)
      exp_mean.i <- c(exp_mean.i,exp_mean.i.j)
    }
    names(exp_per.i) <- unique.cluster
    
    exp_mean.i<- as.vector(scale(exp_mean.i,center = TRUE, scale = TRUE))
    names(exp_mean.i) <- unique.cluster
    marker.df.tmp <- data.frame(gene=rep(i,length(exp_mean.i)),
                                cluster = seq(1,length(exp_mean.i),1),
                                exp_per=exp_per.i,
                                exp_mean = exp_mean.i)
    
    marker.df <- rbind(marker.df,marker.df.tmp)
  }
  
}

marker.df$exp_per <- marker.df$exp_per/20
colnames(marker.df) <- c("gene","cluster","expresson_per","scaled_expression")

pdf("4.Marker_genes_across_different_clusters.pdf", width = 12, height = 6)
ggplot()+geom_point(data = marker.df, 
                    mapping = aes(x=cluster,y=gene,
                                  size=expresson_per, color=scaled_expression))+
  scale_colour_gradient(low = "white",high = "red")+
  scale_x_continuous(breaks=seq(1,length(unique.cluster),1),labels=unique.cluster)+
  scale_y_continuous(breaks=seq(1,length(marker_genes),1),labels=marker_genes)+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
dev.off()

saveRDS(PBMC_new, "RCA_PBMC_after_doublet_removal_and_before_QC.rds")
PBMC_new <- readRDS("RCA_PBMC_after_doublet_removal_and_before_QC.rds")

############# start QC
##### QC on each cluster before QC
# plot NODG distribution
PBMC_data <- PBMC_new$raw.data
nGeneVec <- Matrix::colSums(PBMC_data>0)
# Compute nUMI vector
nUMIVec <- Matrix::colSums(PBMC_data)
# Select mito genes
mito.genes = grep(pattern = "^MT-", x = rownames(PBMC_data), value = T)
# Compute percent.mito vector
pMitoVec <- Matrix::colSums(PBMC_data[mito.genes, , drop = FALSE])/Matrix::colSums(PBMC_data)

QC_list <- list()
unique_celltype <- unique(umap.df$cell)
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- umap.df[which(umap.df$cell==cell.i),"barcode"]
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  p.i <- ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = .3,stroke=0)  +
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()+ylim(0,0.2)
  QC_list[[i]] <- p.i
}
p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(unique_celltype))
pdf("5.QC_plot_of_each_cell_type_before_QC_filtering.pdf", 
    width = 3, height = 3*length(unique(umap.df$cell)))
grid.arrange(p.all,ncol=1)
dev.off()


##### QC filtering
nodg.low <- c(600,1000,1200,1200,1000,300,300,2000)
nodg.up <- c(3950,3000,2950,5100,2900,750,6000,4500)
pmt.low <- c(0.005,0.005,0.005,0.015,0.010,0,0,0.01)
pmt.up <- c(0.065,0.05,0.040,0.0625,0.0475,0.125,0.125,0.03)

good.barcode <- c()
QC_list <- list()
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- umap.df[which(umap.df$cell==cell.i),"barcode"]
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  
  nodg.low.i <- nodg.low[i]
  nodg.up.i <- nodg.up[i]
  pmt.low.i <- pmt.low[i]
  pmt.up.i <- pmt.up[i]
  qc.good.i <- qc.i[which(qc.i$NODG >= nodg.low.i & qc.i$NODG <= nodg.up.i & 
                            qc.i$pMt >= pmt.low.i & qc.i$pMt <= pmt.up.i),]
  good.barcode <- c(good.barcode, rownames(qc.good.i))
  
  
  qc.i$quality <- "bad"
  qc.i[rownames(qc.good.i),"quality"] <- "good"
  qc.i$quality <- factor(qc.i$quality, levels = c("bad","good"))
  p.i <- ggplot2::ggplot() +
    ggplot2::geom_point(data = qc.i, ggplot2::aes(x = NODG, y = pMt,color=quality),size = .3,stroke=0) +
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d(data = qc.i, mapping=aes(x = NODG, y = pMt)) +
    ggplot2::scale_color_manual(values = c("grey","black"))+ylim(0,0.2)
  QC_list[[i]] <- p.i
}

p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(QC_list))
pdf("6.QC_plot_of_each_cell_type_before_QC_filtering_wih_target_areas.pdf", 
    width = 4, height = 3*length(unique(umap.df$cell)))
grid.arrange(p.all,ncol=1)
dev.off()

good.barcode.index <- which(colnames(PBMC_new$raw.data) %in% good.barcode)
PBMC_new$raw.data <- PBMC_new$raw.data[, good.barcode.index]
PBMC_new$data <- PBMC_new$data[, good.barcode.index]
PBMC_new$projection.data <- PBMC_new$projection.data[, good.barcode.index]
PBMC_new$cell.Type.Estimate <- PBMC_new$cell.Type.Estimate[good.barcode.index]

PBMC_new$clustering.out$dynamicColorsList[[1]] <- PBMC_new$clustering.out$dynamicColorsList[[1]][good.barcode.index]


# projection using UMAP
PBMC_new <- computeUMAP(PBMC_new)

umap.df <- PBMC_new$umap.coordinates
umap.df$barcode <- rownames(umap.df)
umap.df$cluster <- PBMC_new$clustering.out$dynamicColorsList[[1]]
umap.df$cell <- NA
umap.df[which(umap.df$cluster=="black"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="blue"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="brown"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="cyan"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="darkgreen"),"cell"] <- "B cells"
umap.df[which(umap.df$cluster=="darkgrey"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkred"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="darkturquoise"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="green"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="greenyellow"),"cell"] <- "mDC"
umap.df[which(umap.df$cluster=="grey60"),"cell"] <- "Plasma B"
umap.df[which(umap.df$cluster=="lightcyan"),"cell"] <- "pDC"
umap.df[which(umap.df$cluster=="lightgreen"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="lightyellow"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="magenta"),"cell"] <- "Monocytes"
umap.df[which(umap.df$cluster=="midnightblue"),"cell"] <- "Platelet"
umap.df[which(umap.df$cluster=="pink"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="purple"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="red"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="royalblue"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="salmon"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="tan"),"cell"] <- "NK cells"
umap.df[which(umap.df$cluster=="turquoise"),"cell"] <- "T cells"
umap.df[which(umap.df$cluster=="yellow"),"cell"] <- "T cells"

pdf("7.RCA_UMAP_of_cell_types_after_QC.pdf",width=9,height=7)
ggplot(data = umap.df, mapping = aes(x = UMAP1, y = UMAP2, color=cell)) + 
  geom_point(size = .5) + 
  scale_color_manual(values=c("blue","brown","darkgreen",
                              "darkorange","steelblue","red", "pink", "turquoise")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()


##### QC on each cluster AFTER QC
QC_list <- list()
# unique_celltype <- unique(umap.df$cell)
for (i in seq_along(unique_celltype)){
  cell.i <- unique_celltype[i]
  barcode.i <- umap.df[which(umap.df$cell==cell.i),"barcode"]
  qc.i <- data.frame(NODG=nGeneVec[barcode.i],
                     pMt = pMitoVec[barcode.i])
  p.i <- ggplot2::ggplot(data = qc.i, ggplot2::aes(x = NODG, y = pMt)) +
    ggplot2::geom_point(size = .3, stroke=0) +
    ggplot2::theme_bw() + ggplot2::ggtitle(paste0(cell.i," n=",length(barcode.i)))+
    ggplot2::stat_density_2d()
  QC_list[[i]] <- p.i
}
p.all <- arrangeGrob(grobs=QC_list,
                     nrow=length(unique_celltype))

pdf("8.QC_plot_of_each_cell_type_after_QC_filtering.pdf", 
    width = 3, height = 3*length(unique(umap.df$cell)))
grid.arrange(p.all,ncol=1)
dev.off()


marker_genes <- c("CD3D","CD3E","CD8A","CD4","NCAM1","FCGR3A",
                  "CD14","MS4A1","ITGA2B","HBB",
                  "CD27","CD38","TNFRSF17","ITGAM","CD1C","IL3RA")


# CD74 HLA-DQB1 HLA-DRA DC
# "ELANE", "LTF","MMP8" myelocyte

marker.df <- data.frame()
unique.cluster <- unique(umap.df$cell)
for (i in seq_along(marker_genes)){
  gene_i <- marker_genes[i]
  if (gene_i %in% rownames(PBMC_new$data)){
    marker.i <- PBMC_new$data[gene_i,,drop=FALSE] 
    exp_per.i <- c()
    exp_mean.i <- c()
    for (j in seq_along(unique.cluster)){
      cluster.j <- unique.cluster[j]
      barcode.j <- umap.df[which(umap.df$cell==cluster.j),"barcode"]
      marker.i.j <- marker.i[,barcode.j]
      exp_per.i.j <- sum(marker.i.j>0)*100/length(marker.i.j)
      exp_mean.i.j <- mean(marker.i.j)
      exp_per.i <- c(exp_per.i,exp_per.i.j)
      exp_mean.i <- c(exp_mean.i,exp_mean.i.j)
    }
    names(exp_per.i) <- unique.cluster
    
    exp_mean.i<- as.vector(scale(exp_mean.i,center = TRUE, scale = TRUE))
    names(exp_mean.i) <- unique.cluster
    marker.df.tmp <- data.frame(gene=rep(i,length(exp_mean.i)),
                                cluster = seq(1,length(exp_mean.i),1),
                                exp_per=exp_per.i,
                                exp_mean = exp_mean.i)
    
    marker.df <- rbind(marker.df,marker.df.tmp)
  }
  
}

marker.df$exp_per <- marker.df$exp_per/20

pdf("9.Marker_genes_across_different_cells_after_QC.pdf", width = 12, height = 5)
ggplot()+geom_point(data = marker.df, 
                    mapping = aes(x=cluster,y=gene,
                                  size=exp_per, color=exp_mean))+
  scale_colour_gradient(low = "white",high = "red")+
  scale_x_continuous(breaks=seq(1,length(unique.cluster),1),labels=unique.cluster)+
  scale_y_continuous(breaks=seq(1,length(marker_genes),1),labels=marker_genes)+
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
dev.off()

saveRDS(PBMC_new, "RCA_PBMC_after_doublet_removal_bad_cell_removal_and_after_cluster_QC.rds")
PBMC_new <- readRDS("RCA_PBMC_after_doublet_removal_bad_cell_removal_and_after_cluster_QC.rds")



umap.df <- PBMC_new$umap.coordinates
umap.df$barcode <- rownames(umap.df)

demuxlet_1 <- read.table("C:/Users/tanlm/Desktop/SG_HEL_B001_L001_5GEX/AIDA-0007.scRNA.AIDA_HEL_B001_L001_demuxlet.missingSampleAdded.best",
                       sep = "\t", header = TRUE)

demuxlet_1$BARCODE <- paste(demuxlet_1$BARCODE, "_11", sep = "")

demuxlet_1 <- demuxlet_1[which(demuxlet_1$BARCODE %in% colnames(PBMC_new$raw.data)),]

demuxlet.singlet_1 <- demuxlet_1[which(demuxlet_1$DROPLET.TYPE=="SNG"),"BARCODE"]

#Demuxlet.singlet_1: 15948 singlets

demuxlet_2 <- read.table("C:/Users/tanlm/Desktop/SG_HEL_B001_L002_5GEX/AIDA-0007.scRNA.AIDA_HEL_B001_L002_demuxlet.missingSampleAdded.best",
                       sep = "\t", header = TRUE)

demuxlet_2$BARCODE <- paste(demuxlet_2$BARCODE, "_12", sep = "")

demuxlet_2 <- demuxlet_2[which(demuxlet_2$BARCODE %in% colnames(PBMC_new$raw.data)),]

demuxlet_2$type <- unlist(lapply(demuxlet_2$BEST, 
                               function(x) unlist(strsplit(x,split = "-"))[1] ))

demuxlet.singlet_2 <- demuxlet_2[which(demuxlet_2$type=="SNG"),"BARCODE"]

###--------------------------------------------------------------------------------------
###Le Min own addition of codes
demuxlet_1_LM<-demuxlet_1[,c("BARCODE","BEST.GUESS","DROPLET.TYPE")]
demuxlet_2_LM<-demuxlet_2[,c("BARCODE","BEST","type")]
names(demuxlet_1_LM)[2]<-"BEST"
names(demuxlet_1_LM)[3]<-"type"
demuxlet_1_LM<-demuxlet_1_LM[which(demuxlet_1_LM$type=="SNG"),]

to_replace <- list("24_2414_L3038016_01_GlobalScreening_iA1_1,24_2414_L3038016_01_GlobalScreening_iA1_1,0.00"="1",
                   "1_2414_AIDA-0002-1_01_GlobalScreening_iA1_1,1_2414_AIDA-0002-1_01_GlobalScreening_iA1_1,0.00"="2",
                   "3_2414_AIDA-0003_01_GlobalScreening_iA1_1,3_2414_AIDA-0003_01_GlobalScreening_iA1_1,0.00"="3",
                   "4_2414_AIDA-0004_01_GlobalScreening_iA1_1,4_2414_AIDA-0004_01_GlobalScreening_iA1_1,0.00"="4",
                   "5_2414_AIDA-0005-1_01_GlobalScreening_iA1_1,5_2414_AIDA-0005-1_01_GlobalScreening_iA1_1,0.00"="5",
                   "7_2414_AIDA-0006_01_GlobalScreening_iA1_1,7_2414_AIDA-0006_01_GlobalScreening_iA1_1,0.00"="6",
                   "SNG-AIDA-0007"="7",
                   "9_2414_AIDA-0008-1_01_GlobalScreening_iA1_1,9_2414_AIDA-0008-1_01_GlobalScreening_iA1_1,0.00"="8",
                   "11_2414_AIDA-0009-2-1_01_GlobalScreening_iA1_1,11_2414_AIDA-0009-2-1_01_GlobalScreening_iA1_1,0.00"="9",
                   "13_2414_AIDA-0010_01_GlobalScreening_iA1_1,13_2414_AIDA-0010_01_GlobalScreening_iA1_1,0.00"="10",
                   "14_2414_AIDA-0011-1_01_GlobalScreening_iA1_1,14_2414_AIDA-0011-1_01_GlobalScreening_iA1_1,0.00"="11",
                   "16_2414_AIDA-0012-1_01_GlobalScreening_iA1_1,16_2414_AIDA-0012-1_01_GlobalScreening_iA1_1,0.00"="12",
                   "18_2414_AIDA-0013_01_GlobalScreening_iA1_1,18_2414_AIDA-0013_01_GlobalScreening_iA1_1,0.00"="13",
                   "19_2414_AIDA-0014-2-1_01_GlobalScreening_iA1_1,19_2414_AIDA-0014-2-1_01_GlobalScreening_iA1_1,0.00"="14",
                   "21_2414_AIDA-0015_01_GlobalScreening_iA1_1,21_2414_AIDA-0015_01_GlobalScreening_iA1_1,0.00"="15",
                   "22_2414_AIDA-0016-1_01_GlobalScreening_iA1_1,22_2414_AIDA-0016-1_01_GlobalScreening_iA1_1,0.00"="16")

demuxlet_1_LM$BEST <- gsubfn(paste(names(to_replace),collapse="|"),to_replace,demuxlet_1_LM$BEST)

to_replace <- list("SNG-24_2414_L3038016_01_GlobalScreening_iA1_1"="1",
                   "SNG-1_2414_AIDA-0002-1_01_GlobalScreening_iA1_1"="2",
                   "SNG-3_2414_AIDA-0003_01_GlobalScreening_iA1_1"="3",
                   "SNG-4_2414_AIDA-0004_01_GlobalScreening_iA1_1"="4",
                   "SNG-5_2414_AIDA-0005-1_01_GlobalScreening_iA1_1"="5",
                   "SNG-7_2414_AIDA-0006_01_GlobalScreening_iA1_1"="6",
                   "SNG-AIDA-0007"="7",
                   "SNG-9_2414_AIDA-0008-1_01_GlobalScreening_iA1_1"="8",
                   "SNG-11_2414_AIDA-0009-2-1_01_GlobalScreening_iA1_1"="9",
                   "SNG-13_2414_AIDA-0010_01_GlobalScreening_iA1_1"="10",
                   "SNG-14_2414_AIDA-0011-1_01_GlobalScreening_iA1_1"="11",
                   "SNG-16_2414_AIDA-0012-1_01_GlobalScreening_iA1_1"="12",
                   "SNG-18_2414_AIDA-0013_01_GlobalScreening_iA1_1"="13",
                   "SNG-19_2414_AIDA-0014-2-1_01_GlobalScreening_iA1_1"="14",
                   "SNG-21_2414_AIDA-0015_01_GlobalScreening_iA1_1"="15",
                   "SNG-22_2414_AIDA-0016-1_01_GlobalScreening_iA1_1"="16")

demuxlet_2_LM$BEST <- gsubfn(paste(names(to_replace),collapse="|"),to_replace,demuxlet_2_LM$BEST)

demuxlet <-rbind(demuxlet_1_LM,demuxlet_2_LM)

demuxlet.subset <- demuxlet[which(demuxlet$BARCODE %in% umap.df$barcode), c("BARCODE", "BEST") ]
umap.df.merge <- merge(x=umap.df, y=demuxlet.subset, by.x = "barcode", by.y = "BARCODE")
length(unique(umap.df.merge$BEST))
unique_samples <- unique(umap.df.merge$BEST)

pdf("10.RCA_UMAP_of_cell_types_by_sample.pdf",width=9,height=7)
ggplot(data = umap.df.merge, mapping = aes(x = UMAP1, y = UMAP2, color=BEST)) + 
  geom_point(size = .5, stroke = 0) + 
  scale_color_manual(values=c("darkviolet",  "forestgreen","royalblue", "red","yellow","darkgrey","black","purple","tan","turquoise","pink","lightgreen","magenta","midnightblue","cyan","darkred")) +
  theme_bw(10) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()

###--------------------------------------------------------------------------------------

for (i in seq_along(unique_samples)){
  sample.i <- unique_samples[i]
  umap.df.merge.i <- umap.df.merge[which(umap.df.merge$BEST==sample.i),]
  ggplot() +
    geom_point(data=umap.df.merge, aes(x=UMAP1, y=UMAP2), size = .5, color= "grey") + 
    geom_point(data = umap.df.merge.i, aes(x=UMAP1, y=UMAP2), size = .5, color= "red") +
    ggsave(paste0("UMAP_for_sample_", sample.i, ".pdf", sep=""))
}

#to get number of samples from each individual
table(umap.df.merge$BEST)



