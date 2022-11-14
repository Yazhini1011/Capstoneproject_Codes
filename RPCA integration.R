

batch.combined.final<- readRDS("C:\\Users\\yazhinir\\OneDrive - A STAR\\Documents\\for_integration\\metadatacompleted.rds")

# split the dataset into a list of two seurat objects (stim and CTRL)
batch.list <- SplitObject(merged.batches, split.by = "batch")



# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = batch.list)
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

##########################integration#######################################################
immune.anchors <- FindIntegrationAnchors(object.list = batch.list, anchor.features = features, reduction = "rpca")
# this command creates an 'integrated' data assay
batch.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(batch.combined) <- "integrated"


saveRDS(batch.combined, "RPCAintegrated.rds")

RPCA_integrated<- readRDS("C:\\Users\\yazhinir\\OneDrive - A STAR\\Desktop\\new_integration\\RPCAintegrated.rds")

###########################################################################performing lisi########################################################





install.packages('devtools')


DefaultAssay(RPCA_integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
RPCA_integrated  <- ScaleData(RPCA_integrated, verbose = FALSE)
RPCA_integrated  <- RunPCA(RPCA_integrated , npcs = 30, verbose = FALSE)
ElbowPlot(RPCA_integrated, ndims =30)
RPCA_integrated  <- RunUMAP(RPCA_integrated , reduction = "pca", dims = 1:18)
RPCA_integrated  <- FindNeighbors(RPCA_integrated, reduction = "pca", dims = 1:18)
RPCA_integrated  <- FindClusters(RPCA_integrated , resolution = 2)
RPCA_integrated  <- FindClusters(RPCA_integrated , resolution = 1)

AFTERRPCA <-DimPlot(RPCA_integrated, reduction = "umap", group.by = "batch")
BEFORE | AFTERRPCA
DimPlot(healthy_only, reduction = "umap", group.by = "batch")

umap.df.s <- FetchData(RPCA_integrated, vars = c("UMAP_1","UMAP_2"))
umap.df.s$barcode <- rownames(umap.df.s)
umap.df.s$barcode <- rownames(umap.df.s)

  


# Run the standard workflow for visualization and clustering
batch.combined <- ScaleData(batch.combined)
batch.combined <- RunPCA(batch.combined, npcs = 30)

pdf("1.Elbowplot.pdf")
ElbowPlot(batch.combined, ndims =30)
dev.off()


#UMAP
batch.combined <- RunUMAP(batch.combined, reduction = "pca", dims = 1:18)

# clustering
batch.combined <- FindNeighbors(batch.combined, dims = 1:18)
batch.combined <- FindClusters(batch.combined, resolution = 2)

umap.df <- FetchData(har, vars = c("UMAP_1","UMAP_2"))
umap.df$barcode <- rownames(umap.df)
write.table(umap.df,"umap.df.txt",col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
saveRDS(batch.combined,"batch.combined_final")

batch.combined <- FindNeighbors(object = batch.combined, 
                                 reduction = "pca", dims = 1:17, 
                                 k.param = 500, verbose = TRUE,compute.SNN = FALSE)
saveRDS(batch.combined,"batch.combined.after.QC.rds")
