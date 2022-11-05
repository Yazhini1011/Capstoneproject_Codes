####################################################################################
########################### profile ethnicity enrichment ###########################
####################################################################################

PBMC.integrated.final <- FindNeighbors(object = PBMC.integrated.final, 
                                       reduction = "pca", dims = 1:17, 
                                       k.param = 500, verbose = TRUE,compute.SNN = FALSE)

graph_nnMat_final <- Matrix::Matrix(PBMC.integrated.final@meta.data[,"integrated_nn"], sparse = TRUE)
saveRDS(graph_nnMat_final,"graph_nnMat_final.rds")
graph_nnMat_final <- readRDS("graph_nnMat_final.rds")

umap.df.final <- read.table("./umap.df.final.txt", header = TRUE, sep = "\t")


comp <- c("Japanese", "Korean","Malay","Indian","Caucasian")

ref <- "Chinese"
ref_num <- nrow(umap.df.final[(umap.df.final$ethnicity==ref),,drop=FALSE])
# calculate total number of ref cells in each cell's 500 neighbours
graph_nnMat_final.ref <- graph_nnMat_final[,umap.df.final[(umap.df.final$ethnicity==ref),"barcode"]]
graph_nnMat_final.ref.sum <- Matrix::rowSums(graph_nnMat_final.ref)


for (i in seq_along(comp)){
  comp.i <- comp[i]
  
  
  # calculate number of ref and comp cells in each cell's 300 neighbours
  graph_nnMat_final.comp <- graph_nnMat_final[, umap.df.final[(umap.df.final$ethnicity==comp.i),"barcode"] ]
  graph_nnMat_final.comp.sum <- Matrix::rowSums(graph_nnMat_final.comp)
  
  comp.ref.df <- data.frame(comp = graph_nnMat_final.comp.sum,
                            ref = graph_nnMat_final.ref.sum)
  
  comp_FC <- (nrow(umap.df.final[(umap.df.final$ethnicity==comp.i),,drop=FALSE])) / (ref_num)
  
  comp.ref.df$enrichment.ratio <-  (1 + comp.ref.df$comp) / (1 + comp.ref.df$ref)
  comp.ref.df$log2Ratio <- log2(comp.ref.df$enrichment.ratio/comp_FC)
  
  comp.ref.df[(comp.ref.df$comp == 0 & comp.ref.df$ref == 0),"log2Ratio"] <- 0
  
  comp.ref.df <- comp.ref.df[umap.df.final$barcode,]
  
  umap.df.final$log2Ratio <- comp.ref.df$log2Ratio
  umap.df.final[(umap.df.final$log2Ratio > 3),"log2Ratio"] <- 3
  umap.df.final[(umap.df.final$log2Ratio < -3),"log2Ratio"] <- -3
  
  ggplot()+
    geom_point(umap.df.final, 
               mapping = aes(x=UMAP_1, y=UMAP_2, 
                             color=log2Ratio), 
               size=0.2, stroke=0)+
    scale_color_gradientn(colors = c("darkorange","#e8ded3","#E8E8E8","#d5e0e8","darkblue" ),
                          values=c(1, .54, .5, .46, 0),
                          limits=c(-3,3))+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"))
  ggsave(paste0("3.Comparative_enrichment_plot_of_",
                comp.i,"_vs_",ref,"_discarded_one_JP_sample.png"),device = "png",
         width = 7, height = 6)
  
}

