Single_R_function <- function(Seurat_Object,
                    Organism="Human",
                    Reduction="umap",
                    Legend_Dot_Size=2,
                    Perform_fine_tuning=FALSE,
                    Parallelize_in_Spark=FALSE,
                    Number_of_cells_per_partition=400,
                    Image_type="png",
                    image_out_dir
                    )
                    {
  #image: png
  
  imageType = "png"
  #Need to copy and paste into console to run:
  #BiocManager::install("GenomeInfoDbData", destdir="/scratch")
  
  suppressMessages(library(GenomeInfoDbData))
  suppressMessages(library(Seurat))
  suppressMessages(library(SingleR))
  suppressMessages(library(gridExtra))
  suppressMessages(library(tools))
  suppressMessages(library(grid))
  suppressMessages(library(gridBase))
  suppressMessages(library(cowplot))
  suppressMessages(library(ggplot2))
  suppressMessages(library(RColorBrewer))
  suppressMessages(library(magrittr))
  
  
  so <- Seurat_Object$value    
  species <- Organism
  doFineTuning <- Perform_fine_tuning
  useSpark <-Parallelize_in_Spark
  clusterList <-so@meta.data$seurat_clusters
  
  annotations <- function(so) {
    counts = GetAssayData(object = so)[, colnames(x = so)]
    if (species == "Human") {
      singler = CreateSinglerObject(counts = counts, project.name = "projectDesc", annot = NULL, min.genes = 0, technology = '10x',
                                    species = toTitleCase(species), citation = '', variable.genes = 'de', normalize.gene.length = F, fine.tune = doFineTuning,
                                    numCores = 4,reduce.file.size = T, clusters = so@meta.data$seurat_clusters
                                    , do.signatures = T
      )
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="HPCA_main")
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="annot")
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single$labels, col.name="HPCA")
      so <- AddMetaData(so, singler$singler[[2]]$SingleR.single.main$labels, col.name="BP_encode_main")
      so <- AddMetaData(so, singler$singler[[2]]$SingleR.single$labels, col.name="BP_encode")
    }
    if (species == "Mouse") {
      singler = CreateSinglerObject(counts = counts, project.name = "projectDesc" , annot = NULL, min.genes = 0, technology = '10x',
                                    species = toTitleCase(species), citation = '', variable.genes = 'de', normalize.gene.length = F, fine.tune = doFineTuning,
                                    numCores = 4, reduce.file.size = T ,clusters = so@meta.data$seurat_clusters
                                    , do.signatures = T
      )
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="immgen_main")
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="annot")
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single$labels, col.name="immgen")
      so <- AddMetaData(so, singler$singler[[2]]$SingleR.single.main$labels, col.name="mouseRNAseq_main")
      so <- AddMetaData(so, singler$singler[[2]]$SingleR.single$labels, col.name="mouseRNAseq")
      
    }
    return(so)
  }
  
  annotations_with_spark <- function(so) {
    N <- Number_of_cells_per_partition
    counts = GetAssayData(object = so)[, colnames(x = so)]
    rownames(counts) <- sub(".*?-", "", rownames(counts)) 
    n = ncol(counts)
    s = seq(1, n, by=N)
    doSingleRPerPartition <- function(i) {
      BiocManager::install("GenomeInfoDbData")
      library(SingleR)
      library(tools)
      A = seq(i,min(i+N-1,n))
      singler = CreateSinglerObject(counts = counts[,A], project.name = "projectDesc", annot = NULL, min.genes = 0, technology = '10x',
                                    species = toTitleCase(species), citation = '', variable.genes = 'de', normalize.gene.length = F, fine.tune = FALSE,
                                    numCores = 1,reduce.file.size = T, clusters = NULL, do.signatures = T
      )
    }
    if (species == "Human") {
      singler.objects <- spark.lapply(s, doSingleRPerPartition)
      singler = SingleR.Combine(singler.objects)
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="HPCA_main")
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="annot")
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single$labels, col.name="HPCA")
      so <- AddMetaData(so, singler$singler[[2]]$SingleR.single.main$labels, col.name="BP_encode_main")
      so <- AddMetaData(so, singler$singler[[2]]$SingleR.single$labels, col.name="BP_encode")
      
    }
    if (species == "Mouse") {
      singler.objects <- spark.lapply(s, doSingleRPerPartition)
      singler = SingleR.Combine(singler.objects)
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="immgen_main")
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single.main$labels, col.name="annot")
      so <- AddMetaData(so, singler$singler[[1]]$SingleR.single$labels, col.name="immgen")
      so <- AddMetaData(so, singler$singler[[2]]$SingleR.single.main$labels, col.name="mouseRNAseq_main")
      so <- AddMetaData(so, singler$singler[[2]]$SingleR.single$labels, col.name="mouseRNAseq")
      
    }
    
    return(so)
  }
  
  if (useSpark) {
    so <- annotations_with_spark(so)
    print("done")
  } else {
    so <- annotations(so)
    print("done")
  }
  
  if (species == "Human") {
    numColors = max(length(unique(so@meta.data$BP_encode_main)),length(unique(so@meta.data$HPCA_main)))
  } else {
    numColors = max(length(unique(so@meta.data$mouseRNAseq_main)),length(unique(so@meta.data$immgen_main)))
  }
  colpaired = colorRampPalette(brewer.pal(12,"Paired"))
  cols=c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4","#42d4f4","#f032e6","#bfef45","#fabebe","#469990","#e6beff","#9A6324","#800000","#aaffc3","#808000","#000075",colpaired(numColors))
  
  imageWidth = 5000
  imageHeight = 3000
  dpi = 300
  
  if (imageType == 'png') {
    png(
      filename=paste0(image_out_dir, "/", "Single_R_result.png"),
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")
  } else {
    library(svglite)
    svglite::svglite(
      filename=paste0(image_out_dir, "/", "Single_R_result.png"),
      width=round(imageWidth/dpi,digits=2),
      height=round(imageHeight/dpi,digits=2),
      pointsize=1,
      bg="white")
  }
  
  if (species =="Human") {
    p1 = DimPlot(so, reduction=Reduction, group.by="HPCA_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=Legend_Dot_Size),colour=guide_legend(ncol=4)) + ggtitle("HPCA Main Cell Type Annotations")
    p2 = DimPlot(so, reduction=Reduction, group.by="BP_encode_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=Legend_Dot_Size),colour=guide_legend(ncol=4))+ ggtitle("BP Encode Main Cell Type Annotations")
  } else {
    p1 = DimPlot(so, reduction=Reduction, group.by="immgen_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=Legend_Dot_Size),colour=guide_legend(ncol=4)) + ggtitle("Immgen Main Cell Type Annotations")
    p2 = DimPlot(so, reduction=Reduction, group.by="mouseRNAseq_main")  + scale_color_manual(values=cols) + theme(legend.position="top") + guides(override.aes = list(size=Legend_Dot_Size),colour=guide_legend(ncol=4))+ ggtitle("Mouse RNAseq Main Cell Type Annotations")
  }
  
  print(plot_grid(p1,p2,nrow=1))
  so@meta.data$Barcode <- rownames(so@meta.data)
  so@meta.data$sample_name <-so@meta.data$orig.ident
  so@meta.data$sample_name <- gsub("-","_",so@meta.data$sample_name)
  so@meta.data$annot <- NULL
  so@meta.data$seurat_clusters <- NULL
  rownames(so@meta.data) <-so@meta.data$Barcode
  return(list(value=so))
}
