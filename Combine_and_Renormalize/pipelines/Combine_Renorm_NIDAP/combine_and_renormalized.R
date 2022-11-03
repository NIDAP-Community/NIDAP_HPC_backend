#' combine_and_renormalized function
#'
#' Developed by CCBR, a NIDAP template tranformed to function
#' Combines samples, rescales and renormalizes, runs Dimensional Reduction, and returns a combined Seurat Object. This template will summarize the multi-dimensionality of your data into a set of "principal components" to allow for easier analysis. There is an option to use this template to perform Integration, as well. This is Step 4 in the canonical Single Cell pipeline after elbow plots and regression.
#' 
#' @param Seurat_dataset_object Input dataset, should be a seurat dataset or a R foundry object (on NIDAP)
#' @param Number_of_components Number, default 15. Select the number of principal components for your analysis. Please see the elbow plot in the previous template to figure out what number of PCs explains your variance cut-off. For example, if the elbow plot has point at (15,0.02), it means that 15 PCs encapsulate 98% of the variance in your data.
#' @param Variables_to_regress Text, defauit "". should be one of (percent.mt, nCount_RNA, S.Score, G2M.Score, CC.Difference), Subtract (‘regress out’) this source of heterogeneity from the data. For example, to regress out mitochondrial effects, input "percent.mt."
#' @param Integrate_Data Boolean, TRUE/FALSE, default FALSE. Perform integration of cells across conditions using the most variant genes to identify cells most similar to each other.  N.B. Always look at cells before deciding whether to perform integration, and seek advice from bioinformatician.
#' @param Cluster_Resolution_low_range Number, default 0.2. Select minimum resolution for clustering plots. The lower you set this, the FEWER clusters will be generated.
#' @param Cluster_Resolution_high_range Number, default 1.2. Select the maximum resolution for clustering. The higher you set this number, the MORE clusters you will produced.
#' @param Cluster_Resolution_range_bins Number, default 0.2. Select the bins for your cluster plots. For example, if you input 0.2 as your bin, and have low/high resolution ranges of 0.2 and 0.6, then the template will produce cluster plots at resolutions of 0.2, 0.4 and 0.6.
#' @param Conserve_memory Boolean, default FALSE. If dataset is larger than ~40k filtered cells, toggle to TRUE. If TRUE, only variable genes will be available for downstream analysis. Default is FALSE.
#' @param Draw_UMAP Boolean, default TRUE. If TRUE, draw UMAP plot.
#' @param Draw_t_SNE Boolean, default TRUE. If TRUE, draw TSNE plot.
#' @param Image_type Text, default "png". Should be png or svg.
#' @param Cell_Hashing_Data Boolean, default FALSE. Toggle "true" if you are using cell-hashed data.
#' @param Project_Name Text, default "scRNAProject". 
#' @param Stop_metadata_from_merging Boolean, default FALSE. Toggle true to stop metadata from merging. Not recommended for standard pipeline.
#' @param seed_for_PCA Number, default 42.
#' @param seed_for_TSNE Number, default 1.
#' @param seed_for_UMAP Number, default 42.
#' @param Run_SCTransform Boolean, default TRUE. Set to TRUE to run SCTransform (recommended current v3 Seurat default). Set to FALSE to run ScaleData (previous v2 Seurat default) instead.
#' @param Exclude_sample Number, default 0. Exclude unwanted samples from the merge step. The number will correspond to the order in which they appear. Leave as 0 if you want to use all samples. If you want to exclude one or several samples, separate each sample number by comma (e.g. 1,2,3,4).
#' @param Number_of_Features Number, default 2000. Number of variable features
#' @param Mean_low_cutoff Number, default 0.1.
#' @param Mean_high_cutoff Number, default 8.
#' @param Dispersion_low_cutoff Number, default 1.
#' @param Dispersion_high_cutoff Number, default 100000. 
#' @param Selection_Method Text, default vst. Should be one of (vst, mean.var.plot, dispersion). Method to choose top variable features.
#' @image_out_dir Text, default at current working directory.
#' 
#' @import Seurat
#' @import ggplot2
#' @import gridExtra
#' @import RColorBrewer
#' @import svglite
#' @import plotly
#' @import jsonlite

#' @return Combined and renormalized dataset
#'
#' @examples
#' combine_and_renormalized(Seurat_dataset_object = dataset)
#'
#'
#'
#' @export
#' 



# combine_and_renormalized function from NIDAP


combine_and_renormalized <- function(Seurat_dataset_object, 
                    Number_of_components, 
                    Variables_to_regress, 
                    Integrate_Data, 
                    Cluster_Resolution_low_range, 
                    Cluster_Resolution_high_range, 
                    Cluster_Resolution_range_bins, 
                    Conserve_memory, 
                    Draw_UMAP, 
                    Draw_t_SNE, 
                    Image_type, 
                    Cell_Hashing_Data, 
                    Project_Name, 
                    Stop_metadata_from_merging, 
                    seed_for_PCA, 
                    seed_for_TSNE, 
                    seed_for_UMAP, 
                    Run_SCTransform, 
                    Exclude_sample, 
                    Number_of_Features, 
                    Mean_low_cutoff, 
                    Mean_high_cutoff, 
                    Dispersion_low_cutoff, 
                    Dispersion_high_cutoff, 
                    Selection_Method,
                    image_out_dir){
  
    if(missing(Seurat_dataset_object)){ stop("No input dataset.") }
  
    if(missing(Number_of_components)){ Number_of_components <- 15 }
    number_of_pc_s_for_pca_umap_t_sne <- Number_of_components
    
    if(missing(Variables_to_regress)){ Variables_to_regress <- "" }
    
    if(missing(Integrate_Data)){ Integrate_Data <- FALSE }
    
    if(missing(Cluster_Resolution_low_range)){ Cluster_Resolution_low_range <- 0.2 }
    cluster_resolution_low_range <- Cluster_Resolution_low_range
    
    if(missing(Cluster_Resolution_high_range)){ Cluster_Resolution_high_range <- 1.2 }
    cluster_resolution_high_range <- Cluster_Resolution_high_range
    
    if(missing(Cluster_Resolution_range_bins)){ Cluster_Resolution_range_bins <- 0.2 }
    cluster_resolution_range_bins <- Cluster_Resolution_range_bins
    
    if(missing(Conserve_memory)){ Conserve_memory <- FALSE }
    Conserve_memory <- Conserve_memory
    
    if(missing(Draw_UMAP)){ Draw_UMAP <- TRUE }
    draw_umap <- Draw_UMAP
    
    if(missing(Draw_t_SNE)){ Draw_t_SNE <- TRUE }
    draw_t_sne <- Draw_t_SNE
    
    if(missing(Image_type)){ Image_type <- "png"}
    
    if(missing(Cell_Hashing_Data)){ Cell_Hashing_Data <- FALSE }
    cell_hashing_data <- Cell_Hashing_Data
    
    if(missing(Project_Name)){ Project_Name <- "scRNAProject" }
    project_name <- Project_Name
    
    if(missing(Stop_metadata_from_merging)){ Stop_metadata_from_merging <- FALSE}
    Stop_metadata_from_merging_ <- Stop_metadata_from_merging
    
    if(missing(seed_for_PCA)){ seed_for_PCA <- 42 }
    
    if(missing(seed_for_TSNE)){ seed_for_TSNE <- 1 }
    
    if(missing(seed_for_UMAP)){ seed_for_UMAP <- 42 }
    
    if(missing(Run_SCTransform)){ Run_SCTransform <- TRUE }
    
    if(missing(Exclude_sample)){ Exclude_sample <- 0}
    
    if(missing(Number_of_Features)){ Number_of_Features <- 2000 }
    
    if(missing(Mean_low_cutoff)){ Mean_low_cutoff <- 0.1 }
    
    if(missing(Mean_high_cutoff)){ Mean_high_cutoff <- 8 }
    
    if(missing(Dispersion_low_cutoff)){ Dispersion_low_cutoff <- 1 }
    
    if(missing(Dispersion_high_cutoff)){ Dispersion_high_cutoff <- 100000 }
    
    if(missing(Selection_Method)){ Selection_Method <- "vst" }
    
    if(missing(image_out_dir)){ image_out_dir <- getwd() }
    
    
    suppressMessages(library(Seurat))
    suppressMessages(library(ggplot2))
    suppressMessages(library(gridExtra))
    suppressMessages(library(RColorBrewer))
    suppressMessages(library(svglite))
    suppressMessages(library(plotly))
    suppressMessages(library(jsonlite))
    
    
    
    if (class(Seurat_dataset_object) == "RFoundryObject"){  
      SO = Seurat_dataset_object$value
      
    }else if(class(Seurat_dataset_object) == "list"){
      SO = Seurat_dataset_object
      
    }else if(class(Seurat_dataset_object) == "Seurat"){
      SO = Seurat_dataset_object
      
    }else{
      stop("Input dataset should be a Seurat object")
      
      }
    

    ## Add commentary on how this toggle works with add.only.var.genes.
    conserve_memory <- Conserve_memory

    # If exclude option is TRUE, filter out undesirable sample
    if (c(0) == 0){
        SO <- SO
    } else {
        SO <- SO[-c(0)]
    }

    #initialize Citeseq functionality as false, 
    #later the template will check for a Protein assay and run if it finds it
    doCiteSeq <- FALSE

    doMergeData <- !FALSE
    dat = vector()
    integratedata = FALSE

    if (length(SO) > 1) {
    for(i in 2:length(SO)){dat=c(dat,SO[[i]]) }
    SO_merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = "scRNAProject", merge.data = TRUE)
    allgenes <- rownames(SO_merge)
    #SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
    } else {
    SO_merge <- SO[[1]]
    allgenes <- rownames(SO_merge)
    #SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)
    }

    if (!("orig.ident" %in% colnames(SO_merge@meta.data))) {
        SO_merge@meta.data$orig.ident <- SO_merge@meta.data$orig_ident
    }

    if ("Protein" %in% names(SO_merge@assays)){
        doCiteSeq <-TRUE
    }

    if(FALSE){
    #SO_merge <- ScaleData(SO_merge, assay = "HTO")
    }
    
    npcs = 15
    Do_SCTransform = TRUE
    vars_to_regress = c()

if (Do_SCTransform){
     if(is.null(vars_to_regress)){
        SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE, conserve.memory = conserve_memory, return.only.var.genes = FALSE)}
     else{       
        SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE,vars.to.regress=vars_to_regress, conserve.memory = conserve_memory, return.only.var.genes = FALSE) 
}
}
else{
    all.genes <- rownames(SO_merge)
    if(is.null(vars_to_regress)){
        SO_merge <- SO_merge
    }
    else{
        #SO_merge <- ScaleData(SO_merge, features=all.genes, assay = "RNA", vars.to.regress=vars_to_regress) 
    }
 DefaultAssay(SO_merge) <- "RNA"   
}

if (length(SO)>1) {
        all_features <- lapply(SO, row.names) %>% Reduce(intersect, .)
    if(integratedata==TRUE){
            integ_features <- SelectIntegrationFeatures(object.list = SO, nfeatures = 3000) 
            if(!is.null(SO[[1]]@assays$SCT)){
                SO <- PrepSCTIntegration(object.list = SO, anchor.features = integ_features)
                k.filter <- min(200, min(sapply(SO, ncol)))
                integ_anchors <- FindIntegrationAnchors(object.list = SO, normalization.method = "SCT", k.filter=k.filter, anchor.features = integ_features)
                SO_merge <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT",features.to.integrate = all_features)
                #SO_merge <- ScaleData(SO_merge,features=all_features)
            }
            else{
                k.filter <- min(200, min(sapply(SO, ncol)))
                integ_anchors <- FindIntegrationAnchors(object.list = SO, k.filter=k.filter, anchor.features = integ_features)
                SO_merge <- IntegrateData(anchorset = integ_anchors,features.to.integrate = all_features)
                #SO_merge <- ScaleData(SO_merge,features=all_features)  
            }}
    }

    SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = 2000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst", verbose = FALSE)
    SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,seed.use = 42)
    SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", dims = 1:npcs, seed.use=42)
    SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = 1)
    SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)
    

    #check for CITE-seq data and if so, run reductions
     if(doCiteSeq) {
        # Moved below integration step. SO_merge is recreated and this information was lost
        #SO_merge <- ScaleData(SO_merge, assay = "Protein")

        print("finding protein variable features...")
        VariableFeatures(SO_merge,assay="Protein") <- rownames(SO_merge$Protein)
        #Support for partial
        if(all(sapply(seq_along(SO),function(i) "Protein" %in% names(SO[[i]]@assays)))){
            print("running protein pca...")
            SO_merge <- RunPCA(object = SO_merge, assay="Protein",npcs = npcs,verbose = FALSE,reduction.name="protein_pca",seed.use = 42)
            SO_merge <- RunUMAP(object = SO_merge, assay="Protein", features=rownames(SO_merge$Protein), reduction.name="protein_umap",seed.use=42)
            SO_merge <- RunTSNE(object = SO_merge, assay="Protein", features=rownames(SO_merge$Protein),seed.use = 1,reduction.name="protein_tsne",check_duplicates=F)
            SO_merge <- FindNeighbors(SO_merge, assay="Protein",graph.name="Protein_snn",features=rownames(SO_merge$Protein))
        }else{
            doCiteSeq <- FALSE #set to false so we don't cluster protein
        }
        
    } else {
        doCiteSeq <- FALSE
    }

    for (i in seq(0.2,1.2,0.2)){
    SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 1)
        if(doCiteSeq){
            SO_merge <- FindClusters(SO_merge, graph.name="Protein_snn",resolution = i, algorithm = 1)
        }
    }
    print("Clustering successful!")
    
    n <- 60
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    qual_col_pals = qual_col_pals[c(7,6,2,1,8,3,4,5),]
    cols = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

    grobsList = list()
    if(TRUE){
    p1 <- DimPlot(SO_merge, reduction = "tsne", group.by = "orig.ident", repel = TRUE,          pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
    guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) +         ggtitle("RNA TSNE")
    grobsList[[length(grobsList)+1]] <- p1
    print("Added RNA TSNE")
    print(length(grobsList))

    }
    if(TRUE){
    p2 <- DimPlot(SO_merge, reduction = "umap", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
    guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("RNA UMAP")
    grobsList[[length(grobsList)+1]] <- p2
    print("Added RNA UMAP")
    print(length(grobsList))

    }
    if(TRUE & doCiteSeq){ 
    p3 <- DimPlot(SO_merge, reduction = "protein_tsne", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
    guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("Antibody TSNE")
    grobsList[[length(grobsList)+1]] <- p3
    print("Added Antibody TSNE")
    print(length(grobsList))

    }
    if(TRUE & doCiteSeq){ 
    p4 <- DimPlot(SO_merge, reduction = "protein_umap", group.by = "orig.ident", repel = TRUE, pt.size=0.02) + theme_classic() + scale_color_manual(values=cols) + theme(legend.position="top", legend.text=element_text(size=5)) +
    guides(colour = guide_legend(ncol=3, override.aes = list(size=1, alpha = 1))) + ggtitle("Antibody UMAP")
    grobsList[[length(grobsList)+1]] <- p4
    print("Added Antibody UMAP")
    print(length(grobsList))

    }

    
    n = ceiling(length(grobsList)^0.5)
    m = ceiling(length(grobsList)/n)
    imageWidth = 1200 * n
    imageHeight = 1200 * m
    dpi = 300

    grobs=arrangeGrob(grobs = grobsList, ncol = n)
    
    
    
    if (Image_type == 'png') {
    png(
      filename=paste0(image_out_dir, "/", "combined_renormalized_result.png"),
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=2,
      bg="white",
      res=dpi,
      type="cairo")
     }else {
        svglite::svglite(
        file=paste0(image_out_dir, "/", "combined_renormalized_result.svg"),
        width=round(imageWidth/dpi,digits=2),
        height=round(imageHeight/dpi,digits=2),
        pointsize=1,
        bg="white")
    }

    plot(grobs)
    
    dev.off()
    slot(SO_merge,"commands") <- list()
    cat("\nPCA Object Checksum:\n")
    #print(digest::digest(SO_merge))
    return(list(value=SO_merge))
}

