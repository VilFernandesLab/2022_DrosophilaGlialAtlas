# Open packages
library(Seurat)
library(ggplot2)

# Select non-glia clusters from embryo dataset
        # open RDS file
        data <- readRDS('embryo_neuronsglia_10percent_jun16.rds')
        data = UpdateSeuratObject(object = data)
        # remove glial clusters (cluster numbers stored in metadata as RNA_snn_res.10)
        data <- subset(data, invert = T, subset = RNA_snn_res.10 %in% c(57,98,108,89,104,106,150,119,84,72,29,79)) # using 'invert=T' will remove all clusters indicated in 'subset = RNA_snn_res.10 %in%'
        # run general Seurat pipeline
        data <- NormalizeData(data)
        data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
        data <- ScaleData(data, verbose = FALSE)
        data <- RunPCA(data, npcs = 30, verbose = FALSE)
        data <- FindNeighbors(data, dims = 1:16)
        data <- FindClusters(data, resolution = 0.5)
        data <- RunUMAP(data, reduction = "pca", dims = 1:20)
        plot <- DimPlot(data, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5)
        plot
        # plot FeaturePlot for UMAP with feature (gene) expression
        gene <- "Rdl"
        featureplot <- FeaturePlot(data, features = gene, label = F, min.cutoff = 0, cols = c('lightgrey','#211752')) 
        featureplot
        # to save FeaturePlot image as png - same code used to save any plot within this script
        ggsave('embryo_neurons_fne.png', featureplot, dpi = 600, device = "png")
        # define all clusters remaining as 1 in metadata 'ALL', to call in 'group.by =' ALL' argument, to give an average of all clusters
        data$all <- c(1)
        # plot RidgePlot (distribution of gene expression) with average of all clusters
        gene <- "Rdl"
        ridgeplot <- RidgePlot(data,features = gene, group.by = 'all', cols = 'lightgrey') + NoLegend() 
                        + scale_x_continuous(limits = c(-1, 5), breaks=c(0,1,2,3,4,5),expand = c(0.001, 0.001))
        ridgeplot
        # to save RDS file analysed
        saveRDS(data, file = "embryo_neuronsOnly_10percent_jun16.rds")

# Select glia clusters from embryo dataset
        # open RDS file
        data <- readRDS('embryo_neuronsglia_10percent_jun16.rds')
        data = UpdateSeuratObject(object = data)
        # select glial clusters (cluster numbers stored in metadata as RNA_snn_res.10) and save
        data <- subset(data, invert = F, subset = RNA_snn_res.10 %in% c(57,98,108,89,104,106,150,119,84,72,29,79))
        saveRDS(data, file = "embryo_glia_10mito_050522.rds")
        # plot UMAP
        plot <- DimPlot(data, reduction = "umap", label = TRUE,repel = TRUE,label.size = 5)
        plot
        # plot FeaturePlot for UMAP with feature (gene) expression
        gene <- "Rdl"
        featureplot <- FeaturePlot(data, features = gene, label = T, min.cutoff = 0,cols = c('lightgrey','#211752')) 
        featureplot

# Clean-up embryo data
        # RidgePlots to assess distribution of expression levels within clusters (used either Rdl, Frq1 or Nckx30C)
        gene <- "Rdl"
        ridgeplot <- RidgePlot(data,features = gene,group.by = 'seurat_clusters') + NoLegend() 
        ridgeplot
        # remove potential 'neuronal contamination' with high expression of validated neuronal genes
        dataNoRdl <- subset(data, subset = Rdl <= 1)
        dataNoRdlFrq1 <- subset(dataNoRdl, subset = Frq1 <= 1)
        dataNoRdlFrq1Nckx30C <- subset(dataNoRdlFrq1, subset = Nckx30C <= 1)
        saveRDS(dataNoRdlFrq1Nckx30C, file = "embryo_glia_10mito_050522_neuroclean1.rds")
        # remove potential hemocyte cells with expression of Hml (hemocyte marker)
        dataNoHemocytes <- subset(dataNoRdlFrq1Nckx30C, subset = Hml <= 0)
        saveRDS(dataNoHemocytes, file = "embryo_glia_10mito_050522_neuroclean1_HmlClean.rds")
        
# Re-analysis of cleaned data using general Seurat pipeline for clustering
        embryo <- readRDS("embryo_glia_10mito_050522_neuroclean1_HmlClean.rds")
        embryo <- NormalizeData(embryo)
        embryo <- FindVariableFeatures(embryo, selection.method = "vst", nfeatures = 2000)
        embryo <- ScaleData(embryo, verbose = FALSE)
        embryo <- RunPCA(embryo, npcs = 30, verbose = FALSE)
        embryo <- JackStraw(embryo, num.replicate = 100,dims = 20)
        embryo <- ScoreJackStraw(embryo, dims = 1:20)
        JackStrawPlot(embryo, dims = 1:15)
        ElbowPlot(embryo)
        embryo <- FindNeighbors(embryo, dims = 1:16)
        embryo <- FindClusters(embryo, resolution = 0.5)
        embryo <- RunUMAP(embryo, reduction = "pca", dims = 1:20)
        plot1 <- DimPlot(embryo, reduction = "umap", label = TRUE,repel = TRUE,label.size = 6)
        plot1
        saveRDS(embryo, file = "embryo_glia_10mito_050522_neuroclean1_HmlClean_reanalysed16dims05res.rds")

# Get csv file with gene markers for each cluster compared to all other clusters
        data <- readRDS("embryo_glia_10mito_050522_neuroclean1_HmlClean_reanalysed16dims05res.rds")
        allmarkers <- FindAllMarkers(data, only.pos = TRUE)
        write.csv(allmarkers, 'AllMarkers_embryo_glia_10mito_050522_neuroclean1_HmlClean_reanalysed16dims05res.csv')

# Plot FeaturePlot for UMAP with feature (gene) expression (clean previous objects)
        data <- readRDS("embryo_glia_10mito_050522_neuroclean1_HmlClean_reanalysed16dims05res.rds")
        gene <- "wrapper"
        plot.gene <- FeaturePlot(data, features = gene, label = TRUE, min.cutoff = 0)
        plot.gene

# Give glial subtype names to clusters (clean previous objects)
        data <- readRDS("embryo_glia_10mito_050522_neuroclean1_HmlClean_reanalysed16dims05res.rds")
        data <- RenameIdents(object = data, `8` = "PNSg_3",`2` = "Unknown_1", 
                             `0` = "PNSg_1", `3` = "Astrocyte",`5` = "Unknown_2",
                             `7` = "Ensheathing",`9` = "Subperineurial",`10` = "MAPg",
                             `6` = "PNSg_2",`4` = "Perineurial",`1` = "Cortex")
        plot <- DimPlot(data, reduction = "umap", label = TRUE,label.size = 10,repel = TRUE) + guides(color = guide_legend(override.aes = list(size=4), ncol=1)) + theme(legend.text=element_text(size=20))# used either group.by = 'FinalIdents' (Desplan) or 'subtype' (Zipursky) or 'AllClusterNrs' (just numbers)
        plot
        # save all glial subtype names in metadata as AllClusterNamesEmbryo
        data[["AllClusterNamesEmbryo"]] <- Idents(object = data)
        saveRDS(data, file = "embryoS17glia.rds")
        
# Expression plot (DotPlot) for marker genes within each cluster (clean previous objects)
        # open packages
        library(tidyr)
        library(dplyr)
        library(magrittr)
        library(ggtext)
        # open file
        data <- readRDS('embryoS17glia.rds')
        # indication to use gene expression on the RNA assay
        DefaultAssay(data) <- "RNA"
        # plot DotPlot with 2 y axis (left axis with cluster numbers and right axis with corresponding glial subtype names)
                # run the DotPlot for marker genes
                plot <- DotPlot(data, scale= T, features = c('CG4797','CG6126','CG5080','PRL-1',
                                            'ltl','moody','wrapper','ana',
                                            'CG9449','Eaat2','CG9657','pum','Fur1',
                                            'Tet','Tre1','alrm','Oaz')) # rna_
                # retrieve data from DotPlot
                plotData <- plot$data
                # create two y axes
                        # define the clusters for the right y axis
                        plotData$id <- factor(plotData$id, levels = rev(c("Perineurial", "MAPg","Subperineurial","Cortex",
                                                  "Ensheathing","Astrocyte","Unknown_1", "Unknown_2",
                                                  "PNSg_1","PNSg_2","PNSg_3")))
                ylabs2 <- levels(plotData$id)
                # define the clusters for the left y axis
                ylabs1 <- as.character(c(8,6,0,5,2,3,7,1,9,10,4))
        # plot the data from DotPlot in graph with 2 y axis
        final.plot <- plotData %>% filter(avg.exp.scaled > 0) %>% ggplot(aes(x = features.plot, y = as.numeric(id), color = avg.exp.scaled, size = pct.exp)) +
                        geom_point() +
                        scale_y_continuous( breaks = 1:length(ylabs1),
                                            labels = ylabs1,
                                            sec.axis = sec_axis(trans = ~.,
                                                                breaks = 1:length(ylabs2),
                                                                labels = ylabs2, name="Annotated Clusters",)) +
                        scale_x_discrete() + theme_classic(14) +
                        theme(axis.text.x = element_text(face = "italic", color="black",angle = 45, hjust = 1,size = 12),
                              axis.text.y = element_text(size = 14,color = "black")) +
                        labs(x= "Marker genes", y = "Clusters", color = "Avg. Exp", size = "% Exp") +
                        scale_color_gradientn(colours=c("white", "#153743")) +
                        scale_size_continuous(range=c(0, 6))
        final.plot
