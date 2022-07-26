# Open packages
library(Seurat)
library(patchwork)
library(ggplot2)

# Open original RDS file and select the glial clusters from the Ozel Adult dataset (GSE142787_Adult)
        data.O <- readRDS('GSE142787_Adult.rds')
        data.O <- subset(data.O, subset = FinalIdents %in% c('189','185','186','187','188','190','196','197','198','199','200','201','202','203','204','205','206','207','208'))
        saveRDS(data.O, file = 'Glia_Adult.rds') # to save RDS file analysed

# Open original 10X data files and select the glial clusters and 96h APF time point from the Kurmangaliyev dataset (MainData)
        # read 10x files (should have barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz)
        tenXdata.K <- Read10X(data.dir = "~/Glia scRNA-seq_Zipursky_2020/Main_dataset/MainData")
        # transform metadata file to dataframe
        metadata <- read.table(file = "~/Glia scRNA-seq_Zipursky_2020/Main_dataset/MainData/metadata.tsv", row.names = 1, header = TRUE)
        metadata.dataframe <- as.data.frame(metadata)
        # create RDS file from 10x files and metadata
        seurat_object_metadata <- CreateSeuratObject(counts = tenXdata.K, meta.data = metadata.dataframe)
        # select glial clusters (G as the first letter of glial clusters, according to Zipursky study) and 96h timepoint
        data.O.96h <- subset(zipursky, subset = time == '96h')
        data.O.96h.glia <- subset(data.O.96h, subset = class == 'G')
        saveRDS(data.O.96h.glia, file = '96h_Glia_metadata.rds')

# Ozel glia data
        # open file and run general Seurat pipeline to run UMAP
        data.O <- readRDS('Glia_Adult.rds')
        data.O <- FindVariableFeatures(data.O, selection.method = "vst", nfeatures = 2000)
        data.O <- ScaleData(data.O, verbose = FALSE)
        data.O <- RunPCA(data.O, npcs = 30, verbose = FALSE)
        data.O <- RunUMAP(data.O, reduction = "pca", dims = 1:20)
        # plot UMAP
        plot.O <- DimPlot(data.O, reduction = "umap", label = TRUE,repel = TRUE,label.size = 5,group.by = 'FinalIdents') # use 'FinalIdents' to show clusters obtained on the original dataset
        plot.O
        # to save DimPlot image as png - same code used to save any plot within this script
        ggsave('UMAP_Ozel.png', plot.O, dpi = 600, device = "png")
        # to save RDS file analysed
        saveRDS(data.O, file = "Glia_Adult_analysed.rds") 

# Kurmangaliyev glia data
        # open file and run general Seurat pipeline to run UMAP
        data.K <- readRDS('96h_Glia_metadata.rds')
        data.K <- NormalizeData(data.K)
        data.K <- FindVariableFeatures(data.K, selection.method = "vst", nfeatures = 2000)
        data.K <- ScaleData(data.K, verbose = FALSE)
        data.K <- RunPCA(data.K, npcs = 30, verbose = FALSE)
        data.K <- RunUMAP(data.K, reduction = "pca", dims = 1:20)
        # plot UMAP
        plot.K <- DimPlot(data.K, reduction = "umap", label = TRUE,repel = TRUE,label.size = 5,group.by = 'subtype') # use 'subtype' to show clusters obtained on the original dataset
        plot.K
        saveRDS(data, file = "96h_Glia_metadata_analysed.rds")

# Integration of Ozel and Kurmangaliyev datasets (clear all previous objects)
        # open and calculate most variable features for Ozel glia dataset
        data.O <- readRDS('Glia_Adult.rds')
        data.O <- FindVariableFeatures(data.O, selection.method = "vst", nfeatures = 2000)
        # open and calculate most variable features for Kurmangaliyev glia dataset
        data.K <- readRDS('96h_Glia_metadata.rds')
        data.K <- NormalizeData(data.K)
        data.K <- FindVariableFeatures(data.K, selection.method = "vst", nfeatures = 2000)
        # integrate the two datasets
        reference.list <- c(data.O,data.K)
        anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        # general Seurat pipeline for clustering
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- JackStraw(integrated, num.replicate = 100,dims = 20)
        integrated <- ScoreJackStraw(integrated, dims = 1:20)
        JackStrawPlot(integrated, dims = 1:20)
        ElbowPlot(integrated)
        integrated <- FindNeighbors(integrated, dims = 1:16) # dimensions used here informed by the ElbowPlot
        integrated <- FindClusters(integrated, resolution = 0.5)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
        # plot UMAP
        plot <- DimPlot(integrated, reduction = "umap", label = TRUE,repel = TRUE,group.by = 'FinalIdents') # used either group.by = 'FinalIdents' (Ozel) or 'subtype' (Kurmangaliyev)
        plot
        saveRDS(integrated, file = "Integrated_Glia_Adult_96h_Glia_metadata_16dim,05res.rds")

        # plot FeaturePlot for UMAP with feature (gene) expression
        gene <- "Rdl" # used either Rdl, Frq1 or Nckx30C
        plot.gene <- FeaturePlot(integrated, features = gene, label = TRUE, min.cutoff = 0) 
        plot.gene
        
# Clean-up integrated data (clear all previous objects)
        # open RDS file of integrated data
        integrated <- readRDS('Integrated_Glia_Adult_96h_Glia_metadata_16dim,05res.rds')
        # indication to use gene expression on the RNA assay
        DefaultAssay(integrated) <- "RNA"
        # RidgePlots to assess distribution of expression levels within clusters (used either Rdl, Frq1 or Nckx30C)
        ridgeplot <- RidgePlot(integrated,features = "Rdl",group.by = 'seurat_clusters') + NoLegend() 
        ridgeplot
        # remove putative neuronal contamination with high expression of validated neuronal genes
        integratedNoRdl <- subset(integrated, subset = Rdl <= 1)
        integratedNoRdlFrq1 <- subset(integratedNoRdl, subset = Frq1 <= 1)
        integratedNoRdlFrq1Nckx30C <- subset(integratedNoRdlFrq1, subset = Nckx30C <= 1)
        saveRDS(integratedNoRdlFrq1Nckx30C, file = "Integrated_Glia_Adult_96h_Glia_metadata_16dim,05res_neuroclean1.rds")
        # remove putative hemocytes with expression of Hml (hemocyte marker)
        integratedNoHemocytes <- subset(integratedNoRdlFrq1Nckx30C, subset = Hml <= 0)
        saveRDS(integratedNoHemocytes, file = "Integrated_Glia_Adult_96h_Glia_metadata_16dim,05res_neuroclean1_HmlClean.rds")
        # plot UMAP after clean-up
        plot <- DimPlot(integratedNoHemocytes, reduction = "umap", label = TRUE,repel = TRUE) # used either group.by = 'FinalIdents' (Ozel) or 'subtype' (Kurmangaliyev)
        plot

# First re-integration (clean previous objects)
        # open files
        data.O <- readRDS('Glia_Adult.rds')
        data.K <- readRDS('96h_Glia_metadata.rds')
        data.integrated <- readRDS('Integrated_Glia_Adult_96h_Glia_metadata_16dim,05res_neuroclean1_HmlClean.rds')
        # retrieve cells left from original datasets
        all.cells <- WhichCells(data.integrated)
        data.O <- subset(data.O, cells = all.cells)
        data.K <- subset(data.K, cells = all.cells)
        # integrate
        data.O <- FindVariableFeatures(data.O, selection.method = "vst", nfeatures = 2000)
        data.K <- NormalizeData(data.K)
        data.K <- FindVariableFeatures(data.K, selection.method = "vst", nfeatures = 2000)
        reference.list <- c(data.O,data.K)
        anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- JackStraw(integrated, num.replicate = 100,dims = 20)
        integrated <- ScoreJackStraw(integrated, dims = 1:20)
        JackStrawPlot(integrated, dims = 1:20)
        ElbowPlot(integrated)
        integrated <- FindNeighbors(integrated, dims = 1:19)
        integrated <- FindClusters(integrated, resolution = 0.5)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:25)
        plot <- DimPlot(integrated, reduction = "umap", label = TRUE,repel = TRUE) # used either group.by = 'FinalIdents' (Ozel) or 'subtype' (Kurmangaliyev)
        plot
        saveRDS(integrated1, file = "Reintegrated_neuroclean1_HmlClean_19dim,05res.rds")

# Remove non-overlapping clusters (<1%)
        cluster <- WhichCells(integrated, idents = "18") # run this for every cluster number
        cluster.data.O <- subset(data.O, cells = cluster)
        cluster.data.K <- subset(data.K, cells = cluster)
        # check how many cells in cluster.data.O and cluster.data.K - if less than 1% of K in total remove cluster
        # clusters 3 and 6 and 8 were removed here
        IntegratedOverlap <- subset(integrated, subset = seurat_clusters %in% c(0,1,2,4,5,7,9,10,11,12,13,14,15,16,17,18))
        saveRDS(IntegratedOverlap, file = "Reintegrated_neuroclean1_HmlClean_19dim,05res_onlyOverlapClusters.rds")
        
# Second re-integration (clean previous objects)
        # open files
        data.O <- readRDS('Glia_Adult.rds')
        data.K <- readRDS('96h_Glia_metadata.rds')
        data.integrated <- readRDS('Reintegrated_neuroclean1_HmlClean_19dim,05res_onlyOverlapClusters.rds')
        # retrieve cells left from original datasets
        all.cells <- WhichCells(data.integrated)
        data.O <- subset(data.O, cells = all.cells)
        data.K <- subset(data.K, cells = all.cells)
        # integrate
        data.O <- FindVariableFeatures(data.O, selection.method = "vst", nfeatures = 2000)
        data.K <- NormalizeData(data.K)
        data.K <- FindVariableFeatures(data.K, selection.method = "vst", nfeatures = 2000)
        reference.list <- c(data.O,data.K)
        anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- JackStraw(integrated, num.replicate = 100,dims = 20)
        integrated <- ScoreJackStraw(integrated, dims = 1:20)
        JackStrawPlot(integrated, dims = 1:20)
        ElbowPlot(integrated)
        integrated <- FindNeighbors(integrated, dims = 1:19)
        integrated <- FindClusters(integrated, resolution = 0.5)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:25)
        plot <- DimPlot(integrated, reduction = "umap", label = TRUE,repel = TRUE,group.by = 'FinalIdents')
        plot
        saveRDS(integrated, file = "Reintegrated_neuroclean1_HmlClean_19dim,05res_onlyOverlapClusters_SecondReintegrated_19dim,05res.rds")
        
# Remove non-overlapping clusters (<1%) (clean previous objects) - After the previous reintegration, one cluster was now only composed of Ozel dataset cells
        data.integrated <- readRDS('Reintegrated_neuroclean1_HmlClean_19dim,05res_onlyOverlapClusters_SecondReintegrated_19dim,05res.rds')
        cluster <- WhichCells(data.integrated, idents = "18")
        cluster.data.O <- subset(data.O, cells = cluster)
        cluster.data.K <- subset(data.K, cells = cluster)
        # check how many cells in cluster.data.O and cluster.data.K - if less than 1% of K in total remove cluster
        # cluster 4 was removed here
        IntegratedOverlap <- subset(data.integrated, subset = seurat_clusters %in% c(0,1,2,3,5,6,7,8,9,10,11,12,13,14,15))
        saveRDS(IntegratedOverlap, file = "Reintegrated_neuroclean1_HmlClean_19dim,05res_onlyOverlapClusters_SecondReintegrated_19dim,05res_NoCluster4.rds")

# Third re-integration (clean previous objects)
        # open files
        data.O <- readRDS('Glia_Adult.rds')
        data.K <- readRDS('96h_Glia_metadata.rds')
        data.integrated <- readRDS('Reintegrated_neuroclean1_HmlClean_19dim,05res_onlyOverlapClusters_SecondReintegrated_19dim,05res_NoCluster4.rds')
        # retrieve cells left from original datasets
        all.cells <- WhichCells(data.integrated)
        data.O <- subset(data.O, cells = all.cells)
        data.K <- subset(data.K, cells = all.cells)
        # integrate
        data.O <- FindVariableFeatures(data.O, selection.method = "vst", nfeatures = 2000)
        data.K <- NormalizeData(data.K)
        data.K <- FindVariableFeatures(data.K, selection.method = "vst", nfeatures = 2000)
        reference.list <- c(data.O,data.K)
        anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- JackStraw(integrated, num.replicate = 100,dims = 20)
        integrated <- ScoreJackStraw(integrated, dims = 1:20)
        JackStrawPlot(integrated, dims = 1:20)
        ElbowPlot(integrated)
        integrated <- FindNeighbors(integrated, dims = 1:18)
        integrated <- FindClusters(integrated, resolution = 0.4)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:25)
        plot <- DimPlot(integrated, reduction = "umap", label = TRUE,repel = TRUE)
        plot
        saveRDS(integrated, file = "ThirdReintegrated_18dim,04res.rds")
        
# Subcluster clusters (clean previous objects)
        # open files
        data.O <- readRDS('Glia_Adult.rds')
        data.K <- readRDS('96h_Glia_metadata.rds')
        data.integrated <- readRDS('ThirdReintegrated_18dim,04res.rds') # using this file to subcluster cluster 8 and 9, that will be separated in the main UMAP
        data.integrated <- readRDS('ThirdReintegrated_18dim,04res_AllClusterNrs.rds') # using this file to subset all other clusters, including the 4 separated clusters (8,9,15,16)
        # retrieve cells left from original datasets
        cluster <- WhichCells(data.integrated, idents = "8")
        cluster.data.O <- subset(data.O, cells = cluster)
        cluster.data.K <- subset(data.K, cells = cluster)
        # integrate
        cluster.data.O <- FindVariableFeatures(cluster.data.O, selection.method = "vst", nfeatures = 2000)
        cluster.data.K <- NormalizeData(cluster.data.K)
        cluster.data.K <- FindVariableFeatures(cluster.data.K, selection.method = "vst", nfeatures = 2000)
        reference.list <- c(cluster.data.O,cluster.data.K)
                # value used is 1 unit less than the minimum number of cells in cluster.data.K or cluster.data.O, for a maximum of 30
                anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:27,k.score = 27,k.filter =27)
                integrated <- IntegrateData(anchorset = anchors, dims = 1:27,k.weight=27)
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- JackStraw(integrated, num.replicate = 100,dims = 20)
        integrated <- ScoreJackStraw(integrated, dims = 1:20)
        JackStrawPlot(integrated, dims = 1:20)
        ElbowPlot(integrated)
        integrated <- FindNeighbors(integrated, dims = 1:3)
        integrated <- FindClusters(integrated, resolution = 0.2)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:15)
        plot <- DimPlot(integrated, reduction = "umap", label = TRUE,repel = TRUE,label.size = 7) + NoLegend()
        plot
        saveRDS(integrated, file = "ThirdReintegrated_Cluster8_3dim02res.rds")
        # get csv file with gene markers for each cluster comparison
                # number of the clusters to compare; use c("2","1") if comparing a group of clusters
                clusterA <- "1" 
                clusterB <- "0" 
                markers <- FindMarkers(integrated, ident.1 = clusterA, ident.2 = clusterB, min.pct = 0.25) 
                write.csv(markers, '1vs0_ThirdReintegrated_Cluster8_3dim02res.csv')
        # plot FeaturePlot for UMAP with feature (gene) expression
        gene <- "Obp18a"
        plot.gene <- FeaturePlot(integrated, features = gene, label = TRUE, min.cutoff = 0)
        plot.gene
        
# Add identity of subclusters in main dataset (clean previous objects)
        # open files
        data.integrated <- readRDS('ThirdReintegrated_18dim,04res.rds')
        subperineurial <- readRDS('ThirdReintegrated_Cluster8_3dim02res.rds')
        marginal <- readRDS('ThirdReintegrated_Cluster9_5dim05res.rds')
        # retrieve cluster 1 cell identities to mark them as a new cluster in the main UMAP   
        subperineurial <- WhichCells(subperineurial, idents = "1")
        marginal <- WhichCells(marginal, idents = "1")
        # add new identities to the main dataset (active.ident)
        Idents(data.integrated, cells = subperineurial) <- "15"
        Idents(data.integrated, cells = marginal) <- "16"
        # save all 17 cluster numbers in metadata as AllClusterNrs
        data.integrated[["AllClusterNrs"]] <- Idents(object = data.integrated)
        saveRDS(data.integrated, file = "ThirdReintegrated_18dim,04res_AllClusterNrs.rds")
        
# Select markers for validation (clean previous objects)
        # open RDS
        integrated.adult <- readRDS("ThirdReintegrated_18dim,04res_AllClusterNrs.rds")
        # get csv file with gene markers for each cluster comparison
        clusterA <- "1" # number of the clusters to compare; use c("2","1") if comparing a group of clusters
        clusterB <- "0" 
        markers <- FindMarkers(integrated.adult, ident.1 = clusterA, ident.2 = clusterB, min.pct = 0.25)
        write.csv(markers, 'Markers_Cluster0_ThirdReintegrated_18dim,04res_AllClusterNrs.csv')
        # plot FeaturePlot for UMAP with feature (gene) expression
        gene <- "Tet"
        plot.gene <- FeaturePlot(integrated.adult, features = gene, label = T, min.cutoff = 0,cols = c('lightgrey','#211752'),label.size = 6)
        plot.gene
        
# Give glial subtype names to clusters (clean previous objects)
        integrated.adult <- readRDS('ThirdReintegrated_18dim,04res_AllClusterNrs.rds')
        integrated.adult <- RenameIdents(object = integrated.adult, `8` = "Fenestrated",`2` = "Astrocyte_2", 
                     `0` = "Astrocyte_1", `3` = "Ensheathing_3",`5` = "Ensheathing_2",
                     `7` = "Ensheathing_1",`9` = "Epithelial",`10` = "Chiasm",
                     `6` = "Distal satellite",`14` = "Proximal satellite",
                     `11` = "Cortex_2",`13` = "Cortex_1",`15` = "Subperineurial",
                     `4` = "Chalice",`1` = "Perineurial",`12` = 'Pseudo-cartridge',`16` = "Marginal")
        plot <- DimPlot(integrated.adult, reduction = "umap", label = TRUE,label.size = 10,repel = TRUE) + guides(color = guide_legend(override.aes = list(size=4), ncol=1)) + theme(legend.text=element_text(size=20))# used either group.by = 'FinalIdents' (Desplan) or 'subtype' (Zipursky) or 'AllClusterNrs' (just numbers)
        plot
        # save all glial subtype names in metadata as AllClusterNames
        integrated.adult[["AllClusterNames"]] <- Idents(object = integrated.adult)
        saveRDS(integrated.adult, file = "adultOLglia.rds")
        
# Expression plot (DotPlot) for marker genes within each cluster (clean previous objects)
        # open packages
        library(tidyr)
        library(dplyr)
        library(ggtext)
        library(magrittr)
        # open file
        data <- readRDS('adultOLglia.rds')
        # indication to use gene expression on the RNA assay
        DefaultAssay(data) <- "RNA"
        # plot DotPlot with 2 y axis (left axis with cluster numbers and right axis with corresponding glial subtype names)
                # run the DotPlot for marker genes
                plot <- DotPlot(data, scale= T,features = c('CG6126','CG9743','Tret1-1',
                                                        'JhI-21','Cyp311a1','CG7135','Fas3','CG13003',
                                                        'GstE9','ltl','NimB4','DAT','wrapper','Rfabg',
                                                        'Dop1R1','ana','Ork1','Obp18a','Dtg','Ndae1',
                                                        'CG43795','axo','List','GstT4',
                                                        'htl','mbc','Tet','wun2','Act79B'))
                # retrieve data from DotPlot
                plotData <- plot$data
                # create two y axes
                        # define the clusters for the right y axis
                        plotData$id <- factor(plotData$id, levels = rev(c("Perineurial", "Chalice", "Fenestrated", "Subperineurial",
                                                                          "Pseudo-cartridge", "Chiasm", "Cortex_1", "Cortex_2",
                                                                          "Distal satellite", "Proximal satellite", "Epithelial",
                                                                          "Ensheathing_1", "Ensheathing_2","Ensheathing_3",
                                                                          "Marginal","Astrocyte_1","Astrocyte_2")))
                        ylabs2 <- levels(plotData$id)
                        # define the clusters for the left y axis
                        ylabs1 <- as.character(c(16,3,5,7,9,2,0,14,6,13,11,10,12,15,8,4,1))
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
                                      axis.text.y = element_text(size = 14, color = "black")) +
                                labs(x= "Marker genes", y = "Clusters", color = "Avg. Exp", size = "% Exp") +
                                scale_color_gradientn(colours=c("white", "#153743")) +
                                scale_size_continuous(range=c(2, 6))
                final.plot

# Plot and save image of gene expression UMAPs for naz+CG31235 (same gene, different names) (clean previous objects)
data <- readRDS('ThirdReintegrated_18dim,04res_AllClusterNrs.rds')
plot <- FeaturePlot(data, min.cutoff = 0,max.cutoff = 2, pt.size = 0.5, features = c('naz',"CG31235"), cols = c('lightgrey','#211752','#211752'),blend = TRUE)
plot[[3]]
ggsave('featureplot_naz.png', plot[[3]], dpi = 600, device = "png")

# Plot and save RidgePlot for distribution of gene expression (clean previous objects)
data <- readRDS('Integrated_Glia_Adult_96h_Glia_metadata_16dim,05res.rds')
gene <- "Frq1"
plot <- RidgePlot(data, features = gene, group.by = 'seurat_clusters') + NoLegend() 
        + scale_x_continuous(limits = c(-1, 5), breaks=c(0,1,2,3,4,5),expand = c(0.001, 0.001))
plot
ggsave('RidgePlot_Adult_Nckx30C.png', plot, dpi = 600, device = "png")

