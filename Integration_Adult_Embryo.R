# Open packages
library(Seurat)
library(ggplot2)

# We used the RDS object of late embryo VNC and young adult optic lobe before renaming the clusters so we can rename them now and distinguish the clusters in the integration as being embryo or adult
# Add cluster names with 'a_' to identify adult clusters in integration
dataA <- readRDS('ThirdReintegrated_18dim,04res_AllClusterNrs.rds')
dataA <- RenameIdents(object = dataA, `8` = "a_Fenestrated",`2` = "a_Astrocyte_2", 
                     `0` = "a_Astrocyte_1", `3` = "a_Ensheathing_3",`5` = "a_Ensheathing_2",
                     `7` = "a_Ensheathing_1",`9` = "a_Epithelial",`10` = "a_Chiasm",
                     `6` = "a_Distal satellite",`14` = "a_Proximal satellite",
                     `11` = "a_Cortex_2",`13` = "a_Cortex_1",`15` = "a_Subperineurial",
                     `4` = "a_Chalice",`1` = "a_Perineurial",`12` = 'a_Pseudo-cartridge',`16` = "a_Marginal")
dataA[["NamesForIntegration"]] <- Idents(object = dataA)
saveRDS(dataA, file = "ThirdReintegrated_18dim,04res_AllClusterNrs_NamesForIntegration.rds")

# Add cluster names with 'e_' to identify embryo clusters in integration
dataE <- readRDS("embryo_glia_10mito_050522_neuroclean1_HmlClean_reanalysed16dims05res.rds")
dataE <- RenameIdents(object = dataE, `8` = "e_PNSg_3",`2` = "e_Unknown_1", 
                     `0` = "e_PNSg_1", `3` = "e_Astrocyte",`5` = "e_Unknown_2",
                     `7` = "e_Ensheathing",`9` = "e_Subperineurial",`10` = "e_Channel perineurial",
                     `6` = "e_PNSg_2",`4` = "e_Surface-only perineurial",`1` = "e_Cortex")
dataE[["NamesForIntegrationEmbryo"]] <- Idents(object = dataE)
saveRDS(dataE, file = "embryo_glia_10mito_050522_neuroclean1_HmlClean_reanalysed16dims05res_NamesForIntegration2.rds")

# Integration (clean previous objects)
        # open and calculate most variable features for young adult glial dataset
        adult <- readRDS('ThirdReintegrated_18dim,04res_AllClusterNrs_NamesForIntegration.rds')
        adult <- FindVariableFeatures(adult, selection.method = "vst", nfeatures = 2000)
        # open and calculate most variable features for embryonic glial dataset
        embryo <- readRDS('embryo_glia_10mito_050522_neuroclean1_HmlClean_reanalysed16dims05res_NamesForIntegration2.rds')
        embryo <- NormalizeData(embryo)
        embryo <- FindVariableFeatures(embryo, selection.method = "vst", nfeatures = 2000)
        # integrate the two datasets
        reference.list <- c(adult,embryo)
        anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
        integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
        # general Seurat pipeline for clustering
        integrated <- ScaleData(integrated, verbose = FALSE)
        integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
        integrated <- JackStraw(integrated, num.replicate = 100,dims = 20)
        integrated <- ScoreJackStraw(integrated, dims = 1:20)
        JackStrawPlot(integrated, dims = 1:20)
        ElbowPlot(integrated)
        integrated <- FindNeighbors(integrated, dims = 1:16)
        integrated <- FindClusters(integrated, resolution = 0.5)
        integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
        saveRDS(integrated1, file = "Integrated_Adult_Embryo_16dim,05res.rds")

# Plot with coloured clusters as in Fig.9B - Embryo clusters (clean previous objects)
        # open RDS
        data <- readRDS('Integrated_Adult_Embryo_16dim,05res.rds')
        # define cluster colours
        color_list_new2 <- c('#0f87a4','#015c27','#8e5198','#f7a26a','#f08d41','#dde336','#bec257','#a19f3b','#e2391e','#454d47','#6e7570')
        plot <- DimPlot(data, label.size = 5,reduction = "umap", label = F,repel = TRUE,group.by = 'NamesForIntegrationEmbryo',cols= color_list_new2,na.value='grey80') + theme(legend.text = element_text(size = 14))
        plot

# Plot with coloured clusters as in Fig.9C - Adult clusters (clean previous objects)
        # open RDS
        data <- readRDS('Integrated_Adult_Embryo_16dim,05res.rds')
        # define cluster colours
        color_list_new <- c('#0f87a4','#64b2c4','#f7a26a','#f3d033','#015c27','#048a3c','#1ad668','#8e5198','#ccb2d3','#ac80b4','#885690','#d96b16','#6be0fa','#f08d41','#26bd65','#e66c5e','#e2391e')
        plot <- DimPlot(data, label.size = 5,reduction = "umap", label = F,repel = TRUE,group.by = 'NamesForIntegration',cols= color_list_new,na.value='grey80') + theme(legend.text = element_text(size = 14))
        plot

# (clean previous objects) Plot UMAP as in Fig.9A, with embryo clusters in #B52828 colour and adult clusters in #1D54A1 colour
# Re-do embryo data with NamesForIntegration for cluster names instead of NamesForIntegrationEmbryo (not saved)
        data <- readRDS('Integrated_Adult_Embryo_16dim,05res.rds')
        color_list_new2 <- c('#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#1D54A1','#B52828','#B52828','#B52828','#B52828','#B52828','#B52828','#B52828','#B52828','#B52828','#B52828','#B52828')
        plot <- DimPlot(data, label.size = 5,reduction = "umap", label = F,repel = TRUE,group.by = 'NamesForIntegration',cols= color_list_new2,na.value='grey80') + NoLegend()
        plot

# Expression plot (DotPlot) for marker genes within each cluster (clean previous objects)
        # open packages
        library(tidyr)
        library(dplyr)
        library(ggtext)
        library(magrittr)
        # open file
        integrated <- readRDS('Integrated_Adult_Embryo_16dim,05res.rds')
        # subset clusters used in Fig9C
        data.subset.adult <- subset(integrated, subset = NamesForIntegration == c("a_Perineurial","a_Chalice","a_Subperineurial",
                                                                "a_Ensheathing_1", "a_Ensheathing_2","a_Ensheathing_3","a_Cortex_1", 
                                                                 "a_Cortex_2","a_Astrocyte_1","a_Astrocyte_2"))
        data.subset.embryo <- subset(integrated, subset = NamesForIntegrationEmbryo == c("e_Surface-only perineurial","e_Channel perineurial","e_Subperineurial","e_Ensheathing", "e_Cortex","e_Astrocyte"))
        # indication to use gene expression on the RNA assay
        DefaultAssay(data.subset.adult) <- "RNA"
        DefaultAssay(data.subset.embryo) <- "RNA"
        # run the DotPlot for marker genes
        plot <- DotPlot(data.subset.adult,features = c('CG6126','ltl','ana','wrapper','Tet'),group.by = 'NamesForIntegration')
        plot1 <- DotPlot(data.subset.embryo,features = c('CG6126','ltl','ana','wrapper','Tet'),group.by = 'NamesForIntegrationEmbryo')
        plot
        plot1
        # retrieve data from DotPlot
        plotData <- plot$data
        plotData1 <- plot1$data
        # combine the two datasets
        plotData.Combined <- bind_rows(plotData, plotData1)
        # define the clusters for y axis
        plotData.Combined$id <- factor(plotData.Combined$id, levels = rev(c("e_Surface-only perineurial","a_Perineurial", "e_Channel perineurial","a_Chalice",
                                                          "e_Subperineurial","a_Subperineurial","e_Cortex", 
                                                          "a_Cortex_1", "a_Cortex_2","e_Ensheathing",
                                                          "a_Ensheathing_1","a_Ensheathing_2","a_Ensheathing_3",
                                                          "e_Astrocyte","a_Astrocyte_1","a_Astrocyte_2")))
        # plot the data in combined DotPlot
        final.plot <- plotData.Combined %>% filter(avg.exp > 0) %>% ggplot(aes(x = features.plot, y = id, color = avg.exp, size = pct.exp)) +
                                geom_point() +
                                scale_x_discrete() + theme_classic(16) +
                                theme(axis.text.x = element_text(face = "italic", color="black",angle = 45, hjust = 1,size = 14),
                                      axis.text.y = element_text(color="black",size = 14),
                                      axis.title.y = element_text(margin = margin(t = 0, r = 40, b = 0, l = 0))) +
                                labs(x= "Marker genes", y = "Annotated Clusters", color = "Avg. Exp", size = "% Exp") +
                                scale_color_gradientn(colours=c("white", "#153743")) +
                                scale_size_continuous(range=c(0, 6))
        final.plot
        ggsave('DotPlot_adult_embryo.png', final.plot, dpi = 600, device = "png") # save plot as png image
