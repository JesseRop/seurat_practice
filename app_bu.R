library(dplyr)
library(Seurat)
library(d3heatmap)
library(shinythemes)
library(shiny)

# Load the PBMC dataset

#pbmc.data = Read10X(data.dir = "D:/GCRF_UoG/tuitorials/seurat/filtered_gene_bc_matrices/hg19/")
pbmc.data = Read10X(data.dir = "~/scrna_shiny/git_seurat/seurat_practice_files/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc = CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# Lets examine a few genes in the first thirty cells
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
dense.size = object.size(as.matrix(pbmc.data))

sparse.size = object.size(pbmc.data)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")

head(pbmc[["percent.mt"]])

# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 = FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

pbmc = subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc = NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc = FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 = head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 = VariableFeaturePlot(pbmc)
plot2 = LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge =0, ynudge=0)
CombinePlots(plots = list(plot1, plot2))

all.genes = rownames(pbmc)
pbmc = ScaleData(pbmc, features = all.genes)

pbmc_pca = RunPCA(pbmc, features = VariableFeatures(object = pbmc))

print(pbmc_pca[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc_pca, dims = 1:2, reduction = "pca")

DimPlot(pbmc_pca, reduction = "pca")

pbmc = JackStraw(pbmc_pca, num.replicate = 100)
pbmc = ScoreJackStraw(pbmc, dims = 1:20)
#JackStrawPlot(pbmc, dims = 1:15)
#ElbowPlot(pbmc)

pbmc = FindNeighbors(pbmc, dims = 1:10)
pbmc = FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

dim1 = 5
dim2 = 10

top_markers_start = 2 
top_markers_end = 10

##tSNE and clustering
#pbmc.umap.list = lapply(dim1:dim2, function(x) {
# RunUMAP(pbmc, reduction = "pca", dims = 1:x)
#}) 
#saveRDS(pbmc.umap.list, "~/scrna_shiny/git_seurat/seurat_practice_files/pbmc.umap.list")
pbmc.umap.list = readRDS("~/scrna_shiny/git_seurat/seurat_practice_files/pbmc.umap.list")

#top_markers = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#saveRDS(top_markers, "~/scrna_shiny/git_seurat/seurat_practice_files/top_markers")
top_markers = readRDS("~/scrna_shiny/git_seurat/seurat_practice_files/top_markers")

#heat_map_list = lapply(top_markers_start:top_markers_end, function(x) {
# topn <- top_markers %>% group_by(cluster) %>% top_n(n = x, wt = avg_logFC)
#DoHeatmap(pbmc, features = topn$gene ) + NoLegend()
#})
#saveRDS(heat_map_list, "~/scrna_shiny/git_seurat/seurat_practice_files/heat_map_list")
heat_map_list = readRDS("~/scrna_shiny/git_seurat/seurat_practice_files/heat_map_list")



ui = fluidPage(theme = shinytheme("lumen"),
               titlePanel("PBMC interactive displays"),
               sidebarLayout(
                 sidebarPanel(
                   radioButtons(inputId = "graph_type", label = strong("Graph types"), choices = c("Principal Component", "UMAP", "Top cluster marker genes against all other clusters", "Marker genes", "Top Genes by cluster")),
                   
                   # Display depending on which graph type is selected
                   conditionalPanel(condition = "input.graph_type == 'Principal Component'",
                                    # Select PC to plot
                                    numericInput(inputId = "pc", label = strong("Principal component heat map"), value = 1, min = 1, max = 15),
                                    sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                   ),
                   conditionalPanel(condition = "input.graph_type == 'UMAP'",
                                    # Select PC to plot
                                    sliderInput(inputId = "UMAP", label = strong("Number of dimensions from 1"), value = 10, min = 5, max = 10),
                                    #actionButton(inputId = "go", label = "Run")
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Top cluster marker genes against all other clusters'",
                                    #Select PC to plot
                                    #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                                    sliderInput(inputId = "no_of_de_genes", label = strong("Choose number of top differentially expressed genes to display:"), value = 10, min = 1, max = 100, step = 1),
                                    #sliderInput(inputId = "cluster_no", label = strong("Number of cluster to display genes for:"), value = 0, min = 0, max = 8, step = 1),
                                    
                                    #selectInput(inputId = "select_cell_type", label = strong("Select cell types to assess effect of stimulation for:"),choices = cell_types, multiple = T)
                                    
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Marker genes'",
                                    #Select PC to plot
                                    selectInput(inputId = "marker_genes", label = strong("Select marker gene:"), choices = all.genes, multiple = T),
                                    #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Top Genes by cluster'",
                                    #Select PC to plot
                                    sliderInput(inputId = "top_marker_genes", label = strong("Number of most differentiated genes to display per cluster:"), min = top_markers_start, max = top_markers_end, value = top_markers_start),
                                    #actionButton(inputId = "go_heatmap", label = "Run")
                                    #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                   )
                 ),
                 
                 # Output: Description, lineplot, and reference
                 mainPanel(  
                   conditionalPanel(
                     condition = "input.graph_type == 'Principal Component'", 
                     plotOutput("heatmap1")),
                   conditionalPanel(
                     condition = "input.graph_type == 'UMAP'", 
                     plotOutput("heatmap")),
                   conditionalPanel(
                     condition = "input.graph_type == 'Top cluster marker genes against all other clusters'", 
                     tableOutput("cluster_marker_table")),
                   conditionalPanel(
                     condition = "input.graph_type == 'Marker genes'", 
                     plotOutput("vln_plot"), plotOutput("marker_umap") ),
                   conditionalPanel(
                     condition = "input.graph_type == 'Top Genes by cluster'", 
                     plotOutput("marker_heatmap"))
                   
                 )
                 #mainPanel(
                 #plotOutput(outputId = "heatmap1", width = "100%", height = "300px"),
                 #plotOutput(outputId = "heatmap"),
                 #plotOutput("marker_plot")
                 #)
               )
)

server = function(input, output) {
  
  output$heatmap1 = renderPlot({
    
    DimHeatmap(pbmc_pca, dims = input$pc, cells = input$cells, balanced = TRUE)
  })
  
  #umap_dimensions = reactive({paste(input$UMAP, collapse = ":")})
  #umap_data = eventReactive(input$go, {
  #RunUMAP(pbmc, dims = 1:input$UMAP)
  #pbmc.umap.list[[input$UMAP-4]]
  #})
  output$heatmap = renderPlot({
    
    DimPlot(pbmc.umap.list[[input$UMAP-4]], reduction = "umap")
    
  })
  
  #reactive({})
  output$cluster_marker_table = renderTable({
    
    top_markers %>% group_by(cluster) %>% top_n(n = input$no_of_de_genes, wt = avg_logFC)
  })
  
  output$vln_plot = renderPlot({
    
    VlnPlot(pbmc, features = input$marker_genes)
  })
  output$marker_umap = renderPlot({
    
    FeaturePlot(pbmc, features = input$marker_genes)
  })
  
  #umap_cluster_heatmap = eventReactive(input$go_heatmap, {
  
  # heat_map_list[[input$top_marker_genes-1]]
  #})
  
  output$marker_heatmap = renderPlot({
    heat_map_list[[input$top_marker_genes-1]]
  })
  
  
}



shinyApp(ui = ui, server = server)
