library(Seurat)
library(cowplot)
library(SeuratData)
library(ggplot2)
library(shiny)
library(shinythemes)

InstallData("ifnb") 

#data("ifnb.SeuratData::ifnb")
ifnb.list = SplitObject(ifnb.SeuratData::ifnb, split.by = "stim")

ifnb.list = lapply(ifnb.list, function(x) { 
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
} )

immune.anchors = FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)

immune.combined = IntegrateData(anchorset = immune.anchors, dims = 1:20)

#Running a single integrated analysis on all cells from the different conditions

DefaultAssay(immune.combined) = "integrated"

##Running the standard workflow
immune.combined = ScaleData(immune.combined, verbose = F)
immune.combined = RunPCA(immune.combined, npcs = 30, verbose = F)

##tSNE and clustering
immune.combined.umap.list = lapply(10:11, function(x) {
  a = RunUMAP(immune.combined, reduction = "pca", dims = 1:x)
  lapply(10:11, function(x) {
    b = FindNeighbors(a, reduction = "pca", dims = 1:x)
    lapply(seq(0.1,0.5, by = 0.1), function(x) FindClusters(b, resolution = x))
  })
  })

#immune.combined = RunUMAP(immune.combined, reduction = "pca", dims = 1:20)

#immune.combined = FindNeighbors(immune.combined, reduction = "pca", dims = 1:10)
#immune.combined.clusters = FindClusters(immune.combined.umap.list[[1]][[1]][[1]], resolution = 0.5)

##visualization








ui = fluidPage(theme = shinytheme("lumen"),
               titlePanel("PBMC interactive displays"),
               sidebarLayout(
                 sidebarPanel(
                   radioButtons(inputId = "graph_type", label = strong("Graph types"), choices = c("UMAP", "Marker genes", "Top Genes by cluster")),
                   
                   # Display depending on which graph type is selected
                   
                   conditionalPanel(condition = "input.graph_type == 'UMAP'",
                                    # Select PC to plot
                                    sliderInput(inputId = "umap_dim", label = strong("Number of UMAP dimensions from 1"), value = 10, min = 10, max = 11, step = 1),
                                    sliderInput(inputId = "neighbours_dim", label = strong("Number of SNN dimensions from 1"), value = 10, min = 10, max = 11, step = 1),
                                    sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = 0.5, min = 0.1, max = 0.5, step = 0.1),
                                    actionButton(inputId = "go", label = "Run")
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Marker genes'",
                                    #Select PC to plot
                                    selectInput(inputId = "marker_genes", label = strong("Select marker gene:"), choices = all.genes, multiple = T),
                                    #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Top Genes by cluster'",
                                    #Select PC to plot
                                    sliderInput(inputId = "top_marker_genes", label = strong("Number of most differentiated genes to display per cluster:"), min = 1, max = 20, value = 3),
                                    actionButton(inputId = "go_heatmap", label = "Run")
                                    #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                   )
                 ),
                 
                 # Output: Description, lineplot, and reference
                 mainPanel(  
                   
                   conditionalPanel(
                     condition = "input.graph_type == 'UMAP'", plotOutput("all_groups"), plotOutput("stim_vs_ctrl")),
                   conditionalPanel(
                     condition = "input.graph_type == 'Marker genes'", plotOutput("vln_plot"), plotOutput("marker_umap") ),
                   conditionalPanel(
                     condition = "input.graph_type == 'Top Genes by cluster'", plotOutput("marker_heatmap"))
                   
                 )
                 #mainPanel(
                 #plotOutput(outputId = "heatmap1", width = "100%", height = "300px"),
                 #plotOutput(outputId = "heatmap"),
                 #plotOutput("marker_plot")
                 #)
               )
)

server = function(input, output) {
  
  #output$all_groups = eventReactive(input$go, {
  output$all_groups = renderPlot({

      p1 = DimPlot(immune.combined.umap.list[[(input$umap_dim - 9)]][[(input$neighbours_dim - 9)]][[(input$clusters_res * 10)]], reduction = "umap", group.by = "stim")
      p2 = DimPlot(immune.combined.umap.list[[(input$umap_dim - 9)]][[(input$neighbours_dim - 9)]][[(input$clusters_res * 10)]], reduction = "umap", label = T)
      plot_grid(p1, p2)
      
    })
  #})
  output$stim_vs_ctrl = renderPlot({
    
    ##Showing the stimulated and control umaps side by side
    DimPlot(immune.combined.umap.list[[(input$umap_dim - 9)]][[(input$neighbours_dim - 9)]][[(input$clusters_res * 10)]], reduction = "umap", split.by = "stim")
  })
  #})
  
 
}  

  #umap_dimensions = reactive({paste(input$UMAP, collapse = ":")})
  umap_data = eventReactive(input$go, {
    RunUMAP(pbmc, dims = 1:input$UMAP)
  })
  output$heatmap = renderPlot({
    
    DimPlot(umap_data(), reduction = "umap", )
    
  })
  
  #reactive({})
  
  output$vln_plot = renderPlot({
    
    VlnPlot(pbmc, features = input$marker_genes)
  })
  output$marker_umap = renderPlot({
    
    FeaturePlot(pbmc, features = input$marker_genes)
  })
  
  umap_cluster_heatmap = eventReactive(input$go_heatmap, {
    
    top_genes = FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>% group_by(cluster) %>% top_n(n = input$top_marker_genes, wt = avg_logFC)
    top_genes = top_genes$gene
  })
  
  output$marker_heatmap = renderPlot({
    DoHeatmap(pbmc, features = umap_cluster_heatmap()) + NoLegend()
  })
  
  
}



shinyApp(ui = ui, server = server)
