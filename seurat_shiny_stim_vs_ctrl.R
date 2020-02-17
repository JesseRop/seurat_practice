library(Seurat)
library(cowplot)
#library(SeuratData)
library(ggplot2)
library(shiny)
library(shinythemes)
library(plotly)

all_immune_genes = readRDS("/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/all_immune_genes.rds")
immune.combined.umap.list = readRDS("/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/immune.combined.umap.list.rds")
immune.combined.clusters.tables = readRDS("/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/immune.combined.clusters.tables.rds")

if (!(all(exists("all_immune_genes"), exists("immune.combined.umap.list"), exists("immune.combined.clusters.tables")))) {
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
  
  all_immune_genes = row.names(immune.combined)
  saveRDS(all_immune_genes, "/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/all_immune_genes.rds")
  

  ##tSNE and clustering
  immune.combined.umap.list = lapply(dim1:dim2, function(x) {
    a = RunUMAP(immune.combined, reduction = "pca", dims = 1:x)
    lapply(dim1:dim2, function(x) {
      b = FindNeighbors(a, reduction = "pca", dims = 1:x)
      lapply(seq(res1, res2, by = 0.1), function(x) FindClusters(b, resolution = x))
    })
  })
  
  saveRDS(immune.combined.umap.list, "/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/immune.combined.umap.list.rds")
  
  immune.combined.clusters.tables = lapply(immune.combined.umap.list, function(x) { 
    lapply(x, function(x) {
      lapply(x, function(x){
        DefaultAssay(x) = "RNA"
        #cluster.markers = lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
        cluster.markers = lapply(6:7, function(y) {
          FindConservedMarkers(x, ident.1 = y, grouping.var = "stim")
          
        })
      }) 
    }) 
  })
  
  saveRDS(immune.combined.clusters.tables, "/home/jr345y/scrna_shiny/git_seurat/seurat_practice_files/immune.combined.clusters.tables.rds")
  
}

cell_types = c( "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T")
dim1 = 19
dim2 = 20

res1 = 0.4
res2 = 0.5
##visualization


ui = fluidPage(theme = shinytheme("lumen"),
               titlePanel("Interferron beta stimulation of PBMCs"),
               sidebarLayout(
                 sidebarPanel(
                   radioButtons(inputId = "graph_type", label = strong("Graph types"), choices = c("UMAP", "Top cluster marker genes conserved between conditions", "Labelled populations", "Differentially expressed genes after IFNb stimulation", "Visualize genes in IFNb stimulated VS Control cells")),
                   
                   # Display depending on which graph type is selected
                   
                   conditionalPanel(condition = "input.graph_type == 'UMAP'",
                                    # Select PC to plot
                                    sliderInput(inputId = "umap_dim", label = strong("Number of UMAP dimensions from 1"), value = dim2, min = dim1, max = dim2, step = 1),
                                    sliderInput(inputId = "neighbours_dim", label = strong("Number of SNN dimensions from 1"), value = dim2, min = dim1, max = dim2, step = 1),
                                    sliderInput(inputId = "clusters_res", label = strong("Louvain algorithm resolution"), value = res2, min = res1, max = res2, step = 0.1),
                                    actionButton(inputId = "go", label = "Run")
                                    
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Top cluster marker genes conserved between conditions'",
                   
                                    sliderInput(inputId = "marker_genes_no", label = strong("Choose number of top markers to display:"), value = 10, min = 1, max = 100, step = 1),
                                    sliderInput(inputId = "marker_genes_cluster", label = strong("Choose cluster to display markers for"), value = 7, min = 6, max = 7, step = 1),
                                    actionButton(inputId = "go_marker", label = "Search"),
                                    selectInput(inputId = "select_markers_umap", label = strong("Select conserved markers:"), choices = all_immune_genes, multiple = T)
                                    
                    ),
                   conditionalPanel(condition = "input.graph_type == 'Labelled populations'",
                                    #Select PC to plot
                                    actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                                    selectInput(inputId = "select_markers_dotplot", label = strong("Select markers for dotplot:"),choices = all_immune_genes, multiple = T)
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Differentially expressed genes after IFNb stimulation'",
                                    #Select PC to plot
                                    #actionButton(inputId = "go_labelled_umap", label = "View labelled clusters"),
                                    sliderInput(inputId = "no_of_de_genes", label = strong("Choose number of top differentially expressed genes to display:"), value = 10, min = 1, max = 100, step = 1),
                                    selectInput(inputId = "select_cell_type", label = strong("Select cell types to assess effect of stimulation for:"),choices = cell_types, multiple = T)
                                    
                   ),
                   conditionalPanel(condition = "input.graph_type == 'Visualize genes in IFNb stimulated VS Control cells'",
                                    #Select PC to plot
                                    selectInput(inputId = "de_genes", label = strong("Choose gene:"),choices = all_immune_genes, multiple = T),
                                    #actionButton(inputId = "go_heatmap", label = "Run"),
                                    #sliderInput(inputId = "cells", label = strong("Select number of cells in each tail"), min = 1, max = 500, value = 50)
                   )
                 ),
                 
                 # Output: Description, lineplot, and reference
                 mainPanel(  
                   
                   conditionalPanel(
                     condition = "input.graph_type == 'UMAP'", 
                     plotOutput("all_groups"), plotOutput("stim_vs_ctrl")),
                   conditionalPanel(
                     condition = "input.graph_type == 'Top cluster marker genes conserved between conditions'", 
                     tableOutput("top_conserved_genes"), plotOutput("conserved_markers_umap")), 
                   conditionalPanel(
                     condition = "input.graph_type == 'Labelled populations'", 
                     plotOutput("labelled_umap"), plotOutput("marker_dotplot")),
                   conditionalPanel(
                     condition = "input.graph_type == 'Differentially expressed genes after IFNb stimulation'", 
                     tableOutput("top_de_genes"), plotOutput("cell_type_plot", hover = hoverOpts(id ="plot_hover")), verbatimTextOutput("hover_info")),
                   conditionalPanel(
                     condition = "input.graph_type == ''Visualize genes in IFNb stimulated VS Control cells'", plotOutput("de_stim_vs_ctrl_um"), plotOutput("de_stim_vs_ctrl_vp") )
                   #conditionalPanel(
                     #condition = "input.graph_type == 'Top Genes by cluster'", plotOutput("marker_heatmap"))
                   
                 )
                 )
)

server = function(input, output) {
  
  umap_clusters = eventReactive(input$go, {
    immune.combined.umap.list[[(input$umap_dim - 18)]][[(input$neighbours_dim - 18)]][[((input$clusters_res * 10)-3)]]
  })
  
  output$all_groups = renderPlot({

      p1 = DimPlot(umap_clusters(), reduction = "umap", group.by = "stim")
      p2 = DimPlot(umap_clusters(), reduction = "umap", label = T)
      plot_grid(p1, p2)
      
    })
  
  ##Showing the stimulated and control umaps side by side
  output$stim_vs_ctrl = renderPlot({
    DimPlot(umap_clusters(), reduction = "umap", split.by = "stim")
  })

  cluster_markers = eventReactive(input$go_marker, {
    head(immune.combined.clusters.tables[[(input$umap_dim - 18)]][[(input$neighbours_dim - 18)]][[((input$clusters_res * 10)-3)]][[(input$marker_genes_cluster - 5)]], n = input$marker_genes_no)
    })
  
  ##Finding conserved genes in clusters in both conditions to annotate cell types
  output$top_conserved_genes = renderTable({
    #isolate({
    cluster_markers()
  },rownames = T)
  #})
  #})
  
  umap_cluster_modified_rna = reactive({
    umap_cluster_modified_ul = umap_clusters()
    DefaultAssay(umap_cluster_modified_ul) = "RNA"
    umap_cluster_modified_ul
  })
  
  output$conserved_markers_umap = renderPlot({
    FeaturePlot(umap_cluster_modified_rna(), features = input$select_markers_umap, min.cutoff = "q9")
  })

  
  umap_cluster_modified_renamed = reactive({
    umap_cluster_modified = umap_cluster_modified_rna()
    umap_cluster_modified <- RenameIdents(umap_cluster_modified, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", `10` = "Mk", `11` = "pDC", `12` = "Eryth")
  })
  
  umap_cluster_modified_umap = eventReactive(input$go_labelled_umap, {
    DimPlot(umap_cluster_modified_renamed(), label = TRUE)
  })
  output$labelled_umap = renderPlot({
    umap_cluster_modified_umap()
  })
  
  umap_cluster_modified_ren_reo = reactive({
    umap_cluster_modified_reo = umap_cluster_modified_renamed()
    Idents(umap_cluster_modified_reo) <- factor(Idents(umap_cluster_modified_reo), levels = c( "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))
    umap_cluster_modified_reo

  })
  
  output$marker_dotplot = renderPlot({
    DotPlot(umap_cluster_modified_ren_reo(), features = input$select_markers_dotplot, cols = c("blue", "red"), dot.scale = 6, split.by = "stim") + RotatedAxis()

  })
 
  stim_markers = reactive({
    umap_cluster_modified = umap_cluster_modified_ren_reo()
    umap_cluster_modified$celltype.stim <- paste(Idents(umap_cluster_modified), umap_cluster_modified$stim, sep = "_")
    umap_cluster_modified$celltype <- Idents(umap_cluster_modified)
    Idents(umap_cluster_modified) <- "celltype.stim"
    umap_cluster_modified
  })
  
  output$top_de_genes = renderTable({
    
    ##Finding conserved genes in clusters in both conditions to annotate cell types
    #isolate({
    b.interferon.response <- FindMarkers(stim_markers(), ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
    head(b.interferon.response, n = input$no_of_de_genes)
  },rownames = T)
  
  
  cell_type_de = reactive({
    theme_set(theme_cowplot())
  
    cells <- subset(umap_cluster_modified_ren_reo(), idents = input$select_cell_type)
    Idents(cells) <- "stim"
    avg.cells <- log1p(AverageExpression(cells, verbose = FALSE)$RNA)
    avg.cells$gene <- rownames(avg.cells)
    avg.cells
    
  })
  output$cell_type_plot = renderPlot({
    ggplot(data=cell_type_de(), aes(CTRL , STIM)) + geom_point() + ggtitle(input$select_cell_type)
    #ggplotly(p)
  })
  
  displayed_text <- reactive({
    req(input$plot_hover)
    nearPoints(cell_type_de(), input$plot_hover)
    
  })
  
  output$hover_info <- renderPrint({
    req(displayed_text())
    
    cat("Name\n")
    displayed_text()
  })
  
  output$de_stim_vs_ctrl_um = renderPlot({
    
    FeaturePlot(stim_markers(), features = input$de_genes, split.by = "stim", max.cutoff = 3,cols = c("grey", "red"))
  })
  
  output$de_stim_vs_ctrl_vp = renderPlot({
    
    VlnPlot(stim_markers(), features = input$de_genes, split.by = "stim", group.by = "celltype", pt.size = 0, combine = FALSE)  
    
    })

  
}


shinyApp(ui = ui, server = server)
