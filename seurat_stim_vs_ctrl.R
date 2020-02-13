library(Seurat)
library(cowplot)
library(SeuratData)
library(ggplot2)

InstallData("ifnb") 

#data("ifnb.SeuratData::ifnb")
ifnb.list = SplitObject(ifnb.SeuratData::ifnb, split.by = "stim")

ifnb.list = lapply(ifnb.list, function(x) { 
  x = NormalizeData(x)
  x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  } )
alist = lapply(1:2, function(x) {
  #x
  lapply(1:2, function(y) {
    #y+2
    lapply(1:2, function(j){
      j^3
    })
    })
  })

interpretor = lapply(1:2, function(x) {
  a = RunUMAP(immune.combined, reduction = "pca", dims = 1:x)
  lapply(dim1:dim2, function(x) {
    b = FindNeighbors(a, reduction = "pca", dims = 1:x)
    lapply(seq(res1, res2, by = 0.1), function(x) FindClusters(b, resolution = x))
  })
})

immune.anchors = FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)

immune.combined = IntegrateData(anchorset = immune.anchors, dims = 1:20)

#Running a single integrated analysis on all cells from the different conditions

DefaultAssay(immune.combined) = "integrated"

##Running the standard workflow
immune.combined = ScaleData(immune.combined, verbose = F)
immune.combined = RunPCA(immune.combined, npcs = 30, verbose = F)

##tSNE and clustering

immune.combined = RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined = FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined = FindClusters(immune.combined, resolution = 0.5)

##visualization

p1 = DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 = DimPlot(immune.combined, reduction = "umap", label = T)
plot_grid(p1, p2)

##Showing the stimulated and control umaps side by side

DimPlot(immune.combined, reduction = "umap", split.by = "stim")

##Finding genes with conserved expression in both conditions in a particular cell-type (cluster7)
DefaultAssay(immune.combined) = "RNA"
nk.markers = FindConservedMarkers(immune.combined, ident.1 = 7, grouping.var = "stim", verbose = F)
head(nk.markers)

FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), min.cutoff = "q9")

immune.combined = RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

DimPlot(immune.combined, label = TRUE)

Idents(immune.combined)  = factor(Idents(immune.combined), levels = c("Mono/Mk Doublets", "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated", "CD4 Naive T", "CD4 Memory T"))
markers.to.plot = c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5", "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1", "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ")

DotPlot(immune.combined, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, split.by = "stim") + RotatedAxis()

##Differences in stimulated versus unstimulated conditions
theme_set(theme_cowplot())
t.cells = subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) = "stim"
avg.t.cells = log1p(AverageExpression(t.cells, verbose = F)$RNA)
avg.t.cells$gene = row.names(avg.t.cells)

cd14.mono = subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) = "stim"
avg.mono.cells = log1p(AverageExpression(cd14.mono, verbose = F)$RNA)
avg.mono.cells$gene = row.names(avg.mono.cells)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 = ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 = LabelPoints(plot = p1, points = genes.to.label, repel = T)

p2 = ggplot(avg.mono.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 = LabelPoints(plot = p1, points = genes.to.label, repel = T)

plot_grid(p1, p2)

##Differences in stim vs unstim in all the clusters (cell types)
immune.combined$celltype.stim = paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype = Idents(immune.combined)

Idents(immune.combined) = "celltype.stim"

b.interferon.response = FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = F)
head(b.interferon.response, n = 15)

##Visualizations

FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))

plots = VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype", pt.size = 0, combine = F)
CombinePlots(plots = plots, ncol = 1)


