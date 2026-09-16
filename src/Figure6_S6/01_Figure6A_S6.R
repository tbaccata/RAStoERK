source(here::here("src/utils/utils_path.R"))

library(viridis)
library(RColorBrewer)
source(file.path(srcDir(), "utils/visualization_utils.R"))

mergedQuant <- read.table(file.path(interimDataDir(), 'all_interactions_apms_turboid.tsv'),
                          header = T, sep = '\t')

pi3kakt <- 'BRAF
RAF1
PIK3CA
PIK3CD
PIK3R1
PIK3R2
PIK3R3
AKT3
DEPTOR
IRS2'

pi3kakt <- (unlist(strsplit(pi3kakt, '\n') ) )
pi3kakt <- pi3kakt[!duplicated(pi3kakt)]

for (method in c('AP-MS', 'TbID')) {
  
  wntDpData <- prepDotplotData(mergedQuant[grep('ras', mergedQuant$bait ),],
                               pi3kakt, method = method  )
  
  toConsider <- unique(wntDpData$Gene[!is.na(wntDpData$padj) &
                                        wntDpData$padj < 0.05 & wntDpData$log2FC>0])
  wntDpData <- wntDpData[wntDpData$Gene %in% toConsider,]
  
  bait2prey <- apmsPreys
  if (method == 'AP-MS') {
    bait2prey <- apmsPreys
  } else {
    bait2prey <- tidPreys
  }
  
  for (bait in unique(wntDpData$Group)) {
    
    tmp <- wntDpData[wntDpData$Group==bait,]
    tmp$significant <- ifelse(
      tmp$Gene%in%bait2prey[[bait]], '1.5', '0'
    )
    
    wntDpData$significant[wntDpData$Group==bait ] <- tmp$significant
  }
  
  
  gdf <- data.frame(comparison=wntDpData$Comparison, group=wntDpData$Group)
  gdf <- gdf[!duplicated(gdf$comparison),]
  
  idx <- match(baitOrder, gdf$group )
  gdf <- gdf[idx,]
  
  wntDpData$Gene <- factor(wntDpData$Gene,
                           levels = rev(pi3kakt))
  
  p <- plotDotplot(wntDpData,
                   gdf,
                   dotplotColors=colorRampPalette(brewer.pal(8, "Blues"))(50), # viridis(100), #brewer.pal(8, "Blues"),
                   maxColorGradient=5,
                   minSizeGradient=1,
                   maxSizeGradient=5,
                   clusteringMetric = "log2FC",
                   dotplot_distance_metric = "euclidean",
                   dotplot_clustering_method = "complete",
                   dotplot_cluster_rows=F,
                   dotplot_cluster_columns=F,
                   sigCutoffValue="adj.p-value",
                   dotplot_ctrl_substraction=T, 
                   show_legend = F )
  p <- p + ggtitle(method)
  
  ggsave(
    file.path(figureDir(), paste0('Figure6/6A_pi3k_xtalk_', method, '_panras.pdf') ),
    p,
    width = 4,
    height = 3.5,
    device = cairo_pdf
  )
}

toFile <- mergedQuant[grep('Nras', mergedQuant$bait ),]
toFile <- toFile[toFile$Gene.names%in%pi3kakt,]

write.table(
  toFile[, grep("combined", names(toFile), invert = T )],
  file.path(figureDir(), 'Figure6/data/nras_pi3k_xtalk_interactome.tsv'),
  row.names = F, sep = '\t', quote = F
)

#
library(grid)
p <- plotDotplot(wntDpData,
                 gdf,
                 dotplotColors=colorRampPalette(brewer.pal(8, "Blues"))(50),
                 maxColorGradient=5,
                 minSizeGradient=1,
                 maxSizeGradient=5,
                 clusteringMetric = "log2FC",
                 dotplot_distance_metric = "euclidean",
                 dotplot_clustering_method = "complete",
                 dotplot_cluster_rows=F,
                 dotplot_cluster_columns=F,
                 sigCutoffValue="adj.p-value",
                 dotplot_ctrl_substraction=T, 
                 show_legend = T )
legend <- cowplot::get_legend(p)

ggsave(plot = legend,
       filename = file.path(figureDir(), 'Figure6/dotplot_legends.pdf'),
       device = cairo_pdf,
       width = 2.46, height = 3.42
)
