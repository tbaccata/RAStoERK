source(here::here("src/utils/utils_path.R"))

library(RColorBrewer)
source(file.path(srcDir(), "utils/visualization_utils.R"))

mergedQuant <- read.table(file.path(interimDataDir(), 'all_interactions_apms_turboid.tsv'),
                          header = T, sep = '\t')

rasEffectors <- 'ARAF
BRAF
RAF1
PIK3CA
PIK3CD
PIK3R1
PIK3R2
PIK3R3
RALGDS
RGL1
RGL2
AFDN
PLCE1
SNX27
ARAP3
RASSF5
RAPGEF2
RAPGEF4
RAPGEF6
KRIT1
RADIL'

rasEffectors <- strsplit(rasEffectors, '\n') %>% unlist()

for (method in c('AP-MS', 'TbID')) {
  effectorData <- prepDotplotData(mergedQuant[grep('ras', mergedQuant$bait ),],
                                  rasEffectors, method = method  )
  
  lvls <- baitOrder[baitOrder %in% unique(effectorData$Group)]
  effectorData$Gene <- factor(effectorData$Gene, levels = rev(rasEffectors))
  
  ###
  bait2prey <- apmsPreys
  if (method == 'AP-MS') {
    bait2prey <- apmsPreys
  } else {
    bait2prey <- tidPreys
  }
  
  for (bait in unique(effectorData$Group)) {
    
    tmp <- effectorData[effectorData$Group==bait,]
    tmp$significant <- ifelse(
      tmp$Gene%in%bait2prey[[bait]], '1.5', '0'
    )
    
    effectorData$significant[effectorData$Group==bait ] <- tmp$significant
  }
  ###
  
  gdf <- data.frame(comparison=effectorData$Comparison, group=effectorData$Group)
  gdf <- gdf[!duplicated(gdf$comparison),]
  
  idx <- match(lvls, gdf$group )
  gdf <- gdf[idx,]
  
  p <- plotDotplot(effectorData,
                   gdf,
                   dotplotColors=colorRampPalette(brewer.pal(8, "Blues"))(50), #viridis(100), #brewer.pal(8, "Blues"),
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
  p <- p + ggtitle(method)
  
  ggsave(
    file.path(figureDir(), paste0('S6/raseffectors_', method, '_blue.pdf') ),
    p,
    width = 6,
    height = 5,
    device = cairo_pdf
  )
  #
  write.table(
    effectorData,
    file.path(figureDir(), paste0('S6/data/raseffectors_', method, '_data.tsv') ),
    row.names = F,
    sep = '\t',
    quote = F
  )
}
