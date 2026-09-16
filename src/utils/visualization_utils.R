library(tidyr)

#' dotplot
#' 
#' @param dataDotplot data.frame, output from getDotplotData.
#' @param dotplotGroupsDf data.frame, containing group comparisons and corresponding group for dotplot.
#' @param dotplotColors vector of colors.
#' @param minColorGradient numeric.
#' @param maxColorGradient numeric.
#' @param minSizeGradient int.
#' @param maxSizeGradient int.
#' @param clusteringMetric char, one of "AvgIntensity" and "log2FC".
#' @param dotplot_distance_metric char, valid option from hclust.
#' @param dotplot_clustering_method char, valid option from hclust.
#' @param dotplot_cluster_columns boolean, performs clustering on columns.
#' @param sigCutoffValue char, choice between "adj.P.Val" and "p-value".
#' @param dotplot_ctrl_substraction boolean, ignores proteins with negative log2FCs.
#' @return ggplot2 object.
#' @examples
plotDotplot <-
  function(dataDotplot,
           dotplotGroupsDf,
           dotplotColors=NULL,
           minColorGradient=NULL,
           maxColorGradient=NULL,
           minSizeGradient=1,
           maxSizeGradient=7,
           clusteringMetric = "log2FC",
           dotplot_distance_metric = "canberra",
           dotplot_clustering_method = "complete",
           dotplot_cluster_rows=TRUE,
           dotplot_cluster_columns=TRUE,
           sigCutoffValue="adj.p-value",
           dotplot_ctrl_substraction=TRUE,
           show_legend=TRUE) {
    
    if (is.null(dotplotColors))
      dotplotColors <- viridis(20, option = "viridis")
    
    dataDotplot <- dataDotplot[!is.na(dataDotplot[[clusteringMetric]]),]
    
    mat <- dataDotplot %>%
      dplyr::select(Gene, Group, all_of(clusteringMetric)) %>%  # drop unused columns to faciliate widening
      pivot_wider(names_from = Group, values_from = all_of(clusteringMetric)) %>%
      data.frame() # make df as tibbles -> matrix annoying
    
    mat <- mat[!is.na(mat$Gene),]
    row.names(mat) <- mat$Gene  # put gene in `row`
    mat$Gene <- NULL #drop gene column as now in rows
    
    tmat <- mat %>% as.matrix()
    tmat[is.na(tmat)] <- 0
    
    dist_meth <- ifelse(!is.null(dotplot_distance_metric),
                        dotplot_distance_metric, 
                        "canberra")
    clst_meth <- ifelse(!is.null(dotplot_clustering_method),
                        dotplot_clustering_method,
                        "complete")
    
    if (!is.null(dotplot_cluster_columns) && dotplot_cluster_columns) { # order by clustering
      cclust <- hclust(dist(t(tmat), method = dist_meth),
                       method = clst_meth) # hclust with distance matrix
      cclust <- reorder(as.dendrogram(cclust), colMeans(mat))
      dataDotplot$Group <-
        factor(dataDotplot$Group, levels = names(mat)[as.hclust(cclust)$order])
    } else { # order by input
      dataDotplot$Group <-
        factor(dataDotplot$Group, levels = dotplotGroupsDf$group)
    }
    
    if (dotplot_cluster_rows) {
      rclust <- hclust(dist(tmat, 
                            method = dist_meth), 
                       method = clst_meth) # hclust with distance matrix
      rclust <- reorder(as.dendrogram(rclust), rowMeans(mat))
      dataDotplot$Gene <-
        factor(dataDotplot$Gene, levels = rownames(mat)[rev(as.hclust(rclust)$order)])
    } 
    # else {
    #   dataDotplot$Gene <-
    #     factor(dataDotplot$Gene, levels = unique(dataDotplot$Gene))
    # }
    
    
    signficantTitle <- "adj.p-value"
    if (sigCutoffValue == "p-value")
      signficantTitle <- "p-value"
    
    filterValues <- TRUE
    if (clusteringMetric == "log2FC" &&
        !is.null(dotplot_ctrl_substraction)) {
      filterValues <- dotplot_ctrl_substraction
    }
    
    preFiltered <- dataDotplot
    if (filterValues) {
      preFiltered <- dataDotplot %>% filter(!!rlang::sym(clusteringMetric) > 0)
    }
    
    x <- preFiltered[[clusteringMetric]]
    
    legName <- "Relative AvgIntensity"
    if (clusteringMetric == "log2FC") {
      preFiltered$relativeAbundance <- (x-min(x,na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T))
      legName <- "Relative Fold Change"
    } else {
      preFiltered$relativeAbundance <- x/max(x, na.rm = T)
    }
    
    minColorGradient <-
      ifelse(is.null(minColorGradient),
             min(preFiltered$log2FC, na.rm = T),
             minColorGradient)
    
    maxColorGradient <-
      ifelse(is.null(maxColorGradient),
             max(preFiltered$log2FC, na.rm = T),
             maxColorGradient)
    
    p <- preFiltered %>%
      ggplot(aes(
        x = Group,
        y = Gene,
        fill = log2FC,
        color = significant,
        size = relativeAbundance,
        stroke = 1
      )) +
      geom_point(shape = 21) +
      cowplot::theme_cowplot(font_family = 'Helvetica') +
      theme(axis.line  = element_blank()) +
      theme(axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      )) + 
      theme(axis.text.y = element_text(family = "Helvetica"))
    
    if (show_legend) {
      p <- p + theme(legend.justification = "top") 
    } else {
      p <- p + theme(legend.position = "none")
    }
    
    p <- p +
      ylab('') + xlab('') +
      scale_x_discrete(position = "top") +
      scale_size_continuous(
        range = c(minSizeGradient, maxSizeGradient),
        name = legName, #paste("Relative", clusteringMetric),
        guide = guide_legend(order=2, override.aes = list(shape = 19)),
        breaks = c(
          min(preFiltered$relativeAbundance ) + 0.01,
          max(preFiltered$relativeAbundance )
        ),
        labels = c('', '')
      ) +
      scale_fill_gradientn(
        colours = dotplotColors,
        limits = c(
          ifelse(
            minColorGradient < min(preFiltered$log2FC, na.rm = T),
            min(preFiltered$log2FC, na.rm = T),
            minColorGradient
          ),
          maxColorGradient
        ),
        oob = scales::squish,
        name = 'Log2FC',
        guide = guide_colorbar(order = 1)
      ) +
      scale_color_manual(
        signficantTitle,
        values = c('skyblue', 'black'),
        limits = c('0', '1.5'),
        labels = c('> 0.05', '\u2264 0.05')
      )
    p
  }



prepDotplotData <- function(mergedQuant, preys, method = 'AP-MS') {
  
  if (!method%in%c('AP-MS', 'TbID') ) {
    stop('No proper method. Allowed methods: "AP-MS" and "TbID".')
  }
  
  filtData <- mergedQuant[mergedQuant$Gene.names %in% preys,]
  filtData$Group <- filtData$bait
  filtData$Comparison <- paste0(filtData$Group, '__vs__GFP8ng')
  filtData$AvgIntensity <- rnorm(nrow(filtData), 10, 1)
  
  filtData$Gene <- filtData$Gene.names
  # filtData$ProteinID <- make.names(filtData$Gene.names, unique = T)
  
  if (method == 'AP-MS' ) {
    filtData$log2FC <- filtData$logFC_APMS
    filtData$padj <- filtData$padj_APMS
  } else {
    filtData$log2FC <- filtData$logFC_TbID
    filtData$padj <- filtData$padj_TbID
  }
  
  filtData$significant <- ifelse(
    filtData$padj < 0.05, '1.5', '0'
  )
  
  cols <- c("Group","Comparison","log2FC","padj","AvgIntensity", "Gene" ,"significant" )
  
  filtData[, cols]
}
