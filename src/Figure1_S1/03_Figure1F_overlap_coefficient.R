#' @author Sebastian Didusch
#' @date 25.06.11
#' 
#' Calculate pairwise overlaps between prey lists

source(here::here("src/utils/utils_path.R"))

# helper functions
prepareOverlapMatrix <- function(overlap, groupOrder) {
  df <- overlap[, c("sample1", "sample2", "overlap")]
  df1 <-
    data.frame(
      sample1 = df$sample2,
      sample2 = df$sample1,
      overlap = df$overlap
    )
  limits <- c(min(df[, 3], na.rm = T) - 0.05, 1)
  names(df) <- c("sample1", "sample2", "overlap")
  df <- rbind(df, df1)
  overlap_matrix <-
    acast(df, sample1 ~ sample2, value.var = "overlap")
  
  groupOrder <- baitOrder[baitOrder %in% groupOrder ]
  overlap_matrix[groupOrder, groupOrder]
}
# data
apmsOverlap <- calc_pairwise_overlaps(apmsPreys)
tidOverlap <- calc_pairwise_overlaps(tidPreys)

# heatmap annotation df
annotDesign <- apmsDesign[apmsDesign$batch==1 & apmsDesign$type != 'Control', c('groups', 'type', 'node') ]
rownames(annotDesign) <- annotDesign$groups
annotDesign$groups <- NULL

annotDesign <- annotDesign[baitOrder,]

### clustering
apmsOverlapMatrix <- prepareOverlapMatrix(apmsOverlap, 
                                          apmsDesign$groups[apmsDesign$type!='Control'])
tidOverlapMatrix <- prepareOverlapMatrix(tidOverlap, 
                                         apmsDesign$groups[apmsDesign$type!='Control'])
colBreaks <- seq(0.1, 0.8, 0.02)
'#0c607f'

apmsColorSpace <- '#c7e9b4'
# ap-ms
pdf(file.path(figureDir(), 'Figure1/overlapcoeff_apms_clustered.pdf'  ),
    height = 4.85, width = 6.25)
pheatmap(apmsOverlapMatrix,
         main = 'Overlap coefficient (AP-MS)',
         cluster_rows = T,
         cluster_cols = T,
         border_color = 'white',
         scale = 'none',
         treeheight_col = 5,
         treeheight_row = 5,
         annotation_row = annotDesign, annotation_colors = annotation_colors,
         color = colorRampPalette(c('white', '#0c607f'))(length(colBreaks)),
         #color = colorRampPalette(brewer.pal(8, 'BuPu'))(length(colBreaks)),
         breaks = seq(0.1, 0.8, 0.02) )
dev.off()
# tbid
pdf(file.path(figureDir(), 'Figure1/overlapcoeff_turboid_clustered.pdf'  ),
    height = 4.85, width = 6.25)
pheatmap(tidOverlapMatrix,
         main = 'Overlap coefficient (TbID)',
         cluster_rows = T,
         cluster_cols = T,
         border_color = 'white',
         treeheight_col = 5,
         treeheight_row = 5,
         scale = 'none',
         annotation_row = annotDesign, annotation_colors = annotation_colors,
         color = colorRampPalette(c('white', '#f05455'))(length(colBreaks)),
         #color = colorRampPalette(brewer.pal(8, 'Blues'))(length(colBreaks)),
         breaks = seq(0.1, 0.8, 0.02) )
dev.off()

write.table(
  apmsOverlapMatrix,
  file.path(figureDir(), 'Figure1/data/overlapcoeff_apms.tsv'),
  sep = '\t', quote = F
)

write.table(
  tidOverlapMatrix,
  file.path(figureDir(), 'Figure1/data/overlapcoeff_turboid.tsv'),
  sep = '\t', quote = F
)
