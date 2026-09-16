#' @author Sebastian Didusch
#' @date 25.09.10
#' UMAP visualize interactome on Organelle IP data
#' 

## define paths
source(here::here("src/utils/utils_path.R"))

combinedPreys <- list()
#
for (bait in names(apmsPreys) ) {
  combinedPreys[[bait]] <- unique(c(apmsPreys[[bait]], tidPreys[[bait]] ))
}

ips <- read.table(file.path(dataDir(), 'integration/subcell_localization/organelle_ip_umap.tsv'),
                  header = T, sep = '\t')

# only highlight canonical swissprot id
ips$Gene_name_canonical[ips$Majority.protein.IDs=='A0A024RAV5;A0A3G1LBH1;Q14015;Q14014'] <- 'KRAS_trembl'
ips$Gene_name_canonical[ips$Majority.protein.IDs%in%c('H7C560;A0A2R8Y8E0', 'D7RF68;A0A2R8Y492')] <- 'BRAF_trembl'
ips$Gene_name_canonical[ips$Majority.protein.IDs%in%c('A0A384P5S9;A0A024R889;Q9UG16', 'A0A0D9SGF6;B4DGT1')] <- 'SPTAN1_trembl'

# tmpIps <- ips[, c('Gene_name_canonical', 'graph_localization_annotation') ]
# tmpIps <- tmpIps[!duplicated(tmpIps$Gene_name_canonical), ]
# 
# tmp <- merge(interactome, 
#              tmpIps, 
#              by.x = 'Gene.names',
#              by.y = 'Gene_name_canonical' )
# names(tmp)[names(tmp) == 'graph_localization_annotation' ] <- 'subcell_localization'

for (protein in names(annotation_colors$paralog)) {
  baits <- apmsDesign$groups[apmsDesign$paralog == protein]
  baitLocDf <- tmp[tmp$bait%in%baits,]
  sort(table(baitLocDf$subcell_localization[!duplicated(baitLocDf$Gene.names )]) )
}

# write.table(tmp,
#             file.path(resultsDir, 
#                       'subcellular_localization/organelleip_umaps/interactome_in_leonetti.tsv' ),
#             row.names = F, sep = '\t', quote = F
# )

genOrgenellPlot <- function(expDesign, listOfPreys, method, protein, toHighlight=NULL) {
  
  baits <- expDesign$groups[expDesign$paralog==protein]
  baits <- baitOrder[baitOrder%in%baits]
  pattern <- paste(baits, collapse = '|')
  
  umapPlot <- ggplot(ips, aes(x=umap_1, y=umap_2, label=Gene_name_canonical)) +
    geom_point(data = ips[!ips$Gene_name_canonical %in% 
                            unlist(listOfPreys[grep(pattern, names(listOfPreys) )] ),],
               color='grey80', size=0.25, alpha=0.15
    ) +
    geom_point(data = ips[ips$Gene_name_canonical%in% listOfPreys[[baits[1] ]], ],
               color='grey10', alpha=0.33, size=1.75, pch=19 ) +
    geom_point(data = ips[ips$Gene_name_canonical%in% listOfPreys[[baits[2] ]], ],
               color='#1f78b4', alpha=0.33, size=1.75, pch=19  )
  
  if (length(baits) > 3) {
    activePreys <- unique(c(listOfPreys[[baits[3] ]], listOfPreys[[baits[4] ]]))
    umapPlot <- umapPlot + 
      geom_point(data = ips[ips$Gene_name_canonical %in% activePreys ,],
                 color='#e31a1c',alpha=0.33, size=1.5, pch=19  )
  } else {
    umapPlot <- umapPlot + 
      geom_point(data = ips[ips$Gene_name_canonical %in% listOfPreys[[baits[3] ]] ,],
               color='#e31a1c',alpha=0.33, size=1.75, pch=19  )
  }
  
  if (method == 'Combined') {
    
    proteinInteractome <- interactome[grep(pattern, interactome$bait), ]
    coreShared <- unique(proteinInteractome$Gene.names[proteinInteractome$significant%in%c('core', 'shared')] )
    apmsSpec <- unique(proteinInteractome$Gene.names[proteinInteractome$significant=='AP-MS'] )
    tidSpec <- unique(proteinInteractome$Gene.names[proteinInteractome$significant=='TbID'] )
    
    apmsSpec <- setdiff(apmsSpec, c(tidSpec, coreShared) )
    tidSpec <- setdiff(tidSpec, c(apmsSpec, coreShared) )
    
    umapPlot <- ggplot(ips, aes(x=umap_1, y=umap_2, label=Gene_name_canonical)) +
      geom_point(data = ips[!ips$Gene_name_canonical %in% 
                              unlist(listOfPreys[grep(pattern, names(listOfPreys) )] ),],
                 color='grey80', size=0.25, alpha=0.15
      ) +
      geom_point(data = ips[ips$Gene_name_canonical%in% apmsSpec, ],
                 color='#025e7d', alpha=0.5, size=1.5, pch=19 ) +
      geom_point(data = ips[ips$Gene_name_canonical%in% tidSpec, ],
                 color='#f05354', alpha=0.5, size=1.5, pch=19  ) +
      geom_point(data = ips[ips$Gene_name_canonical %in% coreShared, ],
                 color='grey50',alpha=0.5, size=1.5, pch=19 )
  }
  
  if (!is.null(toHighlight)) {
    umapPlot <- umapPlot + 
      geom_text_repel(data = ips[ips$Gene_name_canonical%in%toHighlight,],
                      show.legend = F)
  }
  
  umapPlot
}

#dir.create(
#  file.path(figureDir(), 'S4/'), showWarnings = F)
#interactome
pList <- list()
expDesign <- apmsDesign[apmsDesign$batch==1,]
for (protein in names(annotation_colors$paralog)) {
  pname <- protein
  if (protein == 'ERK1') {
    pname <- 'MAPK3'
  } else if (protein == 'ERK2') {
    pname <- 'MAPK1'
  }  else if (protein == 'MEK1') {
    pname <- 'MAP2K1'
  }  else if (protein == 'MEK2') {
    pname <- 'MAP2K2'
  }

  pApms <- genOrgenellPlot(expDesign, apmsPreys, "AP-MS", protein) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  pTid <- genOrgenellPlot(expDesign, tidPreys, "TbID", protein) +
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  
  p <- plot_grid(pApms, pTid, nrow = 1)
  # ggsave(
  #   plot = p,
  #   filename = file.path(resultsDir, 
  #                        paste0('subcellular_localization/organelleip_umaps/',
  #                               protein, '.pdf') ),
  #   width = 8, height = 4
  # )
  ggsave(
    plot = p,
    filename = file.path(figureDir(), 
                         paste0('S4/',
                                protein, '.png') ),
    width = 7.5, height = 3, dpi = 300
  )
  
  pList[[protein]] <- p

}
### legend
library(dplyr)

df_labels <- ips[ips$graph_localization_annotation!='unclassified',] %>%
  group_by(graph_localization_annotation) %>%
  summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2))

umapColors <- c(
  cytosol='#00a69c', proteasome='#c90d7c', '14-3-3_scaffold'='#58595b',
  actin_cytoskeleton='#c76893', plasma_membrane='#708097', recycling_endosome='#b81c8c',
  early_endosome='#8bc53f', "trans-Golgi"='#404041', Golgi='#26a9e0',ERGIC='#ffdd15',
  ER='#603813', peroxisome='#eb008b', lysosome='#ec1c24', mitochondrion='#ffb216',
  nucleus='#fde2b9', 'p-body'='#404041', stress_granule = '#72be44',translation='#dc4297',
  nucleolus='#7f3f97', centrosome='#0071bb')

ggplot(ips[ips$graph_localization_annotation != 'unclassified',], 
       aes(x=umap_1, y=umap_2, color=graph_localization_annotation)) +
  geom_point(  size=2, alpha=0.25, pch=16, show.legend = F ) +
  scale_color_manual(values = umapColors) +
  theme_cowplot() +
  geom_text_repel(data = df_labels, aes(umap_1, 
                                        umap_2, 
                                        label = graph_localization_annotation)) + 
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())

ggsave(
  filename = file.path(figureDir(), 
                       'S4/umap_legend.pdf' ),
  width = 4, height = 4
)
