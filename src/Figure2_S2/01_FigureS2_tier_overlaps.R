#' @author Sebastian Didusch
#' @date 25.03.24
#' 
#' Visualization script for venn diagrams/euler plots
#' first merged AP-MS + TurboID in WT for pan-RAS, pan-RAF, pan-MEK, pan-ERKs
#' then RAS-RAF cross layer paralogs, RAF-MEK, MEK-ERK
#' then activating mutants 

source(here::here("src/utils/utils_path.R"))
library(eulerr)

mergedPreyList <- list()

prepOverlapMat <- function(filtApmsList) {
  apmsNodePreyList <- list(
    ERK = unique(unlist(filtApmsList[grep('erk', names(filtApmsList), ignore.case = T )])),
    MEK = unique(unlist(filtApmsList[grep('mek', names(filtApmsList), ignore.case = T )])),
    RAF = unique(unlist(filtApmsList[grep('raf', names(filtApmsList), ignore.case = T )])),
    RAS = unique(unlist(filtApmsList[grep('ras', names(filtApmsList), ignore.case = T )]))
  )
  
  preys <- sort(unique(unlist(apmsNodePreyList)))
  
  mat <- t(+sapply(apmsNodePreyList, "%in%", x = preys))  ## matrix output
  colnames(mat) <- preys
  
  mat <- t(mat)
  mat <- as.data.frame(mat)
  mat$Gene.names <- rownames(mat)
  mat
}

overlapBinaryMats <- apmsBinaryMat <- tidBinaryMat <- list()
corePreys <- c('KRAS', 'NRAS', 'HRAS', 
               'ARAF', 'BRAF', 'RAF1',
               'MAP2K1', 'MAP2K2', 'MAPK3', 'MAPK1'
)

for (type in c("WT", "Active", "Inactive" ) ) {
  mergedPreyList <- list()
  bois <- apmsDesign$groups[apmsDesign$type==type & apmsDesign$batch==1]
  
  for (bait in bois ) {
    mergedPreyList[[bait]] <- unique(c(apmsPreys[[bait]], tidPreys[[bait]]))
  }
  
  nodePreyList <- list(
    ERK = unique(unlist(mergedPreyList[grep('erk', names(mergedPreyList), ignore.case = T )])),
    MEK = unique(unlist(mergedPreyList[grep('mek', names(mergedPreyList), ignore.case = T )])),
    RAF = unique(unlist(mergedPreyList[grep('raf', names(mergedPreyList), ignore.case = T )])),
    RAS = unique(unlist(mergedPreyList[grep('ras', names(mergedPreyList), ignore.case = T )]))
  )
  
  preys <- sort(unique(unlist(nodePreyList)))
  
  mat <- t(+sapply(nodePreyList, "%in%", x = preys))  ## matrix output
  colnames(mat) <- preys
  
  mat <- t(mat)
  mat <- as.data.frame(mat)
  mat$Gene.names <- rownames(mat)
  overlapBinaryMats[[type]] <- mat
  
  
  binMat <- UpSetR::fromList(nodePreyList)
  binMat <- binMat[, c('RAS', 'RAF', 'MEK', 'ERK')]
  
  # pdf(file.path(resultsDir, 'node_paralog_specificity/overlap/',
  #               paste0('upset_node_', type, '.pdf') ),
  #     width = 5, height = 4, onefile = F)
  # print(UpSetR::upset(binMat,
  #                     sets = c("ERK", "MEK", "RAF","RAS"),
  #                     keep.order = T,
  #                     intersections = list(list('ERK'), list('MEK'),
  #                                          list('RAF'), list('RAS'),
  #                                          list("RAS", "RAF", "MEK"),
  #                                          list("RAS", "RAF"),
  #                                          list("RAF", "MEK"),
  #                                          list("MEK", "ERK")
  #                                          ),
  #                     text.scale = c(1.5, 1.5, 1.5, 1, 1.5, 1.5),
  #                     mainbar.y.label = "# of overlapping preys",  # Change "Intersection size"
  #                     mainbar.y.max = 600,
  #                     sets.x.label = "# of preys"
  # )
  # )
  # dev.off()
  
  # eulerr
  # pdf(file.path(resultsDir, 'node_paralog_specificity/euler/node_overlap', paste0('euler_node_', type, '.pdf') ),
  #     width = 6, height = 4, onefile = F)
  # print(plot(
  #   euler(nodePreyList),
  #   legend = T,
  #   fills = rev(annotation_colors$node),
  #   quantities = T
  # ))
  # dev.off()
  filtApmsList <- apmsPreys[bois]
  filtTidList <- tidPreys[bois]
  #
  apmsMat <- prepOverlapMat(filtApmsList )
  tidMat <- prepOverlapMat(filtTidList )
  
  overlapBinaryMats[[type]] <- mat
  apmsBinaryMat[[type]] <- apmsMat
  tidBinaryMat[[type]] <- tidMat
}

# dir.create(file.path(figureDir(), 'Figure2/data/overlap/'), showWarnings = F)
write.table(
  overlapBinaryMats$WT,
  file.path(figureDir(), 'S2/data/combined_wt_interactome_binary_matrix.tsv'),
  sep = '\t', quote = F
)
write.table(
  overlapBinaryMats$Inactive,
  file.path(figureDir(), 'S2/data/combined_inactive_interactome_binary_matrix.tsv'),
  sep = '\t', quote = F
)
write.table(
  overlapBinaryMats$Active,
  file.path(figureDir(), 'S2/data/combined_active_interactome_binary_matrix.tsv'),
  sep = '\t', quote = F
)

prepTierOverlapDf <- function(toIt) {
  nodes <- c('RAS', 'RAF', 'MEK', 'ERK')
  toCheck <- list(
    RASRAFMEK=c('RAS', 'RAF', 'MEK'),
    RASRAF=c('RAS', 'RAF'),
    RAFMEK=c('RAF', 'MEK'),
    MEKERK=c('MEK', 'ERK') )
  
  allDfs <- list()
  for (elem in names(toIt)) {
    out <- list()
    mat <- toIt[[elem]]
    
    for (tier in nodes) {
      others <- setdiff(nodes, tier)
      n <- nrow(mat[apply(mat[, others], 1, sum) == 0,] )
      out[[paste0(elem, '_', tier) ]] <- data.frame(tier=tier, prey='tier-specific', n=n, type=elem)
    }
    for (cat in names(toCheck) ) {
      tiers <- toCheck[[cat]]
      others <- setdiff(nodes, tiers )
      tmp <- mat[apply(mat[, others, drop=F], 1, sum) == 0,]
      n <- nrow(tmp[apply(tmp[, tiers], 1, sum) == length(tiers), ] )
      
      for (tier in tiers) {
        out[[paste0(elem, '_', cat, '_', tier) ]] <- data.frame(tier=tier, prey=cat, n=n, type=elem)
      }
    }
    tmpDf <- do.call(rbind, out )
    for (tier in nodes) {
      alloc <- sum(tmpDf$n[tmpDf$tier==tier])
      nOthers <- nrow(mat[mat[[tier]] == 1, ]) - alloc
      tmpDf <- rbind(tmpDf, 
                     data.frame(tier=tier, prey='other', n=nOthers, type=elem)
      )
    }
    allDfs[[elem]] <- tmpDf 
    
  }
  df <- do.call(rbind, allDfs)
  df$tier <- factor(df$tier, levels = nodes)
  df$prey <- factor(df$prey, levels = rev(c(
    "tier-specific", "RASRAFMEK", "RASRAF", "RAFMEK", "MEKERK","other"    
  )))
  df$type <- factor(df$type, levels = c('Inactive', 'WT', 'Active'))
  
  df$Prey <- ifelse(
    df$prey == 'tier-specific',
    'tier-specific',
    'inter-tier'
  )
  
  df$percent <- 0
  for (type in names(toIt)) {
    for (tier in nodes) {
      rnames <- rownames(df[df$type == type &
                              df$tier == tier, ])
      df[rnames,]$percent <- round(df[rnames,]$n/sum(df[rnames,]$n )* 100, 2)
    }
  }
  df
}


#####
wtMat <- overlapBinaryMats$WT
inactivMat <- overlapBinaryMats$Inactive
activMat <- overlapBinaryMats$Active

toIt <- list(WT=wtMat, Inactive=inactivMat, Active=activMat)
df <- prepTierOverlapDf(toIt)

ylim_max <- max(df$n)  # or, for stacked bars:
ylim_max <- max(df %>% group_by(tier, type) %>% summarise(total=sum(n)) %>% pull(total))

df$type <- gsub('Inactive', 'CI', df$type)
df$type <- gsub('Active', 'CA', df$type)
df$type <- factor(df$type, levels = c('CI', 'WT', 'CA' ))

pCmb <- ggplot(df, aes(x=tier, y=n, fill=prey, label=percent )) + 
  geom_bar(stat = 'identity', position='stack') +
  scale_fill_manual(values = c('tier-specific'='#7fcdbb',
                               RASRAFMEK='#084594',
                               RASRAF='#4292c6',
                               RAFMEK='#9ecae1',
                               MEKERK='#deebf7',
                               other='grey80'
                               )) +
  ylab('number of preys') + xlab('') +
  theme_cowplot() + theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ), 
    legend.position = 'bottom'
    ) + facet_wrap(~type)
#
ggsave(plot = pCmb,
  file.path(figureDir(), 'S2/intra_inter_tier_prey_distribution.pdf'),
  width = 4.85, height = 2.9
)

write.table(
  df, 
  file.path(figureDir(), 'S2/data/intra_inter_tier_prey_distribution.tsv'),
  row.names = F, sep = '\t', quote = F
)


##### AP-MS overlap
toIt <- list(WT=apmsBinaryMat$WT, Inactive=apmsBinaryMat$Inactive, Active=apmsBinaryMat$Active)
df <- prepTierOverlapDf(toIt)

df$method <- 'AP-MS'
#
toIt <- list(WT=tidBinaryMat$WT, Inactive=tidBinaryMat$Inactive, Active=tidBinaryMat$Active)
dfTbid <- prepTierOverlapDf(toIt)
dfTbid$method <- 'TbID'

df <- rbind(df, dfTbid)

df$method <- factor(df$method, levels = c('AP-MS', 'TbID', 'Combined'))
#df$type <- factor('Inactive', 'WT', 'Active')
df$type <- gsub('Inactive', 'CI', df$type)
df$type <- gsub('Active', 'CA', df$type)
df$type <- factor(df$type, levels = c('CI', 'WT', 'CA' ))

pMethods <- ggplot(df, aes(x=tier, y=n, fill=prey, label=percent )) + 
  geom_bar(stat = 'identity', position='stack') +
  scale_fill_manual(values = c('tier-specific'='#7fcdbb',
                               RASRAFMEK='#084594',
                               RASRAF='#4292c6',
                               RAFMEK='#9ecae1',
                               MEKERK='#deebf7',
                               other='grey80'
  )) +
  ylab('number of preys') + xlab('') +
  theme_cowplot() + theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )) + facet_wrap(~type + method, nrow=1)
#
ggsave(
  file.path(figureDir(), 'S2/methods_intra_inter_tier_prey_distribution.pdf'),
  width = 6.5, height = 5
)

### same scale and legend
pCmb <- pCmb + ylim(0, ylim_max)
pMethods <- pMethods + ylim(0, ylim_max)



# Remove legends from individual plots
pCmb_no_leg <- pCmb + theme(legend.position = "none")
pMethods_no_leg <- pMethods + theme(legend.position = "none")

# Combine vertically or horizontally
combined_plot <- plot_grid(
  pCmb_no_leg, pMethods_no_leg,
  ncol = 2, align = "h"
)

# Add legend at the bottom
final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave(
  plot = final_plot, 
  file.path(figureDir(), 'S2/FigS2_A_B.pdf'),
  width = 10, height = 4.5
)
####
toIt <- list(WT=overlapBinaryMats$WT, Inactive=overlapBinaryMats$Inactive, Active=overlapBinaryMats$Active)
dfCmbd <- prepTierOverlapDf(toIt)
dfCmbd$method <- 'Combined'

df <- rbind(df, dfCmbd)
write.table(
  df, 
  file.path(figureDir(), 'S2/data/methods_intra_inter_tier_prey_distribution.tsv'),
  row.names = F, sep = '\t', quote = F
)
