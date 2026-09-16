#' @author Sebastian Didusch
#' @date 25.09.02
#' core pathway wiring
#' how are the core pathways connected to each other?
#' edge widths are proportional to log2 FCs
#' 
#' How do other, published and comprehensive studies compare?
#' (Buljan et al. and Kennedy et al. both used the 7 kinases from the
#' RAF-MEK-ERK tiers as bait proteins)
#' 

source(here::here("src/utils/utils_path.R"))

library(igraph)
library(viridis)

prepareWiringPlot <- function(G, minApms, maxApms) {
  fixed_positions <- matrix(NA, nrow=vcount(G), ncol=2)
  fixed_positions[V(G)$name == 'ARAF', ] <- c(-0.75, 1)     # ARAF
  fixed_positions[V(G)$name == 'BRAF', ] <- c(0, 1)   # BRAF
  fixed_positions[V(G)$name == 'KRAS', ] <- c(-1, 2)  # HRAS
  fixed_positions[V(G)$name == 'HRAS', ] <- c(1, 2)     # KRAS
  fixed_positions[V(G)$name == 'MEK1', ] <- c(-0.5, 0)     # MEK1
  fixed_positions[V(G)$name == 'MEK2', ] <- c(0.5, 0)   # MEK2
  fixed_positions[V(G)$name == 'ERK2', ] <- c(0.5, -1)     # ERK2
  fixed_positions[V(G)$name == 'ERK1', ] <- c(-0.5, -1)     # ERK1
  fixed_positions[V(G)$name == 'NRAS', ] <- c(0, 2.75)     # NRAS
  fixed_positions[V(G)$name == 'RAF1', ] <- c(0.75, 1)     # RAF1
  
  
  edge_values <- E(G)$log2FC
  normalized_values <- (edge_values - min(minApms)) / (max(maxApms) - min(minApms))
  
  #colorRampPalette(brew_colors(8, ''))
  edge_colors <- viridis(length(normalized_values))[as.numeric(cut(normalized_values, breaks=length(normalized_values)))]
  E(G)$color <- edge_colors
  
  V(G)$vertex.color <- annotation_colors$paralog[V(G)$name ]
  
  V(G)$id <- V(G)$label <- V(G)$name
  E(G)$curved <- 0.4
  
  V(G)$label.family <- 'Helvetica'
  
  labelPos <- rep(0, length(G))
  labelPos[V(G)$name == 'ARAF' ] <- -3     # ARAF
  labelPos[V(G)$name == 'BRAF' ] <- 2.5   # BRAF
  labelPos[V(G)$name == 'KRAS' ] <- 2.5  # HRAS
  labelPos[V(G)$name == 'HRAS' ] <- 2.5     # KRAS
  labelPos[V(G)$name == 'MEK1' ] <- -2.5     # MEK1
  labelPos[V(G)$name == 'MEK2' ] <- -2.5   # MEK2
  labelPos[V(G)$name == 'ERK2' ] <- -2.5    # ERK2
  labelPos[V(G)$name == 'ERK1' ] <- -2.5     # ERK1
  labelPos[V(G)$name == 'NRAS' ] <- 2.5   # NRAS
  labelPos[V(G)$name == 'RAF1' ] <- -3     # RAF1
  
  V(G)$label.position <- labelPos
  return(list(G=G, pos=fixed_positions))
}

### ERK substrates/known interactors wiring
prepareERKWiringPlot <- function(G, minApms, maxApms) {
  fixed_positions <- matrix(NA, nrow=vcount(G), ncol=2)
  fixed_positions[V(G)$name == 'ERK1', ] <- c(-0.5, 0)     # ERK1
  fixed_positions[V(G)$name == 'ERK2', ] <- c(0.5, 0)     # ERK2
  fixed_positions[V(G)$name == 'RPS6KA1', ] <- c(-1.5, -0.25)     # RPS6KA1
  fixed_positions[V(G)$name == 'RPS6KA2', ] <- c(-0.75, -0.33)     # RPS6KA2
  fixed_positions[V(G)$name == 'RPS6KA3', ] <- c(0, -0.4)     # RPS6KA3
  fixed_positions[V(G)$name == 'RPS6KA4', ] <- c(0.75, -0.33)     # RPS6KA4
  fixed_positions[V(G)$name == 'RPS6KA5', ] <- c(1.5, -0.25)     # RPS6KA5
  fixed_positions[V(G)$name == 'MKNK1', ] <- c(-0.5, -0.75)     # MKNK1
  fixed_positions[V(G)$name == 'MKNK2', ] <- c(0.5, -0.75)     # MKNK2
  #
  fixed_positions[V(G)$name == 'DUSP1', ] <- c(-1.5, -0.5)     # MKNK1
  fixed_positions[V(G)$name == 'DUSP5', ] <- c(-1.4, -0.6)     # MKNK1
  fixed_positions[V(G)$name == 'DUSP6', ] <- c(-1.125, -0.75)     # MKNK1
  fixed_positions[V(G)$name == 'DUSP7', ] <- c(-0.75, -1)     # MKNK1
  fixed_positions[V(G)$name == 'DUSP9', ] <- c(0, -1)     # MKNK2
  fixed_positions[V(G)$name == 'DUSP16', ] <- c(0.75, -1)     # MKNK1
  
  fixed_positions[V(G)$name == 'EIF4EBP1', ] <- c(1.5, -0.5)     # MKNK1
  fixed_positions[V(G)$name == 'TCF3', ] <- c(1.125, -0.75)     # MKNK1
  fixed_positions[V(G)$name == 'TOB1', ] <- c(1, -0.4)     # MKNK1
  
  fixed_positions[V(G)$name == 'MEK2', ] <- c(0.5, 0.5)     # MEK2
  fixed_positions[V(G)$name == 'SGK1', ] <- c(0.25, 0.5)     # MEK2
  #fixed_positions[V(G)$name == 'ILK', ] <- c(-0.25, 0.75)     # ILK
  #fixed_positions[V(G)$name == 'LIMS1', ] <- c(-0.5, 0.75)     # LIMS1
  fixed_positions[V(G)$name == 'PARVA', ] <- c(-0.75, 0.75)     # PARVA
  #fixed_positions[V(G)$name == 'PARVB', ] <- c(-0.75, 1)     # PARVB
  #fixed_positions[V(G)$name == 'RSU1', ] <- c(-1.25, 0.5)     # RSU1
  fixed_positions[V(G)$name == 'FGFR1', ] <- c(1.25, 0.75)     # FGFR1
  fixed_positions[V(G)$name == 'THRB', ] <- c(0, 0.75)     # FGFR1
  fixed_positions[V(G)$name == 'ADAM17', ] <- c(0.5, 0.75)     # FGFR1
  fixed_positions[V(G)$name == 'GAB1', ] <- c(1, 0.75)     # GAB1
  fixed_positions[V(G)$name == 'PEA15', ] <- c(1.25, 0.25)     # PEA15
  
  fixed_positions[V(G)$name == 'ERF', ] <- c(-1.25, 0.25)     # RSU1
  fixed_positions[V(G)$name == 'ETV3', ] <- c(-0.75, 0.25)     # RSU1
  fixed_positions[V(G)$name == 'MRTFA', ] <- c(-0.25, 0.25)     # RSU1
  
  edge_values <- E(G)$log2FC
  normalized_values <- (edge_values - min(minApms)) / (max(maxApms) - min(minApms))
  
  #colorRampPalette(brew_colors(8, ''))
  edge_colors <- viridis(length(normalized_values))[as.numeric(cut(normalized_values, breaks=length(normalized_values)))]
  E(G)$color <- edge_colors
  
  V(G)$vertex.color <- 'grey60'
  V(G)$vertex.color[V(G)$name == 'ERK1'] <- '#a62877'
  V(G)$vertex.color[V(G)$name == 'ERK2'] <- '#dfa0c2'
  # annotation_colors$paralog[V(G)$name ]
  
  V(G)$id <- V(G)$label <- V(G)$name
  E(G)$curved <- 0.4
  
  V(G)$label.family <- 'Helvetica'
  
  labelPos <- rep(-1.75, length(G))
  
  labelDegree <- rep(3*pi/2, length(G))
  
  #labelDegree[V(G)$name == 'RSU1' ] <- pi/2 
  labelDegree[V(G)$name == 'PARVB' ] <- pi/2 
  labelDegree[V(G)$name == 'PARVA' ] <- 0
  labelDegree[V(G)$name == 'ILK' ] <- pi
  labelDegree[V(G)$name == 'LIMS1' ] <- pi/2 
  labelDegree[V(G)$name == 'FGFR1' ] <- pi/2 
  labelDegree[V(G)$name == 'THRB' ] <- 0 
  labelDegree[V(G)$name == 'ADAM17' ] <- pi/2 
  labelDegree[V(G)$name == 'RPS6KA1' ] <- pi/2
  labelDegree[V(G)$name == 'RPS6KA5' ] <- pi/2
  labelDegree[V(G)$name == 'DUSP1' ] <- 0
  labelDegree[V(G)$name == 'DUSP5' ] <- 0
  labelDegree[V(G)$name == 'DUSP6' ] <- 0
  labelDegree[V(G)$name == 'DUSP7' ] <- 0
  labelDegree[V(G)$name == 'DUSP16' ] <- pi
  labelDegree[V(G)$name == 'PEA15' ] <- pi
  labelDegree[V(G)$name == 'TCF3' ] <- pi
  labelDegree[V(G)$name == 'ERK1' ] <- pi/2
  labelDegree[V(G)$name == 'ERK2' ] <- pi/2 
  labelDegree[V(G)$name == 'MEK2' ] <- pi
  labelDegree[V(G)$name == 'SGK1' ] <- pi/2 
  labelDegree[V(G)$name == 'GAB1' ] <- pi
  labelDegree[V(G)$name == 'RPS6KA2' ] <- 0
  labelDegree[V(G)$name == 'RPS6KA4' ] <- pi
  labelDegree[V(G)$name == 'ERF' ] <- 0
  labelDegree[V(G)$name == 'MRTFA' ] <- pi/2
  labelDegree[V(G)$name == 'ETV3' ] <- pi/2
  
  V(G)$label.position <- labelPos
  V(G)$label.degree <- labelDegree
  return(list(G=G, pos=fixed_positions))
}

coreNames <- c(KRAS='KRAS', NRAS='NRAS', HRAS='HRAS',
               ARAF='ARAF', BRAF='BRAF', RAF1='RAF1',
               MEK1='MAP2K1', MEK2='MAP2K2',
               ERK1='MAPK3', ERK2='MAPK1'
)

coreInteractome <- interactome[interactome$prey %in% coreNames, ]
coreInteractome$bait_name <- apmsDesign$paralog[match(coreInteractome$bait, apmsDesign$groups) ]
coreInteractome$prey_name <- coreInteractome$prey

coreInteractome$prey_name[coreInteractome$prey_name=='MAP2K1'] <- 'MEK1'
coreInteractome$prey_name[coreInteractome$prey_name=='MAP2K2'] <- 'MEK2'
coreInteractome$prey_name[coreInteractome$prey_name=='MAPK3'] <- 'ERK1'
coreInteractome$prey_name[coreInteractome$prey_name=='MAPK1'] <- 'ERK2'

coreInteractome$in_apms <- 'no'
coreInteractome$in_turboid <- 'no'

for (bait in unique(coreInteractome$bait )) {
  coreInteractome$in_apms[coreInteractome$bait == bait &
                            coreInteractome$prey %in% apmsPreys[[bait]] ] <- 'yes'
  #
  coreInteractome$in_turboid[coreInteractome$bait == bait &
                               coreInteractome$prey %in% tidPreys[[bait]] ] <- 'yes'
}

coreInteractome <- coreInteractome[coreInteractome$bait_name != coreInteractome$prey_name, ]

minApms <- min(coreInteractome$logFC_APMS[coreInteractome$in_apms == 'yes'])
maxApms <- max(coreInteractome$logFC_APMS[coreInteractome$in_apms == 'yes'])
minTid <- min(coreInteractome$logFC_TbID[coreInteractome$in_turboid == 'yes'])
maxTid <- max(coreInteractome$logFC_TbID[coreInteractome$in_turboid == 'yes'])

minApms <- minTid <- min(minApms, minTid)
maxApms <- maxTid <- max(maxApms, maxTid)

type <- 'WT'

#dev.off()
apmsNetworks <- list()
tidNetworks <- list()

op <- file.path(figureDir(), 'Figure1/core_pathway_interactions/')
dir.create(op, showWarnings = FALSE) 
#'/home/diduschs92/Documents/proj/ErkPathwayReferenceMap/results/rewiring/igraph/'

for (type in c('Inactive', 'WT', 'Active') ) {
  ##
  nw <- type
  ##
  stateCoreInteractome <- 
    coreInteractome[coreInteractome$bait%in%apmsDesign$groups[apmsDesign$type==type], ]
  
  apmsDat <- stateCoreInteractome[stateCoreInteractome$in_apms=='yes', ]
  tidDat <- stateCoreInteractome[stateCoreInteractome$in_turboid=='yes', ]
  
  apmsDat <- apmsDat[order(apmsDat$logFC_APMS, decreasing = T), ]
  tidDat <- tidDat[order(tidDat$logFC_TbID, decreasing = T), ]
  
  apmsDat$dupl <- paste0(apmsDat$bait_name, '_', apmsDat$prey_name)
  tidDat$dupl <- paste0(tidDat$bait_name, '_', tidDat$prey_name)
  
  apmsDat <- apmsDat[!duplicated(apmsDat$dupl),]
  tidDat <- tidDat[!duplicated(tidDat$dupl),]
  
  apmsDat <- apmsDat[, c('bait_name', 'prey_name', 'logFC_APMS')]
  apmsDat <- apmsDat[order(apmsDat$bait_name),]
  #
  tidDat <- tidDat[, c('bait_name', 'prey_name', 'logFC_TbID')]
  tidDat <- tidDat[order(tidDat$bait_name),]
  #
  names(apmsDat) <- names(tidDat) <- c('bait_name', 'prey_name', 'log2FC')
  ##
  
  G_apms <- graph_from_data_frame(apmsDat)
  apmsOut <- prepareWiringPlot(G_apms, minApms, maxApms)
  
  G_apms <- apmsOut$G
  fixed_positions <- apmsOut$pos
  
  apmsNetworks[[nw]] <- G_apms
  
  pdf(paste0(op, 'apms_', type, '_core.pdf'), width = 4.9, height = 5.41)
  plot(G_apms,
       edge.curved=E(G_apms)$curved,
       vertex.size=10,
       vertex.color= V(G_apms)$vertex.color,
       vertex.frame.width=0,
       vertex.label.cex=1.5,
       edge.arrow.size=1,
       edge.width=E(G_apms)$log2FC,
       layout=fixed_positions,
       vertex.label.dist=V(G_apms)$label.position,
       vertex.label.color='black',
       edge.color='grey60' #E(G_apms)$color
       #main=nw
  )
  dev.off()
  legend_values <- seq(min(minApms), max(maxApms), length.out=5)
  legend_colors <- viridis(length(legend_values))
  legend_widths <- 1 + 5 * (legend_values - min(minApms)) / (max(maxApms) - min(minApms))
  
  pdf(paste0(op, 'apms_legend.pdf'), width = 4.9, height = 5.41)
  par(mar=c(4,1,1,1)) 
  plot.new() 
  legend("center", legend=round(legend_values, 2), col='grey60', lwd=legend_widths, title="log2FC", bty="n")
  dev.off()
  ### tbid
  
  G_tid <- graph_from_data_frame(tidDat)
  tidOut <- prepareWiringPlot(G_tid, minTid, maxTid)
  
  G_tid <- tidOut$G
  fixed_positions <- tidOut$pos
  
  tidNetworks[[nw]] <- G_tid
  # E(G_tid)$curved <- 0.4
  pdf(paste0(op, 'tid_', type, '_core.pdf'), width = 4.9, height = 5.41)
  plot(G_tid,
       edge.curved=E(G_tid)$curved,
       vertex.size=10,
       vertex.color=V(G_tid)$vertex.color,
       vertex.frame.width=0,
       vertex.label.cex=1.5,
       edge.arrow.size=1,
       edge.width=E(G_tid)$log2FC,
       layout=fixed_positions,
       vertex.label.dist=V(G_tid)$label.position,
       vertex.label.color='black',
       edge.color='grey60' #E(G_apms)$color
       #main=nw
  )
  dev.off()
  legend_values <- seq(min(minTid), max(maxTid), length.out=5)
  legend_colors <- viridis(length(legend_values))
  legend_widths <- 1 + 5 * (legend_values - min(minTid)) / (max(maxTid) - min(minTid))
  
  # Add the legend for edge colors and widths
  # legend("topleft", legend=round(legend_values, 2), col=legend_colors, lwd=legend_widths, title="log2FC")
  pdf(paste0(op, 'tid_legend.pdf'), width = 4.9, height = 5.41)
  par(mar=c(4,1,1,1))  # Adjust margins for the legend
  plot.new()  # Create a new plot
  legend("center", legend=round(legend_values, 2),
         col='grey60', lwd=legend_widths, title="log2FC", bty="n")
  dev.off()
}

# sapply(apmsNetworks, ecount)
# sapply(tidNetworks, ecount)

# sapply(apmsNetworks, function(x) sum(E(x)$log2FC) )
# sapply(tidNetworks, function(x) sum(E(x)$log2FC) )

############# ERK wiring
erkPreyNames <- c(RPS6KA1='RPS6KA1', RPS6KA2='RPS6KA2', RPS6KA3='RPS6KA3',
                  RPS6KA4='RPS6KA4', RPS6KA5='RPS6KA5', MKNK1='MKNK1',
                  MKNK2='MKNK2', ERK1='MAPK3', ERK2='MAPK1', DUSP7='DUSP7',
                  DUSP9='DUSP9', DUSP16='DUSP16', MEK2='MAP2K2', PARVA='PARVA',
                  FGFR1='FGFR1', GAB1='GAB1', PEA15='PEA15', THRB='THRB', 
                  ADAM17='ADAM17', SGK1='SGK1', DUSP1='DUSP1', DUSP5='DUSP5',
                  TOB1='TOB1', TCF3='TCF3',EIF4EBP1='EIF4EBP1',
                  DUSP6='DUSP6', ERF='ERF', ETV3='ETV3', MRTFA='MRTFA'
)

erkInteractome <- interactome[interactome$prey %in% erkPreyNames, ]
erkInteractome <- erkInteractome[grep('ERK', erkInteractome$bait),]
erkInteractome$bait_name <- apmsDesign$paralog[match(erkInteractome$bait, apmsDesign$groups) ]
erkInteractome$prey_name <- erkInteractome$prey

erkInteractome$prey_name[erkInteractome$prey_name=='MAPK3'] <- 'ERK1'
erkInteractome$prey_name[erkInteractome$prey_name=='MAPK1'] <- 'ERK2'
erkInteractome$prey_name[erkInteractome$prey_name=='MAP2K2'] <- 'MEK2'

erkInteractome$in_apms <- 'no'
erkInteractome$in_turboid <- 'no'

for (bait in unique(erkInteractome$bait )) {
  erkInteractome$in_apms[erkInteractome$bait == bait &
                            erkInteractome$prey %in% apmsPreys[[bait]] ] <- 'yes'
  #
  erkInteractome$in_turboid[erkInteractome$bait == bait &
                               erkInteractome$prey %in% tidPreys[[bait]] ] <- 'yes'
}

### erk targets in signor
comparison <- read.table(
  file.path(figureDir(), 'Figure1/data/euler_erkpreys_intact_erk_as_baits.tsv'),
  sep = '\t', header = T)
targets <- comparison$ERK_prey_protein[comparison$ERK.relationship.in.Signor==1 &
                                         comparison$in.IntAct==0]
targets <- c(targets, c('MEK1', 'MEK2') )

erkInteractome <- erkInteractome[erkInteractome$bait_name != erkInteractome$prey_name, ]
for (method in c('in_apms', 'in_turboid' )) {
  erkNet <- erkInteractome[erkInteractome[[method]] =='yes', ] 
  erkNet <- erkNet[!duplicated(paste0(erkNet$bait_name, '_',
                                      erkNet$prey_name)), ]
  #
  erkNet <- erkNet[, c('bait_name', 'prey_name', method)]
  erkNet$edge <- ifelse(erkNet$prey_name %in% targets,
                        'SIGNOR', 'IntAct & SIGNOR')
  write.table(erkNet, file.path(figureDir(), paste0('S1/data/erktargets_',method,'.tsv' )),
              row.names = F, sep = '\t', quote = F)
  
}
erkNet <- erkInteractome[!duplicated(paste0(erkInteractome$bait_name, '_',
                                            erkInteractome$prey_name)), ]

minApms <- min(erkInteractome$logFC_APMS[erkInteractome$in_apms == 'yes'])
maxApms <- max(erkInteractome$logFC_APMS[erkInteractome$in_apms == 'yes'])
minTid <- min(erkInteractome$logFC_TbID[erkInteractome$in_turboid == 'yes'])
maxTid <- max(erkInteractome$logFC_TbID[erkInteractome$in_turboid == 'yes'])

minApms <- minTid <- min(minApms, minTid)
maxApms <- maxTid <- max(maxApms, maxTid)

type <- 'WT'

#dev.off()
apmsNetworks <- list()
tidNetworks <- list()
op <- file.path(figureDir(), 'Figure1/core_pathway_interactions/ERK/')
dir.create(op, showWarnings = FALSE) 

for (type in c('Inactive', 'WT', 'Active') ) {
  ##
  nw <- type
  ##
  stateerkInteractome <- 
    erkInteractome[erkInteractome$bait%in%apmsDesign$groups[apmsDesign$type==type], ]
  
  apmsDat <- stateerkInteractome[stateerkInteractome$in_apms=='yes', ]
  tidDat <- stateerkInteractome[stateerkInteractome$in_turboid=='yes', ]
  
  apmsDat <- apmsDat[order(apmsDat$logFC_APMS, decreasing = T), ]
  tidDat <- tidDat[order(tidDat$logFC_TbID, decreasing = T), ]
  
  apmsDat$dupl <- paste0(apmsDat$bait_name, '_', apmsDat$prey_name)
  tidDat$dupl <- paste0(tidDat$bait_name, '_', tidDat$prey_name)
  
  apmsDat <- apmsDat[!duplicated(apmsDat$dupl),]
  tidDat <- tidDat[!duplicated(tidDat$dupl),]
  
  apmsDat <- apmsDat[, c('bait_name', 'prey_name', 'logFC_APMS')]
  apmsDat <- apmsDat[order(apmsDat$bait_name),]
  #
  tidDat <- tidDat[, c('bait_name', 'prey_name', 'logFC_TbID')]
  tidDat <- tidDat[order(tidDat$bait_name),]
  #
  names(apmsDat) <- names(tidDat) <- c('bait_name', 'prey_name', 'log2FC')
  ##
  
  G_apms <- graph_from_data_frame(apmsDat)
  apmsOut <- prepareERKWiringPlot(G_apms, minApms, maxApms)
  
  G_apms <- apmsOut$G
  fixed_positions <- apmsOut$pos
  
  apmsNetworks[[nw]] <- G_apms
  
  E(G_apms)$color <- 'grey60'
  edgesToColor <- which(
    (tail_of(G_apms, E(G_apms))$name %in% c("ERK1", "ERK2")) & 
      (head_of(G_apms, E(G_apms))$name %in% targets)
  )
  if (length(edgesToColor) > 0) E(G_apms)$color[edgesToColor] <- '#41AB5D'
  
  pdf(paste0(op, 'apms_', type, '_erk.pdf'), width = 4.9, height = 5.41)
  plot(G_apms,
       edge.curved=E(G_apms)$curved,
       vertex.size=10,
       vertex.color= V(G_apms)$vertex.color,
       vertex.frame.width=0,
       vertex.label.cex=1.5,
       edge.arrow.size=1,
       edge.width=E(G_apms)$log2FC,
       layout=fixed_positions,
       vertex.label.dist=V(G_apms)$label.position,
       vertex.label.degree=V(G_apms)$label.degree,
       vertex.label.color='black',
       edge.color=E(G_apms)$color
       #main=nw
  )
  dev.off()
  legend_values <- seq(min(minApms), max(maxApms), length.out=5)
  legend_colors <- viridis(length(legend_values))
  legend_widths <- 1 + 5 * (legend_values - min(minApms)) / (max(maxApms) - min(minApms))
  
  pdf(paste0(op, 'apms_legend.pdf'), width = 4.9, height = 5.41)
  par(mar=c(4,1,1,1))  # Adjust margins for the legend
  plot.new()  # Create a new plot
  legend("center", legend=round(legend_values, 2), col='grey60', lwd=legend_widths, title="log2FC", bty="n")
  dev.off()
  ### tbid
  
  G_tid <- graph_from_data_frame(tidDat)
  tidOut <- prepareERKWiringPlot(G_tid, minTid, maxTid)
  
  G_tid <- tidOut$G
  fixed_positions <- tidOut$pos
  
  E(G_tid)$color <- 'grey60'
  edgesToColor <- which(
    (tail_of(G_tid, E(G_tid))$name %in% c("ERK1", "ERK2")) & 
      (head_of(G_tid, E(G_tid))$name %in% targets)
  )
  if (length(edgesToColor) > 0) E(G_tid)$color[edgesToColor] <- '#41AB5D'
  
  tidNetworks[[nw]] <- G_tid
  # E(G_tid)$curved <- 0.4
  pdf(paste0(op, 'tid_', type, '_erk.pdf'), width = 4.9, height = 5.41)
  plot(G_tid,
       edge.curved=E(G_tid)$curved,
       vertex.size=10,
       vertex.color=V(G_tid)$vertex.color,
       vertex.frame.width=0,
       vertex.label.cex=1.5,
       edge.arrow.size=1,
       edge.width=E(G_tid)$log2FC,
       layout=fixed_positions,
       vertex.label.dist=V(G_tid)$label.position,
       vertex.label.color='black',
       edge.color= E(G_tid)$color
       #main=nw
  )
  dev.off()
  legend_values <- seq(min(minTid), max(maxTid), length.out=5)
  legend_colors <- viridis(length(legend_values))
  legend_widths <- 1 + 5 * (legend_values - min(minTid)) / (max(maxTid) - min(minTid))
  
  # Add the legend for edge colors and widths
  # legend("topleft", legend=round(legend_values, 2), col=legend_colors, lwd=legend_widths, title="log2FC")
  pdf(paste0(op, 'tid_legend.pdf'), width = 4.9, height = 5.41)
  par(mar=c(4,1,1,1))  # Adjust margins for the legend
  plot.new()  # Create a new plot
  legend("center", legend=round(legend_values, 2),
         col='grey60', lwd=legend_widths, title="log2FC", bty="n")
  dev.off()
}

data.frame(
  name   = names(apmsNetworks),
  vcount = sapply(apmsNetworks, vcount),
  ecount = sapply(apmsNetworks, ecount)
)

data.frame(
  name   = names(tidNetworks),
  vcount = sapply(tidNetworks, vcount),
  ecount = sapply(tidNetworks, ecount)
)

######
p <- (file.path(dataDir(), 'integration/apms_pdl_datasets') )

kolch <- read.table(file.path(p, 'kolch_egfr_interactome_mt_lo.tsv'),
                    header = T, sep = '\t'
)

kolch <- kolch[kolch$t.test...................................P.value<0.05 &
                 kolch$Significance.A < 0.05,]
kolch <- kolch[, c('Bait', 'Prey',  "Log2..SILAC.ratio.")]
names(kolch) <- c('bait_name', 'prey_name', 'log2FC')

coreKolch <- kolch[kolch$bait_name %in% coreNames, ]
coreKolch <- coreKolch[coreKolch$prey_name %in% coreNames, ]
coreKolch <- coreKolch[coreKolch$bait_name != coreKolch$prey_name,]

gstaiger <- read.table(file.path(p, 'Kinase_interactome.tsv'),
                       header = T, sep = '\t'
)

gstaiger <- gstaiger[, c('Bait_id', 'Protein_id', 'gfpratio')]
names(gstaiger) <- c('bait_name', 'prey_name', 'log2FC')
coreGstaiger <- gstaiger[gstaiger$bait_name%in%coreNames, ]
coreGstaiger <- coreGstaiger[coreGstaiger$prey_name%in%coreNames, ]
coreGstaiger <- coreGstaiger[coreGstaiger$bait_name != coreGstaiger$prey_name,]

for (i in seq_along(coreNames)){
  coreKolch$bait_name <-  gsub(coreNames[i], names(coreNames)[i], coreKolch$bait_name )
  coreKolch$prey_name <-  gsub(coreNames[i], names(coreNames)[i], coreKolch$prey_name )
  #
  coreGstaiger$bait_name <-  gsub(coreNames[i], names(coreNames)[i], coreGstaiger$bait_name )
  coreGstaiger$prey_name <-  gsub(coreNames[i], names(coreNames)[i], coreGstaiger$prey_name )
  
}

### 
G_kolch <- graph_from_data_frame(coreKolch)
kolchOut <- prepareWiringPlot(G_kolch, min(coreKolch$log2FC), max(coreKolch$log2FC) )

G_tid <- kolchOut$G
fixed_positions <- kolchOut$pos

#E(G_tid)$curved <- 0.4
pdf(paste0(op, 'kolch_mikras_lo_core.pdf'), width = 4.9, height = 5.41)
plot(G_tid,
     edge.curved=E(G_tid)$curved,
     vertex.size=10,
     vertex.color=V(G_tid)$vertex.color,
     vertex.frame.width=0,
     vertex.label.cex=1.5,
     edge.arrow.size=1,
     edge.width=3, #E(G_tid)$log2FC,
     layout=fixed_positions,
     vertex.label.dist=V(G_tid)$label.position,
     vertex.label.color='black',
     edge.color='grey60' #E(G_apms)$color
     #main=nw
)
dev.off()

###
coreGstaiger$log2FC <- log2(coreGstaiger$log2FC)
G_gstaiger <- graph_from_data_frame(coreGstaiger)
gstaigerOut <- prepareWiringPlot(G_gstaiger, min(coreGstaiger$log2FC), max(coreGstaiger$log2FC) )

G_tid <- gstaigerOut$G
fixed_positions <- gstaigerOut$pos

#E(G_tid)$curved <- 0.4
pdf(paste0(op, 'gstaiger_core.pdf'), width = 4.9, height = 5.41)
plot(G_tid,
     edge.curved=E(G_tid)$curved,
     vertex.size=10,
     vertex.color=V(G_tid)$vertex.color,
     vertex.frame.width=0,
     vertex.label.cex=1.5,
     edge.arrow.size=1,
     edge.width=3, #E(G_tid)$log2FC/4,
     layout=fixed_positions,
     vertex.label.dist=V(G_tid)$label.position,
     vertex.label.color='black',
     edge.color='grey60' #E(G_apms)$color
     #main=nw
)
dev.off()

wtCoreInteractome <- coreInteractome[grep('WT', coreInteractome$bait ), ]
wtCoreInteractome <- wtCoreInteractome[!duplicated(paste0(wtCoreInteractome$bait_name, '_',
                                                          wtCoreInteractome$prey_name) ), ]

wtCoreInteractome$log2FC <- apply(wtCoreInteractome[, c('logFC_APMS', 'logFC_TbID' )],1, max)
G_wt_collapsed <- graph_from_data_frame(wtCoreInteractome[, c('bait_name', 'prey_name', 'log2FC')])
wtcOut <- prepareWiringPlot(G_wt_collapsed, min(wtCoreInteractome$log2FC), max(wtCoreInteractome$log2FC) )

G_tid <- wtcOut$G
fixed_positions <- wtcOut$pos

#E(G_tid)$curved <- 0.4
pdf(paste0(op, 'wt_core_apms_tbid_collapsed.pdf'), width = 4.9, height = 5.41)
plot(G_tid,
     edge.curved=E(G_tid)$curved,
     vertex.size=10,
     vertex.color=V(G_tid)$vertex.color,
     vertex.frame.width=0,
     vertex.label.cex=1.5,
     edge.arrow.size=1,
     edge.width=3, #E(G_tid)$log2FC,
     layout=fixed_positions,
     vertex.label.dist=V(G_tid)$label.position,
     vertex.label.color='black',
     edge.color='grey60' #E(G_apms)$color
     #main=nw
)
dev.off()
