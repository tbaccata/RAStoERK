#' @author Sebastian Didusch
#' @date 25.03.13
#' 
#' Script for final thresholding of interactome using
#' 2 negative controls for both methods: GFP-high expression (GFP_8ng, n = 9) 
#' and GFP low expression (GFP_2ng, n = 9).
#' HCPPIs had to be significantly (padj < 0.05) enriched over both controls.
#' Final, tier-specific thresholds were applied over GFP high.
#' 
#' GFP low helped us in defining background contaminants (data/raw/bystanders.txt)
#' It was not needed again in dataset B.
#' 
#' For dataset B, preys needed to be sufficiently enriched over GFP high (n = 3).
#' Again, bona-fide known interactors were used to determined log2FC thresholds.
#' 
#' Finally, a list of bystander from the CRAPome, as well as method-specific 
#' background was excluded from the final interactome.
#' 

## define paths
library(here)
dataDir   <- function(...) here("data", ...)
processedDataDir   <- function(...) here("data", "processed", ...)
interimDataDir   <- function(...) here("data", "interim", ...)

# bystanders
bystanders <- readLines(file.path(dataDir(), 'raw/bystanders.txt' ) )
# thresholds
highThresholds <- read.table(file.path(dataDir(), 'raw/thresholds/highctrl_thresholds.txt' ),
                             header = T, sep = '\t')

apmsDesign <- expDesign <- read.table(file.path(processedDataDir(), 'apms/design_apms.txt' ),
                                      header = T, sep = '\t')
tidDesign <- read.table(file.path(processedDataDir(), 'turboid/design_turboid.txt' ),
                        header = T, sep = '\t')
tidDesign$samples <- as.character(tidDesign$samples)

daTid <- read.table(file.path(processedDataDir(), 'turboid/deqms_turboid.tsv' ),
                    header = T, sep = '\t')

daApms <- read.table(file.path(processedDataDir(), 'apms/deqms_apms.tsv' ),
                     header = T, sep = '\t')

### thresholding vs GFP high (GFP_8ng) and GFP low (GFP_2ng)
# use group comparison against both baits as filters (padj < 0.05 and log2FC > 0)
# ultimate, final thresholding only vs GFP high (padj < 0.05 and log2FC > tier_threshold)
# thresholds have been determined by inspecitng log2FCs of bona fide known interactors

tidPreyList <- list()
G_tid <- NULL
for (bait in unique(apmsDesign$groups[apmsDesign$type != 'Control'])) {
  fcCol <- paste0('logFC_', bait, '__vs__GFP8ng')
  padjCol <- paste0('adj.P.Val_', bait, '__vs__GFP8ng')
  fcCol2 <- paste0('logFC_', bait, '__vs__GFP2ng')
  padjCol2 <- paste0('adj.P.Val_', bait, '__vs__GFP2ng')
  filtDaTid <- daTid[daTid[[fcCol]] > 0 &
                       daTid[[padjCol]] < 0.05, c('Gene.names', fcCol, padjCol)]
  
  node <- unique(apmsDesign$node[apmsDesign$groups==bait])
  tidThresh <- highThresholds$log2fc[highThresholds$node==node & highThresholds$method=='TurboID']
  
  filtDaTid <- daTid[daTid[[fcCol]] > tidThresh &
                       daTid[[padjCol]] < 0.05 &
                       daTid[[fcCol2]] > 0 &
                       daTid[[padjCol2]] < 0.05, c('Gene.names', fcCol, padjCol)]
  
  filtDaTid <- filtDaTid[!filtDaTid$Gene.names %in% bystanders, ]
  
  df <- data.frame(bait=bait, prey=filtDaTid$Gene.names)
  if (is.null(G_tid)) {
    G_tid <- df
  } else {
    G_tid <- rbind(G_tid, df)
  }
  
  tidPreyList[[bait]] <- filtDaTid$Gene.names
}

apmsPreyList <- list()
G_apms <- NULL
for (bait in unique(apmsDesign$groups[apmsDesign$type != 'Control'])) {
  node <- unique(apmsDesign$node[apmsDesign$groups==bait])
  apmsThresh <- highThresholds$log2fc[highThresholds$node==node & highThresholds$method=='AP-MS']
  
  fcCol <- paste0('logFC_', bait, '__vs__GFP8ng')
  padjCol <- paste0('adj.P.Val_', bait, '__vs__GFP8ng')
  fcCol2 <- paste0('logFC_', bait, '__vs__GFP2ng')
  padjCol2 <- paste0('adj.P.Val_', bait, '__vs__GFP2ng')
  
  preys <- daApms$Gene.names[daApms[[fcCol]] > apmsThresh &
                               daApms[[padjCol]] < 0.05 &
                               daApms[[fcCol2]] > 0 &
                               daApms[[padjCol2]] < 0.05
  ]

  preys <- preys[!preys %in% bystanders]
  
  df <- data.frame(bait=bait, prey=preys)
  if (is.null(G_apms)) {
    G_apms <- df
  } else {
    G_apms <- rbind(G_apms, df)
  }
  apmsPreyList[[bait]] <- preys
}

apms_MPLID25011_DA <- read.table(file.path(dataDir(), 'processed/apms/deqms_apms_MPLID25011.tsv'),
                         header = T, sep = '\t'
)
turboid_MPLID25011_DA <- read.table(file.path(dataDir(), 'processed/turboid/deqms_turboid_MPLID25011.tsv'),
                        header = T, sep = '\t'
)

baits <- gsub('__vs__.*', '', grep('logFC', names(apms_MPLID25011_DA), value = T) )
baits <- gsub('logFC_', '', baits)

apms_MPLID25011_df <- list()
tid_MPLID25011_df <- list()

for (bait in baits ) {
  fcCol <- paste0("logFC_", bait, "__vs__GFP8ng")
  pvalCol <- paste0("adj.P.Val_", bait, "__vs__GFP8ng")
  
  # filter on PIK3R2 in G12V for AP-MS
  preys <- apms_MPLID25011_DA$Gene.names[apms_MPLID25011_DA[[fcCol]] > 
                                           highThresholds$log2fc
                                         [highThresholds$method == 'AP-MS' & 
                                             highThresholds$node == 'datasetb'] &
                                           apms_MPLID25011_DA[[pvalCol]] < 0.05]
  apms_MPLID25011_df[[bait]] <- preys[!preys %in% bystanders]
  #
  # filter as in other RAS baits in TurboID
  preys <- turboid_MPLID25011_DA$Gene.names[turboid_MPLID25011_DA[[fcCol]] > 
                                  highThresholds$log2fc[
                                    highThresholds$method=='TurboID' & 
                                      highThresholds$node=='datasetb'] &
                                  turboid_MPLID25011_DA[[pvalCol]] < 0.05 ]
  tid_MPLID25011_df[[bait]] <- preys[!preys %in% bystanders]
}


apmsPreyList <- apmsPreyList[grep('Kras', names(apmsPreyList), invert = T )]
tidPreyList <- tidPreyList[grep('Kras', names(tidPreyList), invert = T )]

for (bait in names(apms_MPLID25011_df)) {
  apmsPreyList[[bait]] <- apms_MPLID25011_df[[bait]]
  tidPreyList[[bait]] <- tid_MPLID25011_df[[bait]]
}

saveRDS(apmsPreyList, file.path(interimDataDir(), 'apms_sn19_multiple_ctrl_preys.rds' ))
saveRDS(tidPreyList, file.path(interimDataDir(), 'turboid_sn19_multiple_ctrl_preys.rds' ))

### combine significant preys into 1 files
combineMethods <- function(bait, apmsPreys, turboidPreys, apmsDA, turboidDA ) {
  fcCol <- paste0("logFC_", bait, "__vs__GFP8ng")
  pvalCol <- paste0("adj.P.Val_", bait, "__vs__GFP8ng")
  
  apmsDf <- apmsDA[apmsDA$Gene.names %in% unique(c(apmsPreys, turboidPreys) ), 
                   c('Gene.names', fcCol, pvalCol)]
  turboidDf <- turboidDA[turboidDA$Gene.names %in% unique(c(apmsPreys, turboidPreys) ), 
                         c('Gene.names', fcCol, pvalCol)]
  
  names(apmsDf) <- c('Gene.names', 'logFC_APMS', 'padj_APMS' )
  names(turboidDf) <- c('Gene.names', 'logFC_TbID', 'padj_TbID' )
  
  mergedDf <- merge(apmsDf, turboidDf, by = 'Gene.names', all = T)
  
  mergedDf$significant <- 'core'
  
  mergedDf$significant <- ifelse(
    mergedDf$Gene.names %in% intersect(apmsPreys, turboidPreys),
    'core',
    ifelse(
      mergedDf$padj_APMS < 0.05 & mergedDf$padj_TbID < 0.05 &
        mergedDf$logFC_APMS > 0 & mergedDf$logFC_TbID > 0,
      'shared',
      ifelse(
        mergedDf$Gene.names %in% apmsPreys,
        'AP-MS',
        'TbID'
      )
    )
  )
  mergedDf$significant[is.na(mergedDf$logFC_APMS)] <- 'TbID'
  mergedDf$significant[is.na(mergedDf$logFC_TbID)] <- 'AP-MS'
  mergedDf$bait <- bait
  mergedDf
}

outDf <- NULL
for (bait in names(apmsPreyList[grep('Kras', names(apmsPreyList), invert = T )] )) {
  apmsPreys <- apmsPreyList[[bait]]
  turboidPreys <- tidPreyList[[bait]]
  
  cmbnd <- combineMethods(bait, apmsPreys, turboidPreys, daApms, daTid )
  if (is.null(outDf) ) {
    outDf <- cmbnd
  } else {
    outDf <- rbind(outDf, cmbnd)
  }
}

# outDf <- NULL
for (bait in names(apms_MPLID25011_df) ) {
  apmsPreys <- apmsPreyList[[bait]]
  turboidPreys <- tidPreyList[[bait]]
  
  cmbnd <- combineMethods(bait, apmsPreys, turboidPreys, apms_MPLID25011_DA, turboid_MPLID25011_DA )
  if (is.null(outDf) ) {
    outDf <- cmbnd
  } else {
    outDf <- rbind(outDf, cmbnd)
  }
}

names(outDf)[which(names(outDf) == 'Gene.names' )] <- 'prey'

outDf <- outDf[, c(
  'bait', 'prey', 'significant', 'logFC_APMS', 
  'padj_APMS', 'logFC_TbID', 'padj_TbID')]

write.table(outDf, file.path(interimDataDir(), 'erkpathmap_interactome.tsv'),
            row.names = F, sep = '\t', quote = F
)

# combine complete interactomes (also non sign.) into 1 files
# add quant values of both mwethods into one row

combineAllQuantMethods <- function(bait, apmsPreys, turboidPreys, apmsDA, turboidDA ) {
  fcCol <- paste0("logFC_", bait, "__vs__GFP8ng")
  pvalCol <- paste0("adj.P.Val_", bait, "__vs__GFP8ng")
  
  #apmsDf <- apmsDA
  #turboidDf <- turboidDA
  #
  apmsDf <- apmsDA[, c('Gene.names', fcCol, pvalCol)]
  turboidDf <- turboidDA[, c('Gene.names', fcCol, pvalCol)]
  
  if (length(grep('Kras', bait ) > 0 ) ) {
    apmsDf <- apms_MPLID25011_DA[, c('Gene.names', fcCol, pvalCol)]
    turboidDf <- turboid_MPLID25011_DA[, c('Gene.names', fcCol, pvalCol)]
  }
  
  names(apmsDf) <- c('Gene.names', 'logFC_APMS', 'padj_APMS' )
  names(turboidDf) <- c('Gene.names', 'logFC_TbID', 'padj_TbID' )
  
  mergedDf <- merge(apmsDf, turboidDf, by = 'Gene.names', all = T)
  mergedDf$bait <- bait
  
  mergedDf$prey <- ifelse(
    mergedDf$Gene.names %in% union(apmsPreys, turboidPreys),
    'yes', 'no'
  )
  mergedDf$bait <- bait
  mergedDf
}

outDf <- list()
for (bait in names(apmsPreyList) ) {
  apmsPreys <- apmsPreyList[[bait]]
  turboidPreys <- tidPreyList[[bait]]
  
  cmbnd <- combineAllQuantMethods(bait, apmsPreys, turboidPreys, daApms, daTid )
  outDf[[bait]] <- cmbnd
}
mergedDaLong <- do.call(rbind, outDf)

mergedDaLong$bystander <- ifelse(mergedDaLong$Gene.names %in% bystanders, 'yes', 'no')

write.table(mergedDaLong, file.path(interimDataDir(), 'all_interactions_apms_turboid.tsv'),
            row.names = F, sep = '\t', quote = F
)
