#' @author Sebastian Didusch
#' @date 25.06.16
#' 
#' global analysis of AP-MS and TurboID data by PCA

source(here::here("src/utils/utils_path.R"))


### merge intensties into one file
apmsSe <- readRDS(file.path(processedDataDir(), 'apms/se_apms.rds'))
apms_MPLID25011_se <- readRDS(file.path(processedDataDir(), 'apms/se_apms_MPLID25011.rds'))
apmsPreys <- readRDS(file.path(processedDataDir(), 'apms/sn19_multiple_ctrl_preys.rds'))

tidSe <- readRDS(file.path(processedDataDir(), 'turboid/se_turboid.rds'))
tid_MPLID25011_se <- readRDS(file.path(processedDataDir(), 'turboid/se_turboid_MPLID25011.rds'))
tidPreys <- readRDS(file.path(processedDataDir(), 'turboid/sn19_multiple_ctrl_preys.rds'))

mergeMeasurements <- function(se1, se2, preyList1, preyList2, bc1=F, bc2=F, MPLID25011_2=F ) {
  uniquePreys1 <- unique(unlist(preyList1))
  uniquePreys2 <- unique(unlist(preyList2))
  
  fdata1 <- rowData(se1)
  fdata2 <- rowData(se2)
  
  # fdata1 <- fdata1[fdata1$Gene.names %in% union(uniquePreys1, uniquePreys2), ]
  # fdata2 <- fdata2[fdata2$Gene.names %in% union(uniquePreys1, uniquePreys2), ]
  
  impName1 <- grep('Imp', assayNames( se1 ), value = T)
  impName2 <- grep('Imp', assayNames( se2 ), value = T)
  
  impLfqs1 <- assay(se1, impName1)[rownames(fdata1),]
  impLfqs2 <- assay(se2, impName1)[rownames(fdata2),]
  
  names(impLfqs1) <- paste0(names(impLfqs1), '_data1')
  names(impLfqs2) <- paste0(names(impLfqs2), '_data2')
  
  impLfqs1$Gene.names <- fdata1$Gene.names
  impLfqs2$Gene.names <- fdata2$Gene.names
  
  impLfqs1$ProteinId_data1 <- rownames(fdata1)
  impLfqs2$ProteinId_data2 <- rownames(fdata2)
  
  mergedLfqs <- merge(impLfqs1, impLfqs2, by='Gene.names', all=T )
  
  rnames <- ifelse(is.na(mergedLfqs$ProteinId_data1), 
                   mergedLfqs$ProteinId_data2, 
                   mergedLfqs$ProteinId_data1 )
  mergedLfqs$ProteinId_data1 <- NULL
  mergedLfqs$ProteinId_data2 <- NULL
  
  imputedMergedLfqsBc <- mergedLfqs
  
  mLfqs1 <- mergedLfqs[, grep('data1', names(mergedLfqs) )]
  mLfqs2 <- mergedLfqs[, grep('data2', names(mergedLfqs) )]
  
  impMLE1 <- impute.MinProb(as.matrix(mLfqs1), q = 0.01)
  impMLE1 <- as.data.frame(impMLE1)
  
  impMLE2 <- impute.MinProb(as.matrix(mLfqs2), q = 0.01)
  impMLE2 <- as.data.frame(impMLE2)
  
  expDesign <- colData(se1)
  expDesign2 <- colData(se2)
  
  if (bc1) {
    impMLE1 <- removeBatchEffect(impMLE1, batch = expDesign$batch, group = expDesign$groups)
    impMLE1 <- as.data.frame(impMLE1)
  }
  if (bc2) {
    impMLE2 <- removeBatchEffect(impMLE2, batch = expDesign2$batch, group = expDesign2$groups)
    impMLE2 <- as.data.frame(impMLE2)
  }
  
  if (MPLID25011_2) {
    groupsOfInterest <- paste0(expDesign$samples[expDesign$paralog!='KRAS'], '_data1')
    edesignOfInterest <- expDesign[expDesign$paralog!='KRAS', ]
    impMLE1 <- impMLE1[, groupsOfInterest]
    
    imputedMergedLfqs <- cbind(impMLE1, impMLE2)
    tmp1 <- edesignOfInterest[, c("samples", "groups")]
    tmp1$batch <- 1
    tmp2 <- expDesign2[, c("samples", "groups")]
    tmp2$batch <- 2
    mergedDesign <- rbind(
      tmp1,
      tmp2
    )
    imputedMergedLfqsBc <- removeBatchEffect(imputedMergedLfqs,
                                             batch = mergedDesign$batch, group = mergedDesign$groups)
    imputedMergedLfqsBc <- as.data.frame(imputedMergedLfqsBc)
    
  } else {
    imputedMergedLfqsBc <- cbind(impMLE1, impMLE2)
  }
  imputedMergedLfqsBc$Gene.names <- mergedLfqs$Gene.names
  imputedMergedLfqsBc
}


apmsMergedLfqs <- mergeMeasurements(apmsSe, apms_MPLID25011_se, apmsPreys, apmsPreys, bc1=T, bc2=F, MPLID25011_2=T )
tidMergedLfqs <- mergeMeasurements(tidSe, tid_MPLID25011_se, tidPreys, tidPreys, bc1=T, bc2=F, MPLID25011_2=T )

mergeDesigns <- function(apmsSe, apms_MPLID25011_se, apmsMerged ) {
  tmp1 <- colData(apmsSe)
  tmp1 <- tmp1[paste0(tmp1$samples, '_data1') %in% names(apmsMerged), ]
  tmp1$samples <- paste0(tmp1$samples, '_data1')
  tmp2 <- colData(apms_MPLID25011_se)
  tmp2$samples <- paste0(tmp2$samples, '_data2')
  tmp2$Replicate <- NULL
  tmp2$batch <- 100
  tmp2$type <- 'Active'
  tmp2$type[tmp2$groups=='KrasWT'] <- 'WT'
  tmp2$type[tmp2$groups=='KrasS17N'] <- 'Inactive'
  tmp2$type[tmp2$groups=='GFP8ng'] <- 'Control'
  tmp2$node <- 'RAS'
  tmp2$paralog <- 'KRAS'
  
  rbind(tmp1, tmp2)
}


apmsMergedDesign <- mergeDesigns(apmsSe, apms_MPLID25011_se, apmsMerged)
tidMergedDesign <- mergeDesigns(tidSe, tid_MPLID25011_se, tidMerged)

apmsMergedLfqs <- apmsMergedLfqs[apmsMergedLfqs$Gene.names %in% unique(unlist(apmsPreys ) ),
                                 c(apmsMergedDesign$samples[apmsMergedDesign$type != 'Control'],
                                   'Gene.names')]
#
tidMergedLfqs <- tidMergedLfqs[tidMergedLfqs$Gene.names %in% unique(unlist(tidPreys ) ),
                               c(tidMergedDesign$samples[tidMergedDesign$type != 'Control'],
                                 'Gene.names')]

annotation_colors$shape <- c(WT=15, CI=6, CA=17, CA2=14)

tmpDat <- computePCA(apmsMergedLfqs, apmsMergedDesign[apmsMergedDesign$type != 'Control', ] )

apmsPcaDf <- tmpDat$df_out
apmsPcaDf$bait <- apmsMergedDesign[apmsMergedDesign$type!='Control',]$groups
apmsPcaDf$paralog <- apmsMergedDesign[apmsMergedDesign$type!='Control',]$paralog
apmsPcaDf$batch <- apmsMergedDesign[apmsMergedDesign$type!='Control',]$batch
apmsPcaDf$node <- apmsMergedDesign[apmsMergedDesign$type!='Control',]$node
apmsPcaDf$state <- apmsMergedDesign[apmsMergedDesign$type!='Control',]$type
percentage <- tmpDat$percentage
#
x <- gsub('PC', '', percentage)
x <- gsub('.*\\( ', '', x)
x <- gsub('%.*', '', x)

apmsX <- as.numeric(x)

apmsPcaDf$batch <- factor(apmsPcaDf$batch)
trueNodes <- apmsPcaDf$node

### highlight all

to_high <- apmsPcaDf #apmsPcaDf[apmsPcaDf$batch==1,]
to_high <- to_high[to_high$bait %in% c(
  'KrasG12V', 'NrasQ61R', 'HrasQ61R', 'BrafV600E', 'Raf1S259F'
),]

apmsPcaDf$tier <- factor(apmsPcaDf$node, levels = c(
  'RAS',
  'RAF',
  'MEK',
  'ERK'
))

apmsPcaDf$paralog <- factor(apmsPcaDf$paralog, levels = c(
  'KRAS',
  'NRAS',
  'HRAS',
  'ARAF', 'BRAF', 'RAF1',
  'MEK1', 'MEK2', 'ERK1', 'ERK2'
))


secondBaits <- c('KrasG12V', 'NrasQ61R', 'HrasQ61R', 'BrafV600E', 'Raf1S259F')
apmsPcaDf$state[apmsPcaDf$group %in% secondBaits] <- 'CA2'
apmsPcaDf$state[apmsPcaDf$state=='Active'] <- 'CA'
apmsPcaDf$state[apmsPcaDf$state=='Inactive'] <- 'CI'

apmsPcaDf$state <- factor(apmsPcaDf$state,
                          levels = c('CI', 'WT', 'CA', 'CA2'))

p <-
  ggplot(apmsPcaDf, aes(
    x = PC1,
    y = PC2,
    color=paralog,
    shape = state,
    label = bait
  ))
p <-
  p + geom_point(size=3) + 
  xlab(percentage[1]) + 
  ylab(percentage[2]) +
  scale_color_manual(values = annotation_colors$paralog ) +
  scale_shape_manual(values = annotation_colors$shape ) +
  #geom_text_repel(data = to_high, show.legend = F ) +
  ggtitle('AP-MS') +
  #ggtitle(paste0('PCA (', nrow(dat), ' prey proteins)' )) +
  theme_cowplot()
p


ggsave(file.path(
  figureDir(),
  'S1/apms_interactome_pca.pdf'
), width = 6, height = 5 )

write.table(
  apmsPcaDf,
  file.path(
    figureDir(),
    'S1/data/apms_interactome_pca.tsv'
  ), row.names = F, sep = '\t'
)


### TurboID
tmpDat <- computePCA(tidMergedLfqs, tidMergedDesign[tidMergedDesign$type != 'Control', ] )

tidPcaDf <- tmpDat$df_out
tidPcaDf$bait <- tidMergedDesign[tidMergedDesign$type!='Control',]$groups
tidPcaDf$paralog <- tidMergedDesign[tidMergedDesign$type!='Control',]$paralog
tidPcaDf$batch <- tidMergedDesign[tidMergedDesign$type!='Control',]$batch
tidPcaDf$node <- tidMergedDesign[tidMergedDesign$type!='Control',]$node
tidPcaDf$state <- tidMergedDesign[tidMergedDesign$type!='Control',]$type
percentage <- tmpDat$percentage
#
x <- gsub('PC', '', percentage)
x <- gsub('.*\\( ', '', x)
x <- gsub('%.*', '', x)
tidX <- as.numeric(x)

tidPcaDf$batch <- factor(tidPcaDf$batch)
trueNodes <- tidPcaDf$node

### highlight all

to_high <- tidPcaDf[tidPcaDf$batch==1,]
to_high <- to_high[grep('WT|NLS|LYN', to_high$bait),]

tidPcaDf$node <- factor(tidPcaDf$node, levels = c(
  'RAS',
  'RAF',
  'MEK',
  'ERK'
  ##'PM',
  #$'NUC'
))

tidPcaDf$paralog <- factor(tidPcaDf$paralog, levels = c(
  'KRAS',
  'NRAS',
  'HRAS',
  'ARAF', 'BRAF', 'RAF1',
  'MEK1', 'MEK2', 'ERK1', 'ERK2'
  ##'PM',
  #$'NUC'
))

tidPcaDf$state[tidPcaDf$group %in% secondBaits] <- 'CA2'
tidPcaDf$state[tidPcaDf$state=='Active'] <- 'CA'
tidPcaDf$state[tidPcaDf$state=='Inactive'] <- 'CI'

p <-
  ggplot(tidPcaDf, aes(
    x = PC1,
    y = PC2,
    color=paralog,
    shape = state
    #label = bait
  ))
p <-
  p + geom_point(size=3) + 
  xlab(percentage[1]) + 
  ylab(percentage[2]) +
  scale_color_manual(values = annotation_colors$paralog ) +
  scale_shape_manual(values = annotation_colors$shape ) +
  #ggrepel::geom_text_repel(data = to_high, show.legend = F ) +
  ggtitle('TurboID') +
  #ggtitle(paste0('PCA (', nrow(dat), ' prey proteins)' )) +
  theme_cowplot()
p

ggsave(file.path(
  figureDir(),
  'S1/tbid_interactome_pca.pdf'
), width = 6, height = 5 )


write.table(
  tidPcaDf,
  file.path(
    figureDir(),
    'S1/data/tbid_interactome_pca.tsv'
  ), row.names = F, sep = '\t'
)
