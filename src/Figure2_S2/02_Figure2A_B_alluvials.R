#' @author Sebastian Didusch
#' @date 25.08.05
#' 
#' Visualization script for alluvial plots
#' who is in activ., inactiv., wt network with a focus on paralog specificity

source(here::here("src/utils/utils_path.R"))
library(ggalluvial)

mergedPreys <- list()

for (bait in unique(names(apmsPreys))) {
  mergedPreys[[bait]] <-
    unique(c(apmsPreys[[bait]],
             tidPreys[[bait]]) )
  
}

allPreys <- sort(unique(unlist(mergedPreys)))
preyMat <- t(+sapply(mergedPreys, "%in%", x = allPreys))  ## rasWtMatrix output
colnames(preyMat) <- allPreys
preyMat <- t(preyMat)
preyMat <- as.data.frame(preyMat)
preyMat$Gene.names <- rownames(preyMat)

returnStateBinaryMat <- function(preyMat, apmsDesign, tier) {
  #' function that return binary matrices of wt, inactive, and 
  #' active interactome of a pathway tier, condesing paralogs.
  #' In case of multiple baits for a paralog (e.g NRAS Q61R and   Q61K)
  #' these columns are merged such that preys are enriched in either of the baits
  #' 
  expDes <- apmsDesign[apmsDesign$batch == 1, ]
  expDes <- expDes[expDes$node == tier, ]
  #
  wtBaits <- expDes$groups[expDes$type=='WT']
  inactiveBaits <- expDes$groups[expDes$type=='Inactive']
  #
  activeBaits <- expDes$groups[expDes$type=='Active']
  
  toMerge <- list()
  for (paralog in unique(expDes$paralog)) {
    baits <- expDes$groups[expDes$type=='Active' & expDes$paralog == paralog ]
    if (length(baits ) > 1) toMerge[[paralog]] <- baits
  }
  wtmat <- preyMat[, c(wtBaits)]
  names(wtmat) <- expDes$paralog[match(wtBaits, expDes$groups)]
  #
  inactiveMat <- preyMat[, c(inactiveBaits)]
  names(inactiveMat) <- expDes$paralog[match(inactiveBaits, expDes$groups)]
  #
  activeMat <- preyMat[, activeBaits]
  
  for (paralog in names(toMerge)) {
    baits <- toMerge[[paralog ]]
    activeMat[[paralog]] <- apply(activeMat[, baits], 1, sum )
    activeMat[[paralog]][activeMat[[paralog]] > 1 ] <- 1
    activeMat[, baits] <- NULL
  }
  idx <- which(names(activeMat) %in% expDes$groups)
  names(activeMat)[idx] <- expDes$paralog[match(names(activeMat)[idx], expDes$groups) ]
  
  corder <- c()
  if (tier == 'RAS') {
    corder <- c('KRAS', 'NRAS', 'HRAS')
  } else if (tier == 'RAF') {
    corder <- c('ARAF', 'BRAF', 'RAF1')
  } else if (tier == 'MEK') {
    corder <- c('MEK1', 'MEK2')
  } else {
    corder <- c('ERK1', 'ERK2')
  }
  wtmat <- wtmat[, corder]
  inactiveMat <- inactiveMat[, corder]
  activeMat <- activeMat[, corder]
  
  output <- list(WT=wtmat, Inactive=inactiveMat, Active=activeMat)
  fout <- lapply(output, function(x) x[apply(x,1,sum) > 0,])
  fout$WT$GENE.NAMES <- rownames(fout$WT)
  fout$Inactive$GENE.NAMES <- rownames(fout$Inactive)
  fout$Active$GENE.NAMES <- rownames(fout$Active)
  
  fout
}

determineSet <- function(binaryDf) {
  nc <- ncol(binaryDf) - 1
  
  tmp <- binaryDf[, 1:nc]
  
  singletons <- apply(tmp, 1, sum ) == 1
  sBait <- colnames(tmp)[apply(tmp[singletons, ], 1, which.max) ]
  
  common <- apply(tmp, 1, sum ) == nc
  
  binaryDf$set <- ''
  binaryDf$set[common] <- 'common'
  binaryDf$set[singletons] <- sBait
  
  if (nc == 3) {
    cmbs <- combn(names(tmp), 2)
    
    for (idx in 1:ncol(cmbs)) {
      baits <- cmbs[, idx]
      missing <- setdiff(names(tmp), baits )
      label <- paste(baits, collapse = '')
      
      rnames <- which(tmp[[baits[1] ]] == 1 &
                        tmp[[baits[2] ]] == 1 &
                        tmp[[missing ]] == 0)
      
      if (length(rnames) > 0)
        binaryDf$set[rnames] <- label
    }
  }
  binaryDf
}

prepareAlluvialData <- function(WtDf, InactivDf,  ActivDf) {
  WtDf$type <- 'WT'
  InactivDf$type <- 'Inactive'
  ActivDf$type <- 'Active'
  
  alluvDf <- merge(WtDf[, c("GENE.NAMES", "set", "type")],
                   InactivDf[, c("GENE.NAMES", "set", "type")],
                   by = "GENE.NAMES", all = T, suffixes = c('_WT', '_Inactive')
  )
  
  alluvDf <- merge(alluvDf,
                   ActivDf[, c("GENE.NAMES", "set", "type")],
                   by = "GENE.NAMES", all = T, suffixes = c('', '_Active')
  )
  names(alluvDf)[which(names(alluvDf) == 'type')] <- 'type_Active'
  names(alluvDf)[which(names(alluvDf) == 'set')] <- 'set_Active'
  
  alluvDf$set_WT[is.na(alluvDf$set_WT)] <- 'missing'
  alluvDf$type_WT[is.na(alluvDf$type_WT)] <- 'WT'
  #
  alluvDf$set_Inactive[is.na(alluvDf$set_Inactive)] <- 'missing'
  alluvDf$type_Inactive[is.na(alluvDf$type_Inactive)] <- 'Inactive'
  #
  alluvDf$set_Active[is.na(alluvDf$set_Active)] <- 'missing'
  alluvDf$type_Active[is.na(alluvDf$type_Active)] <- 'WT'
  
  rasCombs <- combn(c('KRAS', 'NRAS', 'HRAS'), 2)
  rasCombs <- c(paste0(rasCombs[, 1], collapse = ''), paste0(rasCombs[, 2], collapse = ''),
    paste0(rasCombs[, 3], collapse = '')
    )
  rafCombs <- combn(c('ARAF', 'BRAF', 'RAF1'), 2)
  rafCombs <- c(paste0(rafCombs[, 1], collapse = ''), paste0(rafCombs[, 2], collapse = ''),
                paste0(rafCombs[, 3], collapse = '')
  )
  
  for (elem in c(rasCombs, rafCombs ) ) {
    alluvDf$set_WT <- gsub(elem, 'shared', alluvDf$set_WT)
    #
    alluvDf$set_Inactive <- gsub(elem, 'shared', alluvDf$set_Inactive)
    #
    alluvDf$set_Active <- gsub(elem, 'shared', alluvDf$set_Active)
  }
  alluvDf
}

prepareAlluvialPlot <- function(allCats, tier, placeholder=T, gapSize=20 ) {
  rasLvls <- c('common',
               'placeholder1',
               'shared',
               'placeholder2',
               'KRAS',
               'placeholder3',
               'NRAS', 
               'placeholder4',
               'HRAS',
               'placeholder5',
               'missing')
  
  # RAF
  rafLvls <- c('common',
               'placeholder1',
               'shared',
               'placeholder2',
               'ARAF',
               'placeholder3',
               'BRAF', 
               'placeholder4',
               'RAF1',
               'placeholder5',
               'missing')
  # MEK
  mekLvls <- c('common',
               'placeholder1',
               'MEK1',
               'placeholder2',
               'MEK2', 
               'placeholder3',
               'missing')
  # ERK
  erkLvls <- c('common',
               'placeholder1',
               'ERK1',
               'placeholder2',
               'ERK2', 
               'placeholder3',
               'missing')
  
  if (!placeholder) {
    rasLvls <- grep('placeholder', rasLvls, invert = T, value = T )
    rafLvls <- grep('placeholder', rafLvls, invert = T, value = T )
    mekLvls <- grep('placeholder', mekLvls, invert = T, value = T )
    erkLvls <- grep('placeholder', erkLvls, invert = T, value = T )
  }
  
  idx <- length(unique(allCats$WT)) - 1
  allCats <- rbind(
    allCats,
    data.frame(WT=c(paste0('placeholder', 1:idx)),
               Inactive=c(paste0('placeholder', 1:idx)),
               Active=c(paste0('placeholder', 1:idx)),
               Freq=rep(gapSize, idx))
  )
  
  lvls <- c()
  if (tier == 'RAS') {
    lvls <- rasLvls
  } else if (tier == 'RAF') {
    lvls <- rafLvls
  } else if (tier == 'MEK') {
    lvls <- mekLvls
  } else {
    lvls <- erkLvls
  }
  
  allCats$Inactive <- factor(
    allCats$Inactive,
    levels = rev(lvls )
  )
  allCats$WT <- factor(
    allCats$WT,
    levels = rev(lvls )
  )
  allCats$Active <- factor(
    allCats$Active,
    levels = rev(lvls )
  )
  allCats
}


rasMats <- returnStateBinaryMat(preyMat, apmsDesign, 'RAS')
rafMats <- returnStateBinaryMat(preyMat, apmsDesign, 'RAF')
mekMats <- returnStateBinaryMat(preyMat, apmsDesign, 'MEK')
erkMats <- returnStateBinaryMat(preyMat, apmsDesign, 'ERK')

rasWtDf <- determineSet(rasMats$WT)
rasInactivDf <- determineSet(rasMats$Inactive)
rasActivDf <- determineSet(rasMats$Active)
rasAlluvial <- prepareAlluvialData(rasWtDf, rasInactivDf, rasActivDf)

# rasAlluvial[rasAlluvial$GENE.NAMES %in% c('ARAF', 'BRAF', 'RAF1'), ]
# # GENE.NAMES set_WT type_WT set_Inactive type_Inactive set_Active type_Active
# # 58        ARAF   NRAS      WT      missing      Inactive     common      Active
# # 106       BRAF common      WT      missing      Inactive     common      Active
# # 723       RAF1 common      WT      missing      Inactive     common      Active

rasPlotData <- as.data.frame(table(rasAlluvial$set_WT,
                                         rasAlluvial$set_Inactive,
                                         rasAlluvial$set_Active
))
names(rasPlotData) <- c('WT', 'Inactive', 'Active', 'Freq')
rasPlotData <- prepareAlluvialPlot(rasPlotData, 'RAS', placeholder = T, gapSize = 50)
rasPlotData <- rasPlotData[rasPlotData$Freq > 0, ]


ggplot(rasPlotData,
       aes(y = Freq,
           axis1 = Inactive, axis2 = WT, axis3 = Active)) +
  geom_alluvium(aes(fill = WT), curve_type='cubic',
                width = 1/32, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = c(common = "#756bb1", shared = "#9e9ac8",
                               KRAS = '#a50f15', NRAS = '#08519c', HRAS = '#31a354',
                               missing='grey80'
  )) +
  guides(fill = "none") +
  geom_stratum(width = 1/32, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) + ylab('# of preys') +
  scale_x_continuous(breaks = 1:3, labels = c("Inactive", "WT", "Active")) +
  theme_cowplot() #+ coord_flip()
ggsave(file.path(figureDir(), 
                 'Figure2/2A_ras_alluvial.pdf'), 
       width = 6, height = 5.5)

rasWtAggCounts <- aggregate(rasPlotData, Freq ~ WT, sum)
rasInactivAggCounts <- aggregate(rasPlotData, Freq ~ Inactive, sum)
rasActivAggCounts <- aggregate(rasPlotData, Freq ~ Active, sum)
names(rasWtAggCounts) <-
  names(rasInactivAggCounts) <-
  names(rasActivAggCounts) <- c('set', 'count')

rasWtAggCounts <- rasWtAggCounts[grep('placeholder', rasWtAggCounts$set, invert = T), ]
rasInactivAggCounts <- rasInactivAggCounts[grep('placeholder', rasInactivAggCounts$set, invert = T), ]
rasActivAggCounts <- rasActivAggCounts[grep('placeholder', rasActivAggCounts$set, invert = T), ]
#
rasWtAggCounts$type <- 'WT'
rasInactivAggCounts$type <- 'Inactive'
rasActivAggCounts$type <- 'Active'

rasAggCounts <- do.call(rbind, 
        list(rasWtAggCounts, rasInactivAggCounts, rasActivAggCounts)
        )

dir.create(file.path(
  figureDir(), 
  'Figure2/data/alluvial/'
), showWarnings = F)

write.table(
  rasAggCounts,
  file.path(figureDir(), 
            'Figure2/data/alluvial/ras_alluvial_set_counts.txt'),
  sep = '\t', row.names = F, quote = F
)
write.table(
  rasAlluvial[, grep('type_', names(rasAlluvial ), invert = T )],
  file.path(figureDir(), 
            'Figure2/data/alluvial/ras_flows_per_prey.txt'),
  sep = '\t', row.names = F, quote = F
)

########## RAF
rafWtDf <- determineSet(rafMats$WT)
rafInactivDf <- determineSet(rafMats$Inactive)
rafActivDf <- determineSet(rafMats$Active)
rafAlluvial <- prepareAlluvialData(rafWtDf, rafInactivDf, rafActivDf)

rafAlluvial[rafAlluvial$GENE.NAMES %in% c('KRAS', 'NRAS', 'HRAS',
                                          'MAP2K1', 'MAP2K2'),
            c('GENE.NAMES', 'set_Inactive', 'set_WT', 'set_Active')]

rafPlotData <- as.data.frame(table(rafAlluvial$set_WT,
                                   rafAlluvial$set_Inactive,
                                   rafAlluvial$set_Active
))
names(rafPlotData) <- c('WT', 'Inactive', 'Active', 'Freq')
rafPlotData <- prepareAlluvialPlot(rafPlotData, 'RAF', placeholder = T, gapSize = 50)
rafPlotData <- rafPlotData[rafPlotData$Freq > 0, ]


ggplot(rafPlotData,
       aes(y = Freq,
           axis1 = Inactive, axis2 = WT, axis3 = Active)) +
  geom_alluvium(aes(fill = WT), curve_type='cubic',
                width = 1/32, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = c(common = "#756bb1", shared = "#9e9ac8",
                               ARAF = '#a50f15', BRAF = '#08519c', RAF1 = '#31a354',
                               missing='grey80'
  )) +
  guides(fill = "none") +
  geom_stratum(width = 1/32, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) + ylab('# of preys') +
  scale_x_continuous(breaks = 1:3, labels = c("Inactive", "WT", "Active")) +
  theme_cowplot() #+ coord_flip()
ggsave(file.path(figureDir(), 
                 'Figure2/2A_raf_alluvial.pdf'), 
       width = 6, height = 5.5)

rafWtAggCounts <- aggregate(rafPlotData, Freq ~ WT, sum)
rafInactivAggCounts <- aggregate(rafPlotData, Freq ~ Inactive, sum)
rafActivAggCounts <- aggregate(rafPlotData, Freq ~ Active, sum)
names(rafWtAggCounts) <-
  names(rafInactivAggCounts) <-
  names(rafActivAggCounts) <- c('set', 'count')

rafWtAggCounts <- rafWtAggCounts[grep('placeholder', rafWtAggCounts$set, invert = T), ]
rafInactivAggCounts <- rafInactivAggCounts[grep('placeholder', rafInactivAggCounts$set, invert = T), ]
rafActivAggCounts <- rafActivAggCounts[grep('placeholder', rafActivAggCounts$set, invert = T), ]
#
rafWtAggCounts$type <- 'WT'
rafInactivAggCounts$type <- 'Inactive'
rafActivAggCounts$type <- 'Active'

rafAggCounts <- do.call(rbind, 
                        list(rafWtAggCounts, rafInactivAggCounts, rafActivAggCounts)
)

write.table(
  rafAggCounts,
  file.path(figureDir(), 
            'Figure2/data/alluvial/raf_alluvial_set_counts.txt'),
  sep = '\t', row.names = F, quote = F
)
write.table(
  rafAlluvial[, grep('type_', names(rafAlluvial ), invert = T )],
  file.path(figureDir(), 
            'Figure2/data/alluvial/raf_flows_per_prey.txt'),
  sep = '\t', row.names = F, quote = F
)
######### MEK
mekWtDf <- determineSet(mekMats$WT)
mekInactivDf <- determineSet(mekMats$Inactive)
mekActivDf <- determineSet(mekMats$Active)
mekAlluvial <- prepareAlluvialData(mekWtDf, mekInactivDf, mekActivDf)

mekAlluvial[mekAlluvial$GENE.NAMES %in% c('KRAS', 'NRAS', 'HRAS',
                                          'ARAF', 'BRAF', 'RAF1',
                                          'MAP2K1', 'MAP2K2', 'MAPK3', 'MAPK1'),
            c('GENE.NAMES', 'set_Inactive', 'set_WT', 'set_Active')]

mekPlotData <- as.data.frame(table(mekAlluvial$set_WT,
                                   mekAlluvial$set_Inactive,
                                   mekAlluvial$set_Active
))
names(mekPlotData) <- c('WT', 'Inactive', 'Active', 'Freq')
mekPlotData <- prepareAlluvialPlot(mekPlotData, 'MEK', placeholder = T, gapSize = 50)
mekPlotData <- mekPlotData[mekPlotData$Freq > 0, ]


ggplot(mekPlotData,
       aes(y = Freq,
           axis1 = Inactive, axis2 = WT, axis3 = Active)) +
  geom_alluvium(aes(fill = WT), curve_type='cubic',
                width = 1/32, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = c(common = "#756bb1", shared = "#9e9ac8",
                               MEK1 = '#a50f15', MEK2 = '#31a354',
                               missing='grey80'
  )) +
  guides(fill = "none") +
  geom_stratum(width = 1/32, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) + ylab('# of preys') +
  scale_x_continuous(breaks = 1:3, labels = c("Inactive", "WT", "Active")) +
  theme_cowplot() #+ coord_flip()
ggsave(file.path(figureDir(), 
                 'Figure2/2A_mek_alluvial.pdf'), 
       width = 6, height = 5.5)

mekWtAggCounts <- aggregate(mekPlotData, Freq ~ WT, sum)
mekInactivAggCounts <- aggregate(mekPlotData, Freq ~ Inactive, sum)
mekActivAggCounts <- aggregate(mekPlotData, Freq ~ Active, sum)
names(mekWtAggCounts) <-
  names(mekInactivAggCounts) <-
  names(mekActivAggCounts) <- c('set', 'count')

mekWtAggCounts <- mekWtAggCounts[grep('placeholder', mekWtAggCounts$set, invert = T), ]
mekInactivAggCounts <- mekInactivAggCounts[grep('placeholder', mekInactivAggCounts$set, invert = T), ]
mekActivAggCounts <- mekActivAggCounts[grep('placeholder', mekActivAggCounts$set, invert = T), ]
#
mekWtAggCounts$type <- 'WT'
mekInactivAggCounts$type <- 'Inactive'
mekActivAggCounts$type <- 'Active'

mekAggCounts <- do.call(rbind, 
                        list(mekWtAggCounts, mekInactivAggCounts, mekActivAggCounts)
)

write.table(
  mekAggCounts,
  file.path(figureDir(), 
            'Figure2/data/alluvial/mek_set_counts.txt'),
  sep = '\t', row.names = F, quote = F
)
write.table(
  mekAlluvial[, grep('type_', names(mekAlluvial ), invert = T )],
  file.path(figureDir(), 
            'Figure2/data/alluvial/mek_flows_per_prey.txt'),
  sep = '\t', row.names = F, quote = F
)
######## ERK
erkWtDf <- determineSet(erkMats$WT)
erkInactivDf <- determineSet(erkMats$Inactive)
erkActivDf <- determineSet(erkMats$Active)
erkAlluvial <- prepareAlluvialData(erkWtDf, erkInactivDf, erkActivDf)

erkAlluvial[erkAlluvial$GENE.NAMES %in% c('KRAS', 'NRAS', 'HRAS',
                                          'ARAF', 'BRAF', 'RAF1',
                                          'MAP2K1', 'MAP2K2', 'MAPK3', 'MAPK1'),
            c('GENE.NAMES', 'set_Inactive', 'set_WT', 'set_Active')]

erkPlotData <- as.data.frame(table(erkAlluvial$set_WT,
                                   erkAlluvial$set_Inactive,
                                   erkAlluvial$set_Active
))
names(erkPlotData) <- c('WT', 'Inactive', 'Active', 'Freq')
erkPlotData <- prepareAlluvialPlot(erkPlotData, 'ERK', placeholder = T, gapSize = 50)
erkPlotData <- erkPlotData[erkPlotData$Freq > 0, ]


ggplot(erkPlotData,
       aes(y = Freq,
           axis1 = Inactive, axis2 = WT, axis3 = Active)) +
  geom_alluvium(aes(fill = WT), curve_type='cubic',
                width = 1/32, knot.pos = 0, reverse = FALSE) +
  scale_fill_manual(values = c(common = "#756bb1", shared = "#9e9ac8",
                               ERK1 = '#a50f15', ERK2 = '#31a354',
                               missing='grey80'
  )) +
  guides(fill = "none") +
  geom_stratum(width = 1/32, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) + ylab('# of preys') +
  scale_x_continuous(breaks = 1:3, labels = c("Inactive", "WT", "Active")) +
  theme_cowplot() #+ coord_flip()
ggsave(file.path(figureDir(), 
                 'Figure2/2A_erk_alluvial.pdf'), 
       width = 6, height = 5.5)

erkWtAggCounts <- aggregate(erkPlotData, Freq ~ WT, sum)
erkInactivAggCounts <- aggregate(erkPlotData, Freq ~ Inactive, sum)
erkActivAggCounts <- aggregate(erkPlotData, Freq ~ Active, sum)
names(erkWtAggCounts) <-
  names(erkInactivAggCounts) <-
  names(erkActivAggCounts) <- c('set', 'count')

erkWtAggCounts <- erkWtAggCounts[grep('placeholder', erkWtAggCounts$set, invert = T), ]
erkInactivAggCounts <- erkInactivAggCounts[grep('placeholder', erkInactivAggCounts$set, invert = T), ]
erkActivAggCounts <- erkActivAggCounts[grep('placeholder', erkActivAggCounts$set, invert = T), ]
#
erkWtAggCounts$type <- 'WT'
erkInactivAggCounts$type <- 'Inactive'
erkActivAggCounts$type <- 'Active'

erkAggCounts <- do.call(rbind, 
                        list(erkWtAggCounts, erkInactivAggCounts, erkActivAggCounts)
)

write.table(
  erkAggCounts,
  file.path(figureDir(), 
            'Figure2/data/alluvial/erk_set_counts.txt'),
  sep = '\t', row.names = F, quote = F
)
write.table(
  erkAlluvial[, grep('type_', names(erkAlluvial ), invert = T )],
  file.path(figureDir(), 
            'Figure2/data/alluvial/erk_flows_per_prey.txt'),
  sep = '\t', row.names = F, quote = F
)
#######
######
#####
##########
#########
##########
##########
##########
#########
##########
activeUpset <- read.table(file.path(figureDir(), 
                                    'S2/data/combined_active_interactome_binary_matrix.tsv'),
                          header = T, sep = '\t')

rasrafmekActivPrey <- activeUpset$Gene.names[
  (activeUpset$ERK==0 &  activeUpset$MEK==1 & activeUpset$RAF==1 & activeUpset$RAS==1) |
    (activeUpset$ERK==0 &  activeUpset$MEK==1 & activeUpset$RAF==1 & activeUpset$RAS==0) |
    (activeUpset$ERK==0 &  activeUpset$MEK==0 & activeUpset$RAF==1 & activeUpset$RAS==1)]

##
rasrafmekOnlyActivPrey <- activeUpset$Gene.names[
  (activeUpset$ERK==0 &  activeUpset$MEK==1 & activeUpset$RAF==1 & activeUpset$RAS==1)]

rasrafActivPreys <-  activeUpset$Gene.names[(activeUpset$RAF==1 & activeUpset$RAS==1)]
rafmekActivPreys <-  activeUpset$Gene.names[(activeUpset$MEK==1 & activeUpset$RAF==1)]

rasActivVert <- rasAlluvial[rasAlluvial$GENE.NAMES %in% rasrafmekActivPrey, ]
rafActivVert <- rafAlluvial[rafAlluvial$GENE.NAMES %in% rasrafmekActivPrey, ]
mekActivVert <- mekAlluvial[mekAlluvial$GENE.NAMES %in% rasrafmekActivPrey, ]

allDfVert <- merge(rasActivVert[, c('GENE.NAMES', 'set_Active' )],
      rafActivVert[, c('GENE.NAMES', 'set_Active' )],
      by = 'GENE.NAMES', suffixes = c('_ras', '_raf'), all = T)
allDfVert <- merge(allDfVert,
                   mekActivVert[, c('GENE.NAMES', 'set_Active' )],
                   by = 'GENE.NAMES', suffixes = c('_ras', '_raf'),
                   all = T)

names(allDfVert)[which(names(allDfVert)=='set_Active')] <- 'set_Active_mek'
allDfVert$set_Active_ras[is.na(allDfVert$set_Active_ras)] <- 'missing'
allDfVert$set_Active_raf[is.na(allDfVert$set_Active_raf)] <- 'missing'
allDfVert$set_Active_mek[is.na(allDfVert$set_Active_mek)] <- 'missing'
allDfVert$set_Active_raf[allDfVert$set_Active_raf == 'common'] <- 'common2'
allDfVert$set_Active_raf[allDfVert$set_Active_raf == 'shared'] <- 'shared2'
allDfVert$set_Active_mek[allDfVert$set_Active_mek == 'common'] <- 'common3'
allDfVert$set_Active_mek[allDfVert$set_Active_mek == 'missing'] <- 'missing3'

allCats <- as.data.frame(table(allDfVert$set_Active_ras,
                               allDfVert$set_Active_raf,
                               allDfVert$set_Active_mek
))
names(allCats) <- c('RAS', 'RAF', 'MEK', 'Freq')
allCats <- allCats[allCats$Freq > 0, ]

tmp <- data.frame(RAS=c('placeholder1', 'placeholder2','placeholder3', 'placeholder4', 'placeholder5'),
           RAF=c('placeholder1a', 'placeholder2a', 'placeholder3a', 'placeholder4a', 'placeholder5a'),
           MEK=c('placeholder1b','placeholder2b', 'placeholder3b', 'trash1', 'trash2'),
           Freq=c(50, 50, 50, 50, 50)
           )
allCats <- rbind(allCats, tmp)


allCats$RAS <- factor(allCats$RAS, levels = rev(
  c(
    'common',
    'placeholder1',
    #'placeholder2',
    'shared',
    'placeholder2',
    #'placeholder3',
    'KRAS',
    'placeholder3',
    #'placeholder4',
    'NRAS',
    'placeholder4',
    'HRAS',
    'placeholder5',
    'missing'
  )
))
#
allCats$RAF <- factor(
  allCats$RAF,
  levels = rev(c('common2',
                 'placeholder1a',
                 'shared2',
                 'placeholder2a',
                 'ARAF',
                 'placeholder3a',
                 'BRAF',
                 'placeholder4a',
                 'RAF1',
                 'placeholder5a'
  ) )
)

allCats$MEK <- factor(allCats$MEK, levels = rev(
  c(
    'common3',
    'placeholder1b',
    'MEK1',
    'placeholder2b',
    #'placeholder2',
    'MEK2',
    'placeholder3b',
    'trash2',
    'trash1',
    'missing3'
  )
))
##
ggplot(allCats,
       aes(y = Freq,
           axis1 = RAS, axis2 = RAF, axis3 = MEK)) +
  geom_alluvium(aes(fill = RAS), curve_type='cubic',
                width = 1/32, knot.pos = 0, reverse = F) +
  # scale_fill_manual(values = c(common = "#756bb1", shared = "#9e9ac8",
  #                              MEK1 = '#a50f15', MEK2 = '#31a354',
  #                              missing='grey80'
  # )) +
  guides(fill = "none") +
  geom_stratum(width = 1/32, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) + ylab('# of preys') +
  scale_x_continuous(breaks = 1:3, labels = c("RAS", "RAF", "MEK")) +
  theme_cowplot() #+ coord_flip()
ggsave(file.path(figureDir(), 
                 'Figure2/2B_rasrafmek_active_alluvial.pdf'), 
       width = 6, height = 5.5)

rasAggCounts <- aggregate(allCats, Freq ~ RAS, sum)
rafAggCounts <- aggregate(allCats, Freq ~ RAF, sum)
mekAggCounts <- aggregate(allCats, Freq ~ MEK, sum)
names(rasAggCounts) <-
  names(rafAggCounts) <-
  names(mekAggCounts) <- c('set', 'count')

rasAggCounts <- rasAggCounts[grep('placeholder', rasAggCounts$set, invert = T), ]
rafAggCounts <- rafAggCounts[grep('placeholder', rafAggCounts$set, invert = T), ]
mekAggCounts <- mekAggCounts[grep('placeholder', mekAggCounts$set, invert = T), ]
#
rasAggCounts$type <- 'RAS'
rafAggCounts$type <- 'RAF'
mekAggCounts$type <- 'MEK'

rasrafmekAggCounts <- do.call(rbind, 
                        list(rasAggCounts, rafAggCounts, mekAggCounts)
)
write.table(
  rasrafmekAggCounts[grep('missing|trash', rasrafmekAggCounts$set, invert = T), ],
  file.path(figureDir(), 
            'Figure2/data/alluvial/rasrafmek_active_set_counts.txt'),
  sep = '\t', row.names = F, quote = F
)
write.table(
  allDfVert[, grep('type_', names(allDfVert ), invert = T )],
  file.path(figureDir(), 
            'Figure2/data/alluvial/rasrafmek_active_flows_per_prey.txt'),
  sep = '\t', row.names = F, quote = F
)



### 



