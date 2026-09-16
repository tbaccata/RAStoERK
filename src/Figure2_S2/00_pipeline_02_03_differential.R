#' @author Sebastian Didusch
#' @date 25.04.17
#' 
#' define perturbations by comparing mutant vs. wt

source(here::here("src/utils/utils_path.R"))
library(ggpubr)

wtDesign <- apmsDesign[apmsDesign$batch==1 & apmsDesign$type=='WT',]
activDesign <- apmsDesign[apmsDesign$batch==1 & apmsDesign$type=='Active',]
inactivDesign <- apmsDesign[apmsDesign$batch==1 & apmsDesign$type=='Inactive',]

apmsMutVsWt <- read.table(file.path(dataDir(), 'processed/apms/deqms_apms_mut_vs_wt.tsv'),
                          header = T, sep = '\t')
apmsKrasMutWt <- read.table(file.path(dataDir(), 'processed/apms/kras/deqms_apms_kras_mut_vs_wt.tsv'),
                            header = T, sep = '\t')
##
tidMutVsWt <- read.table(file.path(dataDir(), 'processed/turboid/deqms_turboid_mut_vs_wt.tsv'),
                         header = T, sep = '\t')
tidKrasMutWt <- read.table(file.path(dataDir(), 'processed/turboid/kras/deqms_turboid_kras_mut_vs_wt.tsv'),
                           header = T, sep = '\t')

wtComps <- grep('logFC.*WT', names(apmsMutVsWt), value = T)

activVsWtComps <- grep(paste0(activDesign$groups, collapse = '|'), wtComps, value = T )

createPerturbationData <- function(activVsWtComps, onlyWt=F) {
  apmsPertDf <- tidPertDf <- NULL
  
  pctWtSpec <- pctMutSpec <- mutation <- c()
  
  for (elem in activVsWtComps) {
    pvalCol <- gsub('logFC','adj.P.Val', elem)
    
    groups <- unlist(strsplit(elem, '__vs__'))
    mut <- gsub('logFC_', '', groups[1] )
    wt <- groups[2]
    
    apmsSelectedPreys <- unique(unlist(apmsPreys[c(wt, mut)]))
    tidSelectedPreys <- unique(unlist(tidPreys[c(wt, mut)]))
    
    wtPreys <- interactome$prey[interactome$bait==wt]
    mutPreys <- interactome$prey[interactome$bait==mut]
    
    wtSpec <- setdiff(wtPreys, mutPreys )
    mutSpec <- setdiff(mutPreys, wtPreys )
    
    n_wtSpec <- length(wtSpec)/length(union(wtPreys, mutPreys))
    n_mutSpec <- length(mutSpec)/length(union(wtPreys, mutPreys))
    
    pctWtSpec <- c(pctWtSpec, n_wtSpec)
    pctMutSpec <- c(pctMutSpec, n_mutSpec)
    mutation <- c(mutation, mut)
    
    if (onlyWt) {
      apmsSelectedPreys <- unique(unlist(apmsPreys[c(wt)]))
      tidSelectedPreys <- unique(unlist(tidPreys[c(wt)]))
    }
    
    apmsPertData <- apmsMutVsWt[apmsMutVsWt$Gene.names %in% apmsSelectedPreys,
                                c('Gene.names', elem, pvalCol) ] 
    tidPertData <- tidMutVsWt[tidMutVsWt$Gene.names %in% tidSelectedPreys,
                              c('Gene.names', elem, pvalCol) ]
    
    if (length(grep('Kras', wt ) ) > 0 ) {
      apmsPertData <- apmsKrasMutWt[apmsKrasMutWt$Gene.names %in% apmsSelectedPreys,
                                    c('Gene.names', elem, pvalCol) ] 
      tidPertData <- tidKrasMutWt[tidKrasMutWt$Gene.names %in% tidSelectedPreys,
                                  c('Gene.names', elem, pvalCol) ]
    }
    
    names(apmsPertData) <- c('Gene.names', 'logFC', 'padj' )
    names(tidPertData) <- c('Gene.names', 'logFC', 'padj' )
    
    apmsPertData$Comp <- elem
    apmsPertData$Mutant <- mut
    
    apmsPertData$perturbed <-  ifelse(apmsPertData$logFC > 1 & 
                                        apmsPertData$padj < 0.05, 
                                      'gain', ifelse(
                                        apmsPertData$logFC < -1 & 
                                          apmsPertData$padj < 0.05,
                                        'loss',
                                        'unchanged'
                                      ))
    
    tidPertData$Comp <- elem
    tidPertData$Mutant <- mut
    
    tidPertData$perturbed <-  ifelse(tidPertData$logFC > 1 & 
                                       tidPertData$padj < 0.05, 
                                     'gain', ifelse(
                                       tidPertData$logFC < -1 & 
                                         tidPertData$padj < 0.05,
                                       'loss',
                                       'unchanged'
                                     ))
    
    if (is.null(apmsPertDf)) {
      apmsPertDf <- apmsPertData
      tidPertDf <- tidPertData
    } else {
      apmsPertDf <- rbind(apmsPertDf, apmsPertData) 
      tidPertDf <- rbind(tidPertDf, tidPertData)
    }
  }
  gofLofDf <- data.frame(Mutation=mutation,
                         WT_specific=pctWtSpec,
                         MUT_specific=pctMutSpec )
  gofLofDf$shared <- gofLofDf$WT_specific + gofLofDf$MUT_specific
  gofLofDf$shared <- 1 - gofLofDf$shared
  
  return(list(apms=apmsPertDf, tid=tidPertDf, quant=gofLofDf))
}

tmp <- createPerturbationData(wtComps, onlyWt = F)
apmsPertDf <- tmp$apms
tidPertDf <- tmp$tid

write.table(apmsPertDf, file.path(figureDir(), 'Figure2/data/apms_perturbations.tsv' ), row.names = F, sep = '\t', quote = F )
write.table(tidPertDf, file.path(figureDir(), 'Figure2/data/turboid_perturbations.tsv' ), row.names = F, sep = '\t', quote = F )

########### Combine Methods in data.frame
### merge AP-MS and TurboID perturbaton by max log2FC method (minP)
combineMethods <- function(baits, apsmSelectedPreys, turboidPreys, apmsDA, turboidDA ) {
  fcCol <- paste0("logFC_", baits[1], "__vs__", baits[2])
  pvalCol <- paste0("adj.P.Val_", baits[1], "__vs__", baits[2])
  
  apmsDf <- apmsDA[apmsDA$Gene.names %in% unique(c(apsmSelectedPreys, turboidPreys) ), 
                   c('Gene.names', fcCol, pvalCol)]
  turboidDf <- turboidDA[turboidDA$Gene.names %in% unique(c(apsmSelectedPreys, turboidPreys) ), 
                         c('Gene.names', fcCol, pvalCol)]
  
  if (length(grep('Kras', fcCol ) ) > 0 ) {
    apmsDf <- apmsKrasMutWt[apmsKrasMutWt$Gene.names %in% unique(c(apsmSelectedPreys, turboidPreys) ),
                            c('Gene.names', fcCol, pvalCol)]
    turboidDf <- tidKrasMutWt[tidKrasMutWt$Gene.names %in% unique(c(apsmSelectedPreys, turboidPreys) ),
                              c('Gene.names', fcCol, pvalCol)]
  }
  
  names(apmsDf) <- c('Gene.names', 'logFC_APMS', 'padj_APMS' )
  names(turboidDf) <- c('Gene.names', 'logFC_TbID', 'padj_TbID' )
  
  mergedDf <- merge(apmsDf, turboidDf, by = 'Gene.names', all = T)
  
  mergedDf$significant <- 'core'
  
  mergedDf$significant <- ifelse(
    (mergedDf$padj_APMS < 0.05 & mergedDf$padj_TbID < 0.05 &
       mergedDf$logFC_APMS > 0 & mergedDf$logFC_TbID > 0) |
      mergedDf$padj_APMS < 0.05 & mergedDf$padj_TbID < 0.05 &
      mergedDf$logFC_APMS < 0 & mergedDf$logFC_TbID < 0,
    'congruent',
    ifelse(
      (mergedDf$padj_APMS < 0.05 & mergedDf$padj_TbID < 0.05 &
         mergedDf$logFC_APMS > 0 & mergedDf$logFC_TbID < 0) |
        mergedDf$padj_APMS < 0.05 & mergedDf$padj_TbID < 0.05 &
        mergedDf$logFC_APMS < 0 & mergedDf$logFC_TbID > 0,
      'opposite',
      ifelse(
        (is.na(mergedDf$logFC_APMS) & mergedDf$padj_TbID < 0.05) |
          (mergedDf$padj_APMS > 0.05 & mergedDf$padj_TbID < 0.05),
        'TbID',
        ifelse(
          (is.na(mergedDf$logFC_TbID) & mergedDf$padj_APMS < 0.05) |
            (mergedDf$padj_APMS < 0.05 & mergedDf$padj_TbID > 0.05),
          'AP-MS',
          'nsig'
        )
      )
    )
  )
  mergedDf$significant[is.na(mergedDf$logFC_APMS)] <- 'TbID'
  mergedDf$significant[is.na(mergedDf$logFC_TbID)] <- 'AP-MS'
  mergedDf$Mutant <- baits[1]

  mergedDf
}

out <- list()
for(idx in seq_along(wtComps) ) {
  baits <- comps[[idx]]
  apmsSelectedPreys <- unique(unlist(apmsPreys[baits ]))
  tidSelectedPreys <- unique(unlist(tidPreys[baits ]))
  
  cmbnd <- combineMethods(baits, apmsSelectedPreys, tidSelectedPreys, apmsMutVsWt, tidMutVsWt )
  out[[baits[1]]] <- cmbnd
}

mergedPerturbations <- do.call(rbind, out)

mergedPerturbations$perturbed <- ifelse(
  (!is.na(mergedPerturbations$logFC_APMS) & 
            abs(mergedPerturbations$logFC_APMS) > 1 &
                    mergedPerturbations$padj_APMS < 0.05
            ) |
    (!is.na(mergedPerturbations$logFC_TbID) & 
       abs(mergedPerturbations$logFC_TbID) > 1 &
       mergedPerturbations$padj_TbID < 0.05
    ),
  'yes',
  'no'
   )

mergedPerturbations$direction <- ifelse(
  mergedPerturbations$perturbed=='no',
  '-',
  ifelse(
    (!is.na(mergedPerturbations$logFC_APMS) &
       mergedPerturbations$logFC_APMS > 1 &
       mergedPerturbations$padj_APMS < 0.05
       ) |
    (!is.na(mergedPerturbations$logFC_TbID) &
       mergedPerturbations$logFC_TbID > 1 &
       mergedPerturbations$padj_TbID < 0.05
    ),
   'gain',
   'loss')
)
mergedPerturbations$direction[mergedPerturbations$perturbed=='yes' &
                                mergedPerturbations$significant=='opposite'] <- 'opposite'

mergedPertList <- list()
for (mutant in unique(mergedPerturbations$Mutant )) {
  pertGained <- mergedPerturbations$Gene.names[mergedPerturbations$Mutant == mutant &
                                                 mergedPerturbations$direction=='gain']
  pertGained <- gsub(';.*', '', pertGained)
  #
  pertLoss <- mergedPerturbations$Gene.names[mergedPerturbations$Mutant == mutant &
                                               mergedPerturbations$direction=='loss']
  pertLoss <- gsub(';.*', '', pertLoss)
  
  mergedPertList[[paste0(mutant, '_gain' )]] <- pertGained
  mergedPertList[[paste0(mutant, '_loss' )]] <- pertLoss
}

background <- unique(c(apmsKrasMutWt$Gene.names, apmsMutVsWt$Gene.names, 
tidKrasMutWt$Gene.names, tidMutVsWt) )
background <- gsub(';.*', '', background)

library(gprofiler2)
mergedPerturbedGostRes <- gost(mergedPertList, sources = c('REAC', 'KEGG', 'GO:BP'),
                               evcodes = T,
                               custom_bg = background)

write.table(mergedPerturbations, file.path(figureDir(), 'Figure2/data/merged_perturbations.tsv' ),
            row.names = F, sep = '\t', quote = F )
toFile <- mergedPerturbedGostRes$result
toFile$parents <- NULL

write.table(mergedPerturbations, file.path(figureDir(), 'Figure2/data/merged_perturbations.tsv' ),
            row.names = F, sep = '\t', quote = F )
write.table(toFile, file.path(figureDir(), 'Figure2/data/gost_merged_perturbations.tsv' ),
            row.names = F, sep = '\t', quote = F )
