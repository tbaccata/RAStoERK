#' @author Sebastian Didusch
#' date 25.03.11
#' 
#' AP-MS filtering, quantification and differential abundance analysis pipeline.
#' 
#' There are two datasets: 
#' First processed, MPLID23168 which contains the majority (31 out of 35)
#' of bait samples. 
#' 
#' The second dataset is MPLID25011, which contains 4 baits that had to be
#' remeasured because of technical reasons. However, they were 1) measured with
#' the same facility under almost the same experimental conditions and setup, 
#' and 2) processed exactly the same in Spectronaut and in this script.


## define paths
source(here::here("src/utils/utils_path.R"))

# Spectronaut Protein LFQ report
fnameMPLID23168 <- file.path(dataDir(), 'raw/apms_SN19.8/sn_protein_lfq/MPLID23168_APMS_SN198_2025fasta_search_Report_MPL_Protein_LFQ_Report (Pivot).tsv')

fDesign <- file.path(dataDir(), 'raw/apms_SN19.8/sn_protein_lfq/MPLID23168_APMS_SN198_2025fasta_search_ConditionSetup.tsv')
protData <- read.delim(fnameMPLID23168, sep = "\t", header = T, stringsAsFactors = F)


#### 1 define exp. expDesign
conditions <- read.table(fDesign, sep = '\t', header = T, comment.char = '')
expDesign <- conditions[, c("Run.Label", "Condition", "Replicate")]

names(expDesign) <- c('samples', 'groups', 'batch')
expDesign$type <- ''
expDesign$type[grep('WT', expDesign$groups, ignore.case = T )] <- 'WT'
expDesign$type[grep('DD', expDesign$groups, ignore.case = T )] <- 'Active'
expDesign$type[grep('AA|S17N|D594N|D486A|D447A|D447N', expDesign$groups, ignore.case = T )] <- 'Inactive'
expDesign$type[grep('G12|V600E|Q61|G13|S214P|S257L|S259F|G469A', expDesign$groups, ignore.case = T )] <- 'Active'
expDesign$type[grep('LYN|NLS|GFP', expDesign$groups, ignore.case = T )] <- 'Control'
expDesign$groups <- gsub('GFP2ng.*', 'GFP2ng', expDesign$groups)
expDesign$groups <- gsub('GFP8ng.*', 'GFP8ng', expDesign$groups)

expDesign$node <- ''
expDesign$node[grep('RAS', expDesign$groups, ignore.case = T)] <- 'RAS'
expDesign$node[grep('RAF', expDesign$groups, ignore.case = T)] <- 'RAF'
expDesign$node[grep('MEK', expDesign$groups, ignore.case = T)] <- 'MEK'
expDesign$node[grep('ERK', expDesign$groups, ignore.case = T)] <- 'ERK'

expDesign$node[grep('LYN', expDesign$groups, ignore.case = T)] <- 'PM'
expDesign$node[grep('GFP', expDesign$groups, ignore.case = T)] <- 'GFP'
expDesign$node[grep('NLS', expDesign$groups, ignore.case = T)] <- 'NUC'

expDesign$paralog <- 'Control'
expDesign$paralog[grep('kras', expDesign$groups, ignore.case = T)] <- 'KRAS'
expDesign$paralog[grep('nras', expDesign$groups, ignore.case = T)] <- 'NRAS'
expDesign$paralog[grep('hras', expDesign$groups, ignore.case = T)] <- 'HRAS'
expDesign$paralog[grep('araf', expDesign$groups, ignore.case = T)] <- 'ARAF'
expDesign$paralog[grep('braf', expDesign$groups, ignore.case = T)] <- 'BRAF'
expDesign$paralog[grep('raf1', expDesign$groups, ignore.case = T)] <- 'RAF1'

expDesign$paralog[grep('mek1', expDesign$groups, ignore.case = T)] <- 'MEK1'
expDesign$paralog[grep('mek2', expDesign$groups, ignore.case = T)] <- 'MEK2'
expDesign$paralog[grep('erk1', expDesign$groups, ignore.case = T)] <- 'ERK1'
expDesign$paralog[grep('erk2', expDesign$groups, ignore.case = T)] <- 'ERK2'
####


#### 2 DIFF. ABUNDANCE ANALYSIS

method <- 'VSN'

protData <- read.delim(fnameMPLID23168, sep = "\t", header = T, stringsAsFactors = F)

raws <- protData[, grep('Quantity', names(protData) )]

# exclude controls and baits not used in this dataset
baitGroups <- unique(grep('GFP|NLS|LYN11|Kras', expDesign$groups, value = T, invert = T))
se <- preProcess(protData, expDesign, 32, baitGroups, minMSMS=2, minRazor=2, minValueGroup=3)

normDf <- assay(se, 'LFQIntensity')[isQuantRnames(se),]
normDf[is.na(normDf)] <- NA # NaN != NA
#lfqs <- normDf[apply(normDf, 1, function(x) !any(is.na(x)) ),]
lfqs <- renormalizeIntensities(normDf, method)

imputationDf <- lfqs

se <- setAssay(x = se,
               assay = lfqs,
               assayName = paste0('Norm_', method) )


impMLE <- impute.MinProb(as.matrix(lfqs), q = 0.01)
impMLE <- as.data.frame(impMLE)


imputationDf[row.names(lfqs), ] <- impMLE

se <-
  setAssay(x = se,
           assay = imputationDf,
           assayName = paste0('Imp_Norm_', method ) )


# filtered row data as df
filtData <-
  rowData(se)[isQuantRnames(se),]

countIdx <-
  grep("razorUniqueCount.", colnames(filtData))
pep.count.table = data.frame(
  count = apply(filtData[, countIdx], 1, FUN = min),
  row.names = rownames(filtData)
)
pep.count.table$count = pep.count.table$count + 1


contrastMatrix <- data.frame(group1=baitGroups, group2='GFP8ng')

contrastMatrix <- rbind(
  contrastMatrix,
  data.frame(group1=c('LYN11', 'NLS','LYN11', 'NLS', 'GFP8ng'), group2=c('GFP8ng', 'GFP8ng', 'GFP2ng', 'GFP2ng', 'GFP2ng'))
)


groupMut <- c()
groupWt <- c()
for (bait in unique(expDesign$groups)) {
  # ignore WTs, and baits not used in this dataset, and controls
  if (length(grep('WT|Kras|LYN|NLS|GFP', bait )) > 0 ) next
  
  para <- unique(expDesign$paralog[expDesign$groups==bait])
  
  wt <- unique( grep('WT', expDesign$groups[expDesign$paralog==para], value = T )  )
  inactive <- unique( expDesign$groups[expDesign$paralog==para & expDesign$type=='Inactive']  )
  
  # if (inactive == bait) next
  
  # print(paste(bait, inactive))
  
  groupMut <- c(groupMut, bait)
  groupWt <- c(groupWt, wt)
  
}
wtContrasts <- data.frame(group1=groupMut, group2=groupWt)

groupActive <- c()
groupInactive <- c()
for (elem in unique(expDesign$paralog)) {
  if (elem == 'Control') next
  activeMut <- unique(expDesign$groups[expDesign$paralog==elem & expDesign$type == 'Active'])
  inactiveMut <- unique(expDesign$groups[expDesign$paralog==elem & expDesign$type == 'Inactive'])
  
  for (a in activeMut) {
    print(paste(a, inactiveMut))
    groupActive <- c(groupActive, a)
    groupInactive <- c(groupInactive, inactiveMut)
  }
  if (length(activeMut) == 2) {
    #print(paste(activeMut[1], activeMut[2] ))
    groupActive <- c(groupActive, activeMut[1] )
    groupInactive <- c(groupInactive, activeMut[2] )
  }
}

wtContrasts <- rbind(
  wtContrasts,
  data.frame(group1=groupActive, group2=groupInactive)
)


contrastMatrix <- rbind(
  contrastMatrix,
  data.frame(group1=baitGroups, group2='GFP2ng')
)



dataMatrix <- as.matrix(assay(se, paste0('Imp_Norm_', method ))[isQuantRnames(se),])

tmpOut <- groupComparisonsBatchCorrected(dataMatrix,
                                         contrastMatrix,
                                         expDesign,
                                         pep.count.table)
tmpOut$Gene.names <- filtData[, geneName, drop = T]

wtOut <- groupComparisonsBatchCorrected(dataMatrix,
                                        wtContrasts,
                                        expDesign,
                                        pep.count.table)
wtOut$Gene.names <- filtData[, geneName, drop = T]


dir.create(file.path(dataDir(), 'processed/apms'), showWarnings = FALSE)
write.table(tmpOut,
            file.path(dataDir(), 'processed/apms/deqms_apms.tsv'),
            row.names = F, sep = '\t', quote = F
)
write.table(wtOut,
            file.path(dataDir(), 'processed/apms/deqms_apms_mut_vs_wt.tsv'),
            row.names = F, sep = '\t', quote = F
)
write.table(expDesign,
            file.path(dataDir(), 'processed/apms/design_apms.txt'),
            row.names = F, sep = '\t', quote = F
)
saveRDS(se, file.path(dataDir(), 'processed/apms/se_apms.rds'))

#### process data set B the same way (vsn, pre-processing, imputation, etc.)
MPLID25011Design <- read.table(file.path(dataDir(),
                                   'raw/apms_SN19.8/sn_protein_lfq/MPLID25011_SN198_2025fasta_APMS_search_ConditionSetup.tsv'),
                         header = T, sep = '\t', comment.char = "")

MPLID25011Design <- MPLID25011Design[, c("Run.Label", "Condition", "Replicate")]
names(MPLID25011Design)[1] <- 'samples'
names(MPLID25011Design)[2] <- 'groups'
MPLID25011Design$groups[MPLID25011Design$groups=='GFP'] <- 'GFP8ng'
MPLID25011Design$samples <- make.names(MPLID25011Design$groups, unique = T)
MPLID25011Design$groups <- gsub('KRAS_', 'Kras', MPLID25011Design$groups)

fnameDatasetB <- file.path(dataDir(), 'raw/apms_SN19.8/sn_protein_lfq/MPLID25011_SN198_2025fasta_APMS_search_Report_MPL_Protein_LFQ_Report (Pivot).tsv')

protDataKras <- read.delim(fnameDatasetB, sep = "\t", header = T, stringsAsFactors = F)
MPLID25011SummExp <- preProcess(protDataKras, MPLID25011Design, 32, 
                     unique(MPLID25011Design$groups[MPLID25011Design$groups != 'GFP8ng']),
                     minMSMS=2, minRazor=2, minValueGroup=3)

normDf <- assay(MPLID25011SummExp, 'LFQIntensity')[isQuantRnames(MPLID25011SummExp),]
normDf[is.na(normDf)] <- NA # NaN != NA
#lfqs <- normDf[apply(normDf, 1, function(x) !any(is.na(x)) ),]
lfqs <- renormalizeIntensities(normDf, "VSN")

imputationDf <- lfqs

mse <- setAssay(x = MPLID25011SummExp,
                assay = lfqs,
                assayName = paste0('Norm_VSN') )

impMLE <- impute.MinProb(as.matrix(lfqs), q = 0.01)
impMLE <- as.data.frame(impMLE)


imputationDf[row.names(lfqs), ] <- impMLE

MPLID25011SummExp <-
  setAssay(x = MPLID25011SummExp,
           assay = imputationDf,
           assayName = paste0('Imp_Norm_VSN' ) )

# filtered row data as df
filtData <-
  rowData(MPLID25011SummExp)[isQuantRnames(MPLID25011SummExp),]

countIdx <-
  grep("razorUniqueCount.", colnames(filtData))
pep.count.table = data.frame(
  count = apply(filtData[, countIdx], 1, FUN = min),
  row.names = rownames(filtData)
)
pep.count.table$count = pep.count.table$count + 1

baitGroups <- unique(grep('GFP|NLS|LYN11', colData(MPLID25011SummExp)$groups, value = T, invert = T))

contrastMatrix <- data.frame(group1=baitGroups, group2='GFP8ng')

dataMatrix <- as.matrix(assay(MPLID25011SummExp, 'Imp_Norm_VSN' )[isQuantRnames(MPLID25011SummExp),])

tmpOut <- groupComparisons(dataMatrix,
                                         contrastMatrix,
                                         colData(MPLID25011SummExp),
                                         pep.count.table)
tmpOut$Gene.names <- filtData[, geneName, drop = T]

MPLID25011_differential_contrasts <- data.frame(
  group1=c(
    'KrasG12V', 'KrasG12D', 'KrasS17N'
  ), group2=c(
    'KrasWT', 'KrasWT', 'KrasWT'
  )
)
MPLID25011_differential <- groupComparisons(dataMatrix,
                              MPLID25011_differential_contrasts,
                           colData(MPLID25011SummExp),
                           pep.count.table)
MPLID25011_differential$Gene.names <- filtData[, geneName, drop = T]

write.table(tmpOut,
            file.path(dataDir(), 'processed/apms/deqms_apms_MPLID25011.tsv'),
            row.names = F, sep = '\t', quote = F
)
write.table(MPLID25011_differential,
            file.path(dataDir(), 'processed/apms/deqms_apms_MPLID25011_mut_vs_wt.tsv'),
            row.names = F, sep = '\t', quote = F
)
write.table(MPLID25011Design,
            file.path(dataDir(), 'processed/apms/kras/design_apms_MPLID25011.txt'),
            row.names = F, sep = '\t', quote = F
)
# saveRDS(MPLID25011SummExp, file.path(dataDir(), 'processed/apms/se_apms_MPLID25011.rds'))
