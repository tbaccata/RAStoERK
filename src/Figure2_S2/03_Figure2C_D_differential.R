#' @author Sebastian Didusch
#' @date 25.07.02
#' 
#' define groups to highlight for perturbation network
#' 

library(igraph)
library(ggpubr)

source(here::here("src/utils/utils_path.R"))

mergedPerturbedGres <- read.table(
  file.path(figureDir(), 'Figure2/data/gost_merged_perturbations.tsv' ),
  header = T, sep = '\t'
)
mergedPerturbations <- read.table(
  file.path(figureDir(), 'Figure2/data/merged_perturbations.tsv'  ),
  header = T, sep = '\t'
)
tidPerturbations <- read.table(
  file.path(figureDir(), 'Figure2/data/turboid_perturbations.tsv' ),
  header = T, sep = '\t')
  
apmsPerturbations <- read.table(
  file.path(figureDir(), 'Figure2/data/apms_perturbations.tsv' ),
  header = T, sep = '\t')


rasRegulated <- read.table(file.path(dataDir(), 'integration/RASassocProteins.tsv'),
                           header = T, sep = '\t')


rasSignaling <- mergedPerturbedGres$intersection[mergedPerturbedGres$term_name == 
                                                   'Ras signaling pathway'] %>% 
  strsplit(",") %>% unlist() %>% unique()
rasSignaling <- c(rasRegulated$Protein, rasSignaling)

# only interested in direct ras signaling effectors and close-by regulators
toExlude <- c("KRAS", "NRAS","HRAS", "ARAF", 'BRAF', 
              'RAF1', "MAP2K1","MAP2K2", 'MAPK3', 'MAPK1', "KSR1","KSR2")

chaps1433s <- c('YWHAB','YWHAE','YWHAG','YWHAH','YWHAQ','YWHAZ','SFN')

rtkSignaling <- mergedPerturbedGres$intersection[mergedPerturbedGres$term_name == 
                                                   'Signaling by Receptor Tyrosine Kinases'] %>% 
  strsplit(",") %>% unlist() %>% unique()

pi3kSignaling <- mergedPerturbedGres$intersection[mergedPerturbedGres$term_name == 
                                                    'PI3K-Akt signaling pathway'] %>% 
  strsplit(",") %>% unlist() %>% unique()

rasSignaling <- rasSignaling[!rasSignaling %in% c(chaps1433s, toExlude)]
rtkSignaling <- rtkSignaling[!rtkSignaling %in% c(chaps1433s, toExlude) ]
pi3kSignaling <- pi3kSignaling[!pi3kSignaling %in% c(chaps1433s, toExlude) ]

toExclude <- intersect(pi3kSignaling, rtkSignaling)
rtkSignaling <- rtkSignaling[!rtkSignaling %in% toExclude]

pi3kSignaling <- pi3kSignaling[!pi3kSignaling %in% c("CDC37", "MAP2K1",
                                                     "MAP2K2","HSP90AA1",
                                                     "HSP90B1", "HSP90AB1",
                                                     "HRAS", "KRAS", 'RAF1'
) ]

annotProc <- list('signaling by RTKs' = rtkSignaling,
                  'PI3K-Akt signaling pathway' = pi3kSignaling,
                  '14-3-3 proteins' = chaps1433s
)

alreadyAnnot <- annotProc %>% unlist() %>% unique()

annotProc[['RAS signaling']] <- setdiff(rasSignaling, alreadyAnnot)

#######################
#### Perturbation Plots

lvls <- c(baitOrder[baitOrder %in% unique(apmsDesign$groups[apmsDesign$type=='Active' ])],
          baitOrder[baitOrder %in% unique(apmsDesign$groups[apmsDesign$type=='Inactive' ]) ]
)

mergedPerturbations$Mutant <- factor(mergedPerturbations$Mutant,
                                     levels =  (lvls) )
apmsPerturbations$Mutant <- factor(apmsPerturbations$Mutant,
                                   levels =  (lvls) )
tidPerturbations$Mutant <- factor(tidPerturbations$Mutant,
                                  levels =  (lvls) )

#plotPert <- mergedPerturbations[mergedPerturbations$significant!='incongruent',]
#plotPert$Process <- ''
#
apmsPlotPert <- apmsPerturbations
apmsPlotPert$Process <- ''
#
tidPlotPert <- tidPerturbations
tidPlotPert$Process <- ''
annotList <- list()
for (proc in names(annotProc)) {
  print(proc)
  # plotPert$Process[which(plotPert$Gene.names %in% annotProc[[proc]]) ] <- proc
  #
  apmsPlotPert$Process[which(apmsPlotPert$Gene.names %in% annotProc[[proc]]) ] <- proc
  #
  tidPlotPert$Process[which(tidPlotPert$Gene.names %in% annotProc[[proc]]) ] <- proc
  
  annotList[[proc]] <- data.frame(Protein=annotProc[[proc]], Process=proc)
}
annotDf <- do.call(rbind, annotList)

###
### ap-ms perturbation plot
ggplot(apmsPlotPert, aes(x=Mutant, y=logFC, color=Process, label=Gene.names )) + 
  geom_jitter(data = apmsPlotPert[apmsPlotPert$perturbed=='unchanged',],
              width = 0.2, alpha=0.25, size=0.5, color='grey60') +
  geom_jitter(data = apmsPlotPert[apmsPlotPert$perturbed!='unchanged' &
                                    apmsPlotPert$Process == "",],
              width = 0.2, color='grey60', size=0.5) +
  geom_jitter(data = apmsPlotPert[apmsPlotPert$Process != "",],
              width = 0.2) +
  geom_hline(yintercept = c(-1, 1), lty=2) +
  # geom_vline(xintercept = c(-1, 1), lty=2) +
  # geom_text_repel(data = apmsPlotPert[apmsPlotPert$Gene.names %in% toLabel,],
  #                    color='black', show.legend = F) +
  ylim(min(min(tidPlotPert$logFC), min(apmsPlotPert$logFC) ),
        max(max(tidPlotPert$logFC), max(apmsPlotPert$logFC) ) ) +
  scale_color_manual(values = c('14-3-3 proteins'='#8da0cb',
                                'PI3K-Akt signaling pathway'='#a6d854',
                                'RAS signaling'='#fc8d62',
                                'signaling by RTKs'='#e78ac3')) +
  theme_cowplot() + ylab('logFC(Mut/WT)') + xlab('') +
  ggtitle('AP-MS') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ), legend.position = 'none' )
ggsave(file.path(figureDir(), 'Figure2/2D_apms_differential_annotated_swapped.pdf' ), 
       width = 5.5, height = 4)


### tbid perturbation plot
ggplot(tidPlotPert, aes(x=Mutant, y=logFC, color=Process, label=Gene.names )) + 
  geom_jitter(data = tidPlotPert[tidPlotPert$perturbed=='unchanged',],
              width = 0.2, alpha=0.25, size=0.5, color='grey60') +
  geom_jitter(data = tidPlotPert[tidPlotPert$perturbed!='unchanged' &
                                   tidPlotPert$Process == "",],
              width = 0.2, color='grey60', size=0.5) +
  geom_jitter(data = tidPlotPert[tidPlotPert$Process != "",],
              width = 0.2) +
  geom_hline(yintercept = c(-1, 1), lty=2) + 
  # geom_text_repel(data = tidPlotPert[tidPlotPert$Gene.names %in% toLabel,],
  #                 color='black', show.legend = F) +
  ylim(min(min(tidPlotPert$logFC), min(apmsPlotPert$logFC) ),
        max(max(tidPlotPert$logFC), max(apmsPlotPert$logFC) ) ) +
  scale_color_manual(values = c('14-3-3 proteins'='#8da0cb',
                                'PI3K-Akt signaling pathway'='#a6d854',
                                'RAS signaling'='#fc8d62',
                                'signaling by RTKs'='#e78ac3')) +
  theme_cowplot() + ylab('logFC(Mut/WT)') + xlab("") +
  ggtitle('TbID') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ), legend.position = 'none' )
ggsave(file.path(figureDir(), 'Figure2/2D_turboid_differential_annotated_swapped.pdf' ), 
       width = 5.5, height = 4)

write.table(
  apmsPlotPert,
  file.path(figureDir(), 'Figure2/data/apms_differential_annotated.tsv'),
  row.names = F, sep = '\t', quote = F
)
#
write.table(
  tidPlotPert,
  file.path(figureDir(), 'Figure2/data/turboid_differential_annotated.tsv'),
  row.names = F, sep = '\t', quote = F
)

#### Congruency between AP-MS and TbID
# ### congruence
# aInc <- mergedPerturbations[mergedPerturbations$significant=='opposite' & abs(mergedPerturbations$logFC_APMS) > 1, ]
# tInc <- mergedPerturbations[mergedPerturbations$significant=='opposite' & abs(mergedPerturbations$logFC_TbID) > 1, ]
# perturbedIncongruent <- rbind(tInc, aInc)
# #perturbedIncongruent$perturbed <- 'incongruent'
# 
# apmsInconLoss <- unique(perturbedIncongruent$Gene.names[perturbedIncongruent$logFC_APMS < -1] )
# apmsInconGain <- unique(perturbedIncongruent$Gene.names[perturbedIncongruent$logFC_APMS > 1] )
# tidInconLoss <- unique(perturbedIncongruent$Gene.names[perturbedIncongruent$logFC_TbID < -1] )
# tidInconGain <- unique(perturbedIncongruent$Gene.names[perturbedIncongruent$logFC_TbID > 1] )

nInconruentTid <- sum(abs(mergedPerturbations$logFC_TbID[mergedPerturbations$significant=='opposite']) > 1 )
nIncongruentApms <- sum(abs(mergedPerturbations$logFC_APMS[mergedPerturbations$significant=='opposite']) > 1 )
nIncongruent <- nIncongruentApms + nInconruentTid

nCongruent <- sum(mergedPerturbations$perturbed!='no' & mergedPerturbations$significant=='congruent')
nApms <- sum(mergedPerturbations$perturbed!='no' & mergedPerturbations$significant=='AP-MS')
# nApmsOnly <- sum(mergedPerturbations$perturbed!='unchanged' & 
#                    mergedPerturbations$significant=='AP-MS' & 
#                    is.na(mergedPerturbations$logFC_TbID))
# nApms <- nApms - nApmsOnly
nTid <- sum(mergedPerturbations$perturbed!='unchanged' & mergedPerturbations$significant=='TbID')
# nTidOnly <- sum(mergedPerturbations$perturbed!='unchanged' & 
#                   mergedPerturbations$significant=='TbID' &
#                   is.na(mergedPerturbations$logFC_APMS))
# nTid <- nTid - nTidOnly

dfCongruency <- data.frame(
  differential=c('congruent', 'opposite', 'AP-MS',  'TbID'),
  number=c(nCongruent, nIncongruent, nApms, nTid)
)
# dfCongruency <- data.frame(
#   perturbed=c('congruent', 'opposite', 'AP-MS', 'only id. in AP-MS', 'TbID', 'only id. in TbID'),
#   number=c(nCongruent, nIncongruent, nApms, nApmsOnly, nTid, nTidOnly)
# )
dfCongruency$pct <- 100 * dfCongruency$number/sum(dfCongruency$number)

dfCongruency$Interactome <- 'Differential interactors'

dfCongruency$differential <- factor(dfCongruency$differential,
                                 levels = c('congruent', 'opposite', 
                                            'AP-MS', 
                                            #'only id. in AP-MS',
                                            'TbID'
                                            #,'only id. in TbID'
                                            ) )

ggplot(dfCongruency, aes(y=Interactome, x=pct, fill = differential ) ) +
  geom_bar(stat = 'identity', position = 'stack') + ylab('') + xlab('% of diff. interactions (Mut/WT)') +
  scale_fill_manual(values = c(congruent='#7fc97f', 
                               opposite='#f0027f', 
                               'AP-MS'='#beaed4', 
                               #'only id. in AP-MS'='#decbe4',
                               TbID='#386cb0'
                               #,  'only id. in TbID'='#b3cde3'
                               ) 
                    ) +
  theme_cowplot() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )
ggsave(
  file.path(figureDir(), 'Figure2/2C_differential_congruence_by_method.pdf'),
  width = 6, height = 2.5
)

write.table(
  dfCongruency,
  file.path(figureDir(), 'Figure2/data/differential_congruence_by_method.tsv'),
  row.names = F, sep = '\t', quote = F
)
