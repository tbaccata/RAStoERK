#' @author Sebastian Didusch
#' @date 25.06.07
#' 
#' extract core pathway enrichment
#' where are core (bait - bait) PPIs enriched?
#' does paralog have main downstream targets?

source(here::here("src/utils/utils_path.R"))

pcts <- c()
ns <- c()
methods <- c()
nodes <- c()
baits <- c()
paras <- c()
knownInteractors <- novelPreys <- list()
for (bait in unique(apmsDesign$groups[apmsDesign$type != 'Control']) ) {

  para <- unique(apmsDesign$paralog[apmsDesign$groups==bait])
  node <- unique(apmsDesign$node[apmsDesign$groups==bait])
  
  geneName <- ''
  
  if ( para %in% c('KRAS', 'NRAS', 'HRAS', 'ARAF', 'BRAF', 'RAF1' ) ) {
    geneName <- para
  } else if (para == 'MEK1') {
    geneName <- 'MAP2K1'
  } else if (para == 'MEK2') {
    geneName <- 'MAP2K2'
  } else if (para == 'ERK1') {
    geneName <- 'MAPK3'
  } else if (para == 'ERK2') {
    geneName <- 'MAPK1'
  }
  
  aPreys <- gsub(';.*', '', apmsPreys[[bait]] )
  tPreys <- gsub(';.*', '', tidPreys[[bait]] )
  cmbPreys <- unique(c(aPreys, tPreys))
  #### summary
  
  paras <- c(paras, c(para, para, para ))
  nodes <- c(nodes, c(node, node, node) )
  baits <- c(baits, c(bait, bait, bait) )
  
  interactors <- names(neighbors(ppi, geneName))
  apmsKnown <- length(intersect(interactors, aPreys))
  tidKnown <- length(intersect(interactors, tPreys))
  cmbKnown <- length(intersect(interactors, cmbPreys))
  
  knownInteractors[[paste0('AP-MS_',bait)]] <- c(geneName, intersect(interactors, aPreys))
  knownInteractors[[paste0('TbID_',bait)]] <- c(geneName, intersect(interactors, tPreys))
  novelPreys[[bait]] <- c(geneName, setdiff(c(unique(aPreys, tPreys)), interactors  ) )
  
  nApms <- length(aPreys)
  nTid <- length(tPreys)
  nCmb <- length(cmbPreys)
  ns <-c(ns, c(nApms, nTid, nCmb))
  methods <- c(methods, c('AP-MS', 'TurboID', 'Combined' ))
  
  apmsPctKnown <- apmsKnown/length(aPreys)
  tidPctKnown <- tidKnown/length(tPreys)
  tidPctKnown <- tidKnown/length(tPreys)
  cmbPctKnown <- cmbKnown/length(cmbPreys)
  pcts <- c(pcts, c(apmsPctKnown, tidPctKnown, cmbPctKnown) )
}

df <- data.frame(bait=baits, paralog=paras, node=nodes,method=methods, preys=ns, known=pcts)

tmp <- df

tmp$nknown <- tmp$preys * tmp$known
tmp$new <- tmp$preys - tmp$nknown

x <- tmp[, -which( names(tmp ) == 'nknown' )]
names(x)[which( names(x) == 'new' )] <- 'number'
x$interactor <- 'novel'

y <- tmp[, -which( names(tmp ) == 'new' )]
names(y)[which( names(y) == 'nknown' )] <- 'number'
y$interactor <- 'known'

plotDf <- rbind(x, y)
plotDf <- plotDf[order(plotDf$preys, decreasing = T), ]
ordering <- plotDf[plotDf$method=='AP-MS', 'bait']
ordering <- ordering[!duplicated(ordering)]
plotDf$bait <- factor(plotDf$bait, levels = ordering)
##
# ggplot(plotDf[plotDf$method != 'Combined',],
#        aes(x=bait, y=number, fill=interactor)) +
#   geom_bar(position = 'stack', stat = 'identity') +
#   scale_fill_manual(values=c(known='#A6CEE3', novel='#1F78B4') ) +
#   labs(y = "% of preys") + xlab('') +
#   theme_cowplot() +
# theme(
#   axis.text.x = element_text(
#     angle = 90,
#     hjust = 1,
#     vjust = 0.5
#   )) + facet_wrap(~method, scales = 'free_x', nrow = 2)
# ggsave(file.path(figureDir(), 'S1/stacked_all_baits_preys_known.pdf'),
#        width = 7,
#        height = 6)

plotDf$known <- plotDf$known * 100
#write.table(plotDf[plotDf$method != 'Combined',],
#            file.path(figureDir(), 'S1/data/stacked_all_baits_preys_known.tsv'),
#       row.names = F, sep = '\t', quote = F )

sumPctKnownApms <- sum(plotDf$preys[plotDf$method=='AP-MS'] * 
                         plotDf$known[plotDf$method=='AP-MS']) / sum(plotDf$preys[plotDf$method=='AP-MS'])

sumPctNovelApms <- 100 - sumPctKnownApms

sumPctKnownTid <- sum(plotDf$preys[plotDf$method=='TurboID'] * 
                        plotDf$known[plotDf$method=='TurboID']) / sum(plotDf$preys[plotDf$method=='TurboID'])
sumPctNovelTid <- 100 - sumPctKnownTid

sumPctKnownCmb <- sum(plotDf$preys[plotDf$method=='Combined'] * 
                        plotDf$known[plotDf$method=='Combined']) / sum(plotDf$preys[plotDf$method=='Combined'])
sumPctNovelCmb <- 100 - sumPctKnownCmb

summaryPlotData <- data.frame(
  method = c('AP-MS', 'AP-MS', 'TurboID', 'TurboID', 'Combined', 'Combined'),
  pct = c(
    sumPctKnownApms,
    sumPctNovelApms,
    sumPctKnownTid,
    sumPctNovelTid,
    sumPctKnownCmb,
    sumPctNovelCmb
  ),
  prey = c('known', 'novel', 'known', 'novel',  'known', 'novel')
)
###
ggplot(summaryPlotData, aes(x=method, y=pct, fill=prey) ) +
  geom_bar(position = 'stack', stat = 'identity') + theme_cowplot() + 
  scale_fill_manual(values=c(known='#A6CEE3', novel='#1F78B4') ) +
  ylab('% preys') + xlab('')
ggsave(
  file.path(figureDir(), 'Figure1/method_pct_known_novel.pdf' ),
  width = 4, height = 4.5
)
write.table(
  summaryPlotData,
  file.path(figureDir(), 'Figure1/data/method_pct_known_novel.tsv' ),
  row.names = F, sep = '\t', quote = F
)
######

### collapse mutants into 1 entity per protein
collapsedDf <- list()
for (paralog in unique(apmsDesign$paralog[apmsDesign$type != 'Control'] ) ) {
    paraBaits <- unique(apmsDesign$groups[apmsDesign$paralog == paralog])
    
    apmsParaPreys <- unique(unlist(apmsPreys[paraBaits ]))
    tidParaPreys <- unique(unlist(tidPreys[paraBaits ]))
    
     
    apmsKnowns <- tidKnowns <- c()
    for (bait in paraBaits ) {
      apmsKnowns <- c(apmsKnowns, knownInteractors[[paste0('AP-MS_', bait ) ]] )
      tidKnowns <- c(tidKnowns, knownInteractors[[paste0('TbID_', bait ) ]] )
    }
    apmsKnowns <- unique(apmsKnowns)
    tidKnowns <- unique(tidKnowns)
    
    tmp <- data.frame(protein=paralog,
               n_apms=length(apmsParaPreys), n_apms_known=length(apmsKnowns),
               n_tid=length(tidParaPreys), n_tid_known=length(tidKnowns)
               )
    collapsedDf[[paralog]] <- tmp
}
collapsedDf <- do.call(rbind, collapsedDf)
collapsedDf$pct_apms <- collapsedDf$n_apms_known/collapsedDf$n_apms * 100
collapsedDf$pct_tid <- collapsedDf$n_tid_known/collapsedDf$n_tid * 100

collapsedDf$protein <- factor(
  collapsedDf$protein,
  levels = names(annotation_colors$paralog)
)

collapsedDf$novel_apms <- collapsedDf$n_apms - collapsedDf$n_apms_known
collapsedDf$novel_tid <- collapsedDf$n_tid - collapsedDf$n_tid_known

out <- list()
for (protein in unique(collapsedDf$protein) ) {
  
  apms_novel <- collapsedDf$n_apms[collapsedDf$protein == protein] - collapsedDf$n_apms_known[collapsedDf$protein == protein]
  tid_novel <- collapsedDf$n_tid[collapsedDf$protein == protein] - collapsedDf$n_tid_known[collapsedDf$protein == protein]
  
  tmp1 <- rbind(
    data.frame(protein=protein,
             number=apms_novel,
             interactor='novel',
             known=collapsedDf$pct_apms[collapsedDf$protein == protein],
             method='AP-MS'),
    data.frame(protein=protein,
               number=collapsedDf$n_apms_known[collapsedDf$protein == protein],
               interactor='known',
               known=collapsedDf$pct_apms[collapsedDf$protein == protein],
               method='AP-MS') 
  )
  
  tmp2 <- rbind(
    data.frame(protein=protein,
               number=tid_novel,
               interactor='novel',
               known=collapsedDf$pct_tid[collapsedDf$protein == protein],
               method='TbID'),
    data.frame(protein=protein,
               number=collapsedDf$n_tid_known[collapsedDf$protein == protein],
               interactor='known',
               known=collapsedDf$pct_tid[collapsedDf$protein == protein],
               method='TbID') 
  )
  out[[protein]] <- rbind(tmp1, tmp2)
}
#
collapsedDf <- do.call(rbind, out)
collapsedDf$protein <- factor(collapsedDf$protein,
                              levels = names(annotation_colors$paralog)
                              )

outList <- list()
for (protein in unique(collapsedDf$protein) ) {
  nApms <- sum(collapsedDf$number[collapsedDf$protein==protein & collapsedDf$method=='AP-MS'])
  nTid <- sum(collapsedDf$number[collapsedDf$protein==protein & collapsedDf$method=='TbID'])
  
  pApms <- collapsedDf$known[collapsedDf$protein==protein &
                                collapsedDf$method=='AP-MS' &
                                collapsedDf$interactor=='known']
  pTid <- collapsedDf$known[collapsedDf$protein==protein &
                             collapsedDf$method=='TbID' &
                             collapsedDf$interactor=='known']
  
  tmpDf <- rbind(data.frame(protein=protein,
             preys=nApms,
             known=pApms, method='AP-MS'
             ),
        data.frame(protein=protein,
                   preys=nTid,
                   known=pTid, method='TbID'
        ) )
  outList[[protein]] <- tmpDf
}


outDf <- do.call(rbind, outList)
coeff <- max(outDf$known)/max(outDf$preys)
outDf$protein <- factor(
  outDf$protein, levels = names(annotation_colors$paralog)
)
outDf$node <- apmsDesign$node[match(outDf$protein, apmsDesign$paralog) ]
outDf$Label <- outDf$protein

ggplot(outDf, aes(x=protein)) +
  geom_bar(aes(y=preys, fill=protein), stat = 'identity', position = 'dodge') +
  geom_point( aes(y= known/coeff), color='#fec44f', group=1, size=2 ) +
  geom_line( aes(y= known/coeff), color='#fec44f', size=0.5,
             group=1, 
             show.legend = FALSE) +
  scale_y_continuous(
    # Features of the first axis
    name = "# of preys",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name="% known")
  ) + 
  #scale_fill_brewer(palette = 'Paired') +
  xlab('') +
  theme_cowplot() + 
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.title.y.right = element_text(color = '#fec44f'),
    axis.text.y.right = element_text(colour = '#fec44f')
  ) +
  scale_fill_manual(label=names(annotation_colors$paralog),
                    breaks = names(annotation_colors$paralog),
                    values = annotation_colors$paralog) +
  facet_wrap(~method, scales = 'free_x', ncol = 2)
ggsave(file.path(figureDir(), 'Figure1/stacked_proteins_collapsed_known_2ndaxis.pdf'),
       width = 6,
       height = 2.5
)
#
write.table(outDf,
  file.path(figureDir(), 'Figure1/data/stacked_proteins_collapsed_known_2ndaxis.tsv'),
  row.names = F, sep = '\t', quote = F
  )


#### Check which interactors are known/novel for 
# i, ii) preys specific to methods
# iii) preys significant in boith methods

toIt <- list('AP-MS' = setdiff(unlist(apmsPreys), unlist(tidPreys)),
      'Both' = intersect(unlist(tidPreys) , unlist(apmsPreys)),
      'TbID-MS' = setdiff(unlist(tidPreys), unlist(apmsPreys)) )

library(eulerr)
tmp <- list('AP-MS' = unique(unlist(apmsPreys)), 
            'TbID-MS' = unique(unlist(tidPreys)) )
fit <- euler(tmp)
plot(fit, quantities=T, fills=c('#f05455', '#0c607f') )

out <- list()
for (elem in names(toIt)) {
  pois <- toIt[[elem]]
  #
  olInteractome <- interactome[interactome$prey%in%pois, ]
  
  pcts <- c()
  ns <- c()
  methods <- c()
  nodes <- c()
  baits <- c()
  paras <- c()
  knownInteractors <- novelPreys <- list()
  for (bait in unique(olInteractome$bait )) {
    para <- unique(apmsDesign$paralog[apmsDesign$groups==bait])
    node <- unique(apmsDesign$node[apmsDesign$groups==bait])
    
    geneName <- ''
    
    if ( para %in% c('KRAS', 'NRAS', 'HRAS', 'ARAF', 'BRAF', 'RAF1' ) ) {
      geneName <- para
    } else if (para == 'MEK1') {
      geneName <- 'MAP2K1'
    } else if (para == 'MEK2') {
      geneName <- 'MAP2K2'
    } else if (para == 'ERK1') {
      geneName <- 'MAPK3'
    } else if (para == 'ERK2') {
      geneName <- 'MAPK1'
    }
    
    preys <- unique(gsub(';.*', '', olInteractome$prey[olInteractome$bait==bait] ))
    
    paras <- c(paras, c(para))
    nodes <- c(nodes, c(node) )
    baits <- c(baits, c(bait) )
    
    interactors <- names(neighbors(ppi, geneName))
    
    cmbKnown <- length(intersect(interactors, preys))
    
    nCmb <- length(preys)
    ns <-c(ns, c(nCmb))
    methods <- c(methods, c( 'Combined' ))
    
    cmbPctKnown <- cmbKnown/length(preys)
    pcts <- c(pcts, c(cmbPctKnown) )
  }
  ###
  df <- data.frame(bait=baits, paralog=paras, node=nodes,method=elem, preys=ns, known=pcts)
  tmp <- df
  tmp$nknown <- tmp$preys * tmp$known
  tmp$new <- tmp$preys - tmp$nknown
  
  x <- tmp[, -which( names(tmp ) == 'nknown' )]
  names(x)[which( names(x) == 'new' )] <- 'number'
  x$interactor <- 'novel'
  y <- tmp[, -which( names(tmp ) == 'new' )]
  names(y)[which( names(y) == 'nknown' )] <- 'number'
  y$interactor <- 'known'
  
  plotDf <- rbind(x, y)
  plotDf <- plotDf[order(plotDf$preys, decreasing = T), ]
  #ordering <- plotDf[plotDf$method==elem, 'bait']
  #ordering <- ordering[!duplicated(ordering)]
  #plotDf$bait <- factor(plotDf$bait, levels = ordering)

  plotDf$known <- plotDf$known * 100

  sumPctKnownCmb <- sum(plotDf$preys[plotDf$method==elem] * 
                          plotDf$known[plotDf$method==elem]) / sum(plotDf$preys[plotDf$method==elem])
  sumPctNovelCmb <- 100 - sumPctKnownCmb
  
  summaryPlotData <- data.frame(
    method = c(elem, elem),
    pct = c(
      sumPctKnownCmb,
      sumPctNovelCmb
    ),
    prey = c('known', 'novel' )
  )
  out[[elem]] <- summaryPlotData
}

summaryPlotData <- do.call(rbind, out)
###
ggplot(summaryPlotData, aes(x=method, y=pct, fill=prey) ) +
  geom_bar(position = 'stack', stat = 'identity') + theme_cowplot() + 
  scale_fill_manual(values=c(known='#A6CEE3', novel='#1F78B4') ) +
  ylab('% preys') + xlab('')

### AP-MS vs. TurboID
# what's the fraction of significant preys in either / both of them?
# head(interactome )
# 
# out <- list( )
# for (bait in unique(interactome$bait ) ) {
#   tmp <- interactome[interactome$bait == bait, ]
#   nPreys <- nrow(tmp)
#   
#   df <- data.frame(table(tmp$significant) )
#   names(df) <- c('significant', 'Count' )
#   
#   df$Fraction <- df$Count/nPreys
#   
#   df$bait <- bait
#   
#   out[[bait]] <- df
#   
# }
#
### all baits
# numDf <- do.call(rbind, out)
# numDf$Fraction <- numDf$Fraction * 100
# idxs <- match(numDf$bait, apmsDesign$groups)
# numDf$node <- apmsDesign$node[idxs]
# 
# numDf$bait <- factor(numDf$bait, levels = baitOrder)
# 
# ggplot(numDf, aes(y=bait, x=Fraction, fill=significant)) +
#   geom_bar(position = 'stack', stat = 'identity') +
#   scale_fill_manual(values=c('AP-MS'='pink', TbID='skyblue', shared='#b2df8a', core='#33a02c') ) +
#   labs(x = "% of preys") +
#   facet_wrap(~node, scales = 'free_y') +
#   theme_cowplot()
# ggsave(
#   file.path(figureDir(), 'S1/allbaits_significant_in_method.pdf' ),
#   width = 7.25, height = 5
# )
# write.table(
#   numDf,
#   file.path(figureDir(), 'S1/data/allbaits_significant_in_method.tsv' ),
#   quote = F, sep = '\t'
# )
# 
# avgFractions <- data.frame(table(interactome$significant)/nrow(interactome))
# names(avgFractions) <- c('significant', 'Fraction' )
# avgFractions$dummy <- 'Interactome'
# avgFractions$Fraction <- 100 * avgFractions$Fraction
# # 
# ggplot(avgFractions, aes(x=dummy, y=Fraction, fill=significant)) +
#   geom_bar(position = 'stack', stat = 'identity') +
#   scale_fill_manual(values=c('AP-MS'='pink', TbID='skyblue', shared='#b2df8a', core='#33a02c') ) +
#   labs(y = "% of preys") + xlab('') +
#   theme_cowplot()
# 
# ggsave(
#   file.path(figureDir(), 'S1/avg_pct_significant_in_method.pdf' ),
#   width = 3.5, height = 4.5
# )
# 
# write.table(avgFractions,
#             file.path(figureDir(), 'S1/data/avg_pct_significant_in_method.tsv' ),
#             row.names = F, sep = '\t', quote = F )