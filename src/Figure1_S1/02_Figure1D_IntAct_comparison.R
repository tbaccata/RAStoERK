#' @author Sebastian Didusch
#' @date 31.09.03
#' 
#' Compare other AP-MS and BioID data sets with our data

library(ggpubr)
library(ggbreak)
library(ggrepel)
library(eulerr)
### newest version of ggrepel plots labels twice when using scale_x_break
### remotes::install_github('slowkow/ggrepel@e94776b0a75b5c8f929c2e6ebb4c897135f6510a')
### download verison 0.9.5

source(here::here("src/utils/utils_path.R"))

core <- c('KRAS', 'NRAS', 'HRAS',
          'ARAF', 'BRAF', 'RAF1',
          'MAP2K1', 'MAP2K2',
          'MAPK3', 'MAPK1')

#### 
p <- (file.path(dataDir(), 'integration/apms_pdl_datasets') )
allKnownIntact <- read.table(file.path(p, "intact_core_pathway_as_bait.dat"),
                             header = T, sep = '\t')


#View(allKnownIntact)
no2hybridIntact <- allKnownIntact[grep('two hybrid', allKnownIntact$method, invert = T ), ]
refs <- (unlist(strsplit(no2hybridIntact$reference[no2hybridIntact$expansion_methd == 'psi-mi:MI:1060(spoke expansion)'], '\\|')) )
pmids <- grep('pubmed', (refs), value = T)
tab <- sort(table(pmids) )

pmidsOi <- names(tab[tab>25] )

outList <- list()
for (pmid in pmidsOi) {
  dat <- no2hybridIntact[grep(pmid, no2hybridIntact$reference), c('gene_a', 'gene_b') ]
  #
  
  dataSet <- names(
    tail(
      sort(table(no2hybridIntact$first_author[grep(pmid, no2hybridIntact$reference)])), 
      1))
  
  names(dat) <- c('bait', 'prey')
  dat <- dat[dat$bait%in%core,]
  dat$data <- dataSet
  
  outList[[dataSet]] <- dat
}

allKnownIntact$dupl <- paste0(allKnownIntact$gene_a, '_', allKnownIntact$gene_b)

allIntactAsBait <- allKnownIntact[!duplicated(allKnownIntact$dupl), c('gene_a', 'gene_b') ]
names(allIntactAsBait) <- c('bait', 'prey')
#
table(allIntactAsBait$bait)

allIntactAsBait$data <- 'IntAct'

outList[['IntAct_2025']] <- allIntactAsBait

for (type in c('Inactive', 'WT', 'Active')) {
  tmpDesign <- apmsDesign[apmsDesign$type==type & 
                            apmsDesign$batch==1, ]
  apmsOutlist <- tidOutlist <- list()
  
  paralogs <- unique(apmsDesign$paralog[apmsDesign$type==type & 
                                          apmsDesign$batch==1] )
  
  for (protein in paralogs) {
    pname <- protein
    if (protein == 'MEK1') {
      pname <- 'MAP2K1'
    } else if (protein == 'MEK2') {
      pname <- 'MAP2K2'
    } else if (protein == 'ERK1') {
      pname <- 'MAPK3'
    } else if (protein == 'ERK2') {
      pname <- 'MAPK1'
    }
    
    aPs <- unique(unlist(apmsPreys[ tmpDesign$groups[tmpDesign$paralog==protein]   ]))
    tPs <- unique(unlist(tidPreys[tmpDesign$groups[tmpDesign$paralog==protein] ] ))
    apmsOutlist[[pname]] <- data.frame(bait=pname, 
                                       prey=aPs )
    tidOutlist[[pname]] <- data.frame(bait=pname,
                                      prey=tPs )
  }
  apmsTypeDf <- do.call(rbind, apmsOutlist)
  tidTypeDf <- do.call(rbind, tidOutlist)
  #
  apmsTypeDf$data <- paste0('AP-MS_',type)
  tidTypeDf$data <- paste0('TbID_',type)
  
  outList[[paste0('AP-MS_', type )]] <- apmsTypeDf
  outList[[paste0('TbID_', type )]] <- tidTypeDf
}

dfAllApms <- dfAllTbid <- dfAll <- data.frame(bait='0', prey='0',data='0' )
for (protein in unique(paralogs)) {
  
  tmpDesign <- apmsDesign[apmsDesign$batch == 1, ]
  
  pname <- protein
  if (protein == 'MEK1') {
    pname <- 'MAP2K1'
  } else if (protein == 'MEK2') {
    pname <- 'MAP2K2'
  } else if (protein == 'ERK1') {
    pname <- 'MAPK3'
  } else if (protein == 'ERK2') {
    pname <- 'MAPK1'
  }
  
  aPs <- unique(unlist(apmsPreys[ tmpDesign$groups[tmpDesign$paralog==protein]   ]))
  tPs <- unique(unlist(tidPreys[tmpDesign$groups[tmpDesign$paralog==protein] ] ))
  dfAllApms <- rbind(dfAllApms, data.frame(bait=pname, 
                                           prey=aPs, data='AP-MS_all' ) )
  dfAllTbid <- rbind(dfAllTbid, data.frame(bait=pname,
                                           prey=tPs, data='TbID_all' ) )
  dfAll <- rbind(dfAll,
                 data.frame(bait=pname,
                            prey=unique(c(aPs,tPs)), data='All' )
  )
}
dfAllApms <- dfAllApms[dfAllApms$bait != '0', ]
dfAllTbid <- dfAllTbid[dfAllTbid$bait != '0', ]
dfAll <- dfAll[dfAll$bait != '0', ]
#
outList[['TbID_allbaits']] <- dfAllTbid
outList[['AP-MS_allbaits']] <- dfAllApms
outList[['Both_allbaits']] <- dfAll

df <- do.call(rbind, outList)


write.table(
  df,
  file.path(figureDir(),'Figure1/data/networks_of_other_interactomes.tsv'),
  row.names = F, sep = '\t', quote = F
)

####
nbaits <- unlist(lapply(outList, function(x) length(unique(x$bait)) ))
npreys <- unlist(lapply(outList, function(x) length(unique(x$prey))  ))


df <- data.frame(nbaits, npreys )
df$dataSet <- rownames(df)
df$method <- ''

for (first_author in unique(grep('et',df$dataSet,value = T) )) {
  meth <- tail(
    names(
      sort(
        table(
          no2hybridIntact$method[no2hybridIntact$first_author==first_author ] )
      ) 
    ), 
    1)
  df$method[df$dataSet==first_author] <- meth
}

df$method[grep('TbID', df$dataSet)] <- 'TbID'
df$method[grep('AP-MS', df$dataSet)] <- 'AP-MS'
df$method[grep('Both', df$dataSet, ignore.case = F)] <- 'Combined'

df$nbaits[grep('Active', df$dataSet )] <- 15
df$nbaits[grep('Both|all', df$dataSet, ignore.case = T )] <- 35

df$method[grep('coimmunoprecipitation|pull', df$method)] <- 'AP-MS'
df$method[grep('proximity-dependent', df$method)] <- 'BioID'



x <- no2hybridIntact[no2hybridIntact$first_author%in%df$dataSet, c('first_author', 'reference')]
author2pubmed <- x[!duplicated(x$first_author),]

replDf <- data.frame(
  first_author = c('Cho NH. et al.(2022)', 'Kennedy SA et al.(2019)',
                   'Buljan M. et al.(2020)', 'Swaney DL. et al.(2021)', 'Adhikari H. et al.(2018)',
                   'Hein MY. et al.(2015)', 'Huttlin EL. et al.(2021)', 'So J. et al.(2015)',
                   'Kim M. et al.(2021)', 'Luis F. Iglesias-Martinez et al.(2023)' ),
  label = c('OpenCell', 'EGFR network', 'Buljan et al.', 'Swaney et al.', 'Adhikari et al.',
            'Hein et al.', 'BioPlex3.0', 'So et al.',
            'Kim et al.', 'Iglesias-Martinez et al.')
)

df$label <- df$dataSet
df$reference <- ''
for (fa in replDf$first_author) {
  df$label[df$dataSet==fa] <- replDf$label[replDf$first_author==fa]
  df$reference[df$dataSet==fa] <- author2pubmed$reference[author2pubmed$first_author==fa]
}  
#plotdf <- df[df$nbaits>1 & df$dataSet != 'Both_allbaits',]
plotdf <- df # df[df$nbaits>1,]

plotdf <- plotdf[grep('AP-MS|TbID', plotdf$dataSet, invert = T ), ]
plotdf <- plotdf[plotdf$dataSet != 'IntAct_2025',]
plotdf$label[plotdf$label=='Both_allbaits'] <- 'all'

size_range <- range(plotdf$npreys, na.rm = TRUE)


ggplot(plotdf, 
       aes(
         x = nbaits,
         y = npreys,
         size = npreys,
         color = method,
         label = label
       )) + geom_jitter(width = 0.25) +
  # geom_jitter(data = plotdf[grep('all|AP-MS|TbID', plotdf$dataSet),], width = 0.5) + 
  # geom_jitter(data = plotdf[grep('all|AP-MS|TbID', plotdf$dataSet, invert = T),],
  #               alpha=0.6, width = 0.25) + 
  scale_color_manual(values = c('AP-MS'='#f05455ff',
                                BioID='#8E9DCC',
                                TbID='#0c607eff',
                                Combined='#878E88',
                                Y2H='#B4F5A6'
  )) +
  scale_size_area(limits = size_range, max_size = 12,
                  breaks=c(30, 500, 1500, 2500)
  ) +
  #scale_size_area(max_size = 12) +   # scales by AREA, not radius
  
  #scale_y_break(c(1620, 2520), scales = 0.1) +
  scale_x_continuous(breaks= c(1,3,4,5,6,7,35), limits = c(0.8,38) ) +
  scale_x_break(c(8, 33), scales = 0.1) +
  geom_text_repel(show.legend = FALSE, color = "black", alpha=1, size = 3) +
  xlab("# of baits") +
  ylab("# of preys") +
  theme_cowplot()

ggsave(
  file.path(figureDir(), 'Figure1/comparison_with_other_interactomes.pdf'),
  width = 8, height = 6
)
write.table(
  plotdf,
  file.path(figureDir(),'Figure1/data/comparison_with_other_interactomes.tsv'),
  row.names = F, sep = '\t', quote = F
)

#### Focus on ERK interactome
erkPreys <- unique(c(apmsPreys$ERK2WT, apmsPreys$ERK1WT,
                     tidPreys$ERK1WT, tidPreys$ERK2WT) )
erksIntact <- allKnownIntact[allKnownIntact$gene_a %in% c('MAPK3', 'MAPK1'), ]
erksIntact <- erksIntact[grep('two hybrid', erksIntact$method, invert = T),]
# only high through-put
erksIntact <- erksIntact[erksIntact$expansion_methd == 'psi-mi:MI:1060(spoke expansion)', ]

intactUniqueApmsErkPreys <- unique(erksIntact$gene_b)
interactors <- list('in IntAct'=intactUniqueApmsErkPreys,
                    'this study'=erkPreys  )
erkInteractors <- data.frame('ERK_prey_protein'=unique(unlist(interactors)) )

erkInteractors$`in IntAct` <- ifelse(erkInteractors$ERK_prey_protein %in% intactUniqueApmsErkPreys,1,0)
erkInteractors$`this study` <- ifelse(erkInteractors$ERK_prey_protein %in% erkPreys,1,0)
erkInteractors$`ERK relationship in Signor` <- ifelse(erkInteractors$ERK_prey_protein %in% 
                                                        erkPreys,1,0)
G_signor <- graph_from_data_frame(signorNetwork, directed = F)
erkSignalingInteractors <- unique(c(names(neighbors(G_signor, 'MAPK3')),
  names(neighbors(G_signor, 'MAPK1'))))

interactors[['in Signor']] <- erkSignalingInteractors
### Venn diagram
pdf(file.path(figureDir(), 'S1/venn_signor_erkpreys_intact_erk_as_baits.pdf' ),
    width = 5, height = 5)
plot(venn(interactors),
     legend = T,
     fills = c('in IntAct' = '#D9E5D6',
               'this study' = '#c2649b',
               'in Signor' = '#41AB5D'
               ),
     quantities = T,
     main = 'ERK relationship')
dev.off()
###

erkInteractors$`ERK relationship in Signor` <- ifelse(erkInteractors$ERK_prey_protein %in% 
                                                        erkSignalingInteractors,1,0)
#
write.table(
  erkInteractors,
  file.path(figureDir(),'Figure1/data/euler_erkpreys_intact_erk_as_baits.tsv'),
  row.names = F, sep = '\t', quote = F
)

df <- do.call(rbind, outList)
### opencel  Cho NH. et al.(2022)
### bioplex Huttlin EL. et al.(2021)
### kolch Kennedy SA et al.(2019)
### Buljan M. et al.(2020)

toIt <- list(
  OpenCell='Cho NH. et al.(2022)',
  BioPlex='Huttlin EL. et al.(2021)',
  Kolch='Kennedy SA et al.(2019)',
  Gstaiger='Buljan M. et al.(2020)'
)

erkInteractomes <- df[df$bait%in%c('MAPK1', 'MAPK3'), ]
for (elem in names(toIt)) {
  preys <- erkInteractomes$prey[erkInteractomes$data == toIt[[elem]] ]
  erkInteractors[[elem ]] <- ifelse(erkInteractors$ERK_prey_protein %in% preys, 
                                    1, 0 )
  signorIntactPreys <- preys[preys %in% erkSignalingInteractors]
  print(paste(setdiff(signorIntactPreys, erkInteractors$ERK_prey_protein ), collapse = ',' ))
}

write.table(
  erkInteractors,
  file.path(figureDir(),'Figure1/data/euler_erkpreys_intact_erk_as_baits.tsv'),
  row.names = F, sep = '\t', quote = F
)


### Core pathway network
# get high-throughput studies
htApms <- no2hybridIntact[no2hybridIntact$expansion_methd=='psi-mi:MI:1060(spoke expansion)', ]
# get core pathway
htApmsCore <- htApms[htApms$gene_b %in% coreNames, ]

g <- make_empty_graph(directed = TRUE)
all_nodes <- unique(c(htApmsCore$gene_a, htApmsCore$gene_b))
g <- add_vertices(g, length(all_nodes), name = all_nodes)
for (elem in unique(htApmsCore$gene_a ) ) {
  tmp <- htApmsCore[htApmsCore$gene_a==elem, ]
  
  tab <- table(tmp$gene_b)
  for (idx in 1:length(tab)) {
    g <- add_edges(g, c(elem, names(tab[idx]) ),
                   attr=list(weight= tab[idx]) )
  }
}

write.table(as_data_frame(g),
            file.path(figureDir(),
                      paste0('Figure1/core_pathway_interactions/core_pathway_intact.tsv')
                      ),
            row.names = F, sep = '\t', quote = F)


tmp <- as.data.frame(table(htApmsCore$gene_a))
names(tmp) <- c('bait', 'degree')

write.table(htApmsCore,
            file.path(figureDir(),
                      paste0('Figure1/core_pathway_interactions/core_pathway_intact_mitab.tsv')
            ),
            row.names = F, sep = '\t', quote = F)
