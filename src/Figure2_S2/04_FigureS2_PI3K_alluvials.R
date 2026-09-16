source(here::here("src/utils/utils_path.R"))


rasAlluvial <- read.table(file.path(figureDir(), 
            'Figure2/data/alluvial/ras_flows_per_prey.txt'),
  sep = '\t', header = T
)
#
rafAlluvial <- read.table(file.path(figureDir(), 
                                     'Figure2/data/alluvial/raf_flows_per_prey.txt'),
                           sep = '\t', header = T
)
#
mekAlluvial <- read.table(file.path(figureDir(), 
                                     'Figure2/data/alluvial/mek_flows_per_prey.txt'),
                           sep = '\t', header = T
)
#
erkAlluvial <- read.table(file.path(figureDir(), 
                                     'Figure2/data/alluvial/erk_flows_per_prey.txt'),
                           sep = '\t', header = T
)

allPreys <- unique(c(
  rasAlluvial$GENE.NAMES,
  rafAlluvial$GENE.NAMES,
  mekAlluvial$GENE.NAMES,
  erkAlluvial$GENE.NAMES
))

### PI3K 
library(gprofiler2)
gres <- gost(allPreys,
     sources = c('GO:BP', 'KEGG', 'REAC'),
     significant = F, evcodes = T
     )

gostres <- gres$result


pi3kPreys <- unlist(strsplit(gostres$intersection[gostres$term_name=='PI3K-Akt signaling pathway'], 
                ',') )

rasAlluvial <- rasAlluvial[rasAlluvial$GENE.NAMES%in%pi3kPreys, ]
rafAlluvial <- rafAlluvial[rafAlluvial$GENE.NAMES%in%pi3kPreys, ]
mekAlluvial <- mekAlluvial[mekAlluvial$GENE.NAMES%in%pi3kPreys, ]
erkAlluvial <- erkAlluvial[erkAlluvial$GENE.NAMES%in%pi3kPreys, ]

toIt <- list(ras=rasAlluvial,
             raf=rafAlluvial,
             mek=mekAlluvial,
             erk=erkAlluvial
             )

alluvColors <- list(
  ras = c(
    common = "#00b1e8",
    shared = "#34bcbd",
    KRAS = '#f26924',
    NRAS = "#f68b22",
    HRAS = "#faa61e",
    missing = 'grey80'
  ),
  raf = c(
    common = "#00b1e8",
    shared = "#34bcbd",
    ARAF = "#174a92",
    BRAF = "#4f88c7",
    RAF1 = "#a1bce2",
    missing = 'grey80'
  ),
  mek = c(
    common = "#00b1e8",
    shared = "#34bcbd",
    MEK1 = "#009147",
    MEK2 = "#88be40",
    missing = 'grey80'
  ),
  erk = c(
    common = "#00b1e8",
    shared = "#34bcbd",
    ERK1 = "#a62877", 
    ERK2 = "#dfa0c2",
    missing = 'grey80'
  )
)

dir.create(file.path(
  figureDir(), 
  'S2/pi3k/'
), showWarnings = F)


plotDatas <- list()
for (elem in names(toIt) ) {
  alluv <- toIt[[elem]]
  #
  alluvPlotData <- as.data.frame(table(alluv$set_WT,
                                      alluv$set_Inactive,
                                      alluv$set_Active
  ))
  names(alluvPlotData) <- c('WT', 'Inactive', 'Active', 'Freq')
  alluvPlotData <- prepareAlluvialPlot(alluvPlotData, toupper(elem), placeholder = T, gapSize = 5)
  alluvPlotData <- alluvPlotData[alluvPlotData$Freq > 0, ]
  
  
  ggplot(alluvPlotData,
         aes(y = Freq,
             axis1 = Inactive, axis2 = WT, axis3 = Active)) +
    geom_alluvium(aes(fill = WT), curve_type='cubic',
                  width = 1/32, knot.pos = 0, reverse = FALSE) +
    scale_fill_manual(values = alluvColors[[elem]] ) +
    guides(fill = "none") +
    geom_stratum(width = 1/32, reverse = FALSE) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              reverse = FALSE) + ylab('# of preys') +
    scale_x_continuous(breaks = 1:3, labels = c("Inactive", "WT", "Active")) +
    theme_cowplot() #+ coord_flip()
  ggsave(file.path(figureDir(), 
                   paste0('S2/pi3k/S2A_', elem, '_alluvial.pdf') ), 
         width = 6, height = 5.5)
  
}


### CA overlap
caAlluvial <- read.table(
  file.path(
    figureDir(),
    'Figure2/data/alluvial/rasrafmek_active_flows_per_prey.txt'),
  sep = '\t',
  header = T
)

caAlluvial <- caAlluvial[caAlluvial$GENE.NAMES%in%pi3kPreys, ]
allCats <- as.data.frame(table(caAlluvial$set_Active_ras,
                               caAlluvial$set_Active_raf,
                               caAlluvial$set_Active_mek
))
names(allCats) <- c('RAS', 'RAF', 'MEK', 'Freq')
allCats <- allCats[allCats$Freq > 0, ]

tmp <- data.frame(RAS=c('placeholder1', 'placeholder2','placeholder3', 'placeholder4', 'placeholder5'),
                  RAF=c('placeholder1a', 'placeholder2a', 'placeholder3a', 'placeholder4a', 'placeholder5a'),
                  MEK=c('placeholder1b','placeholder2b', 'placeholder3b', 'trash1', 'trash2'),
                  Freq=c(5, 5, 5, 5, 5)
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
                 'S2/pi3k/2B_rasrafmek_ca_pi3k_alluvial.pdf'), 
       width = 6, height = 5.5)
