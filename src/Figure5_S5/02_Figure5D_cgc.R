#' @author Sebastian Didusch
#' @date 25.10.02
#' 
#' map Cancer Gene Consensus data on interactome
#' What are the cancer modules seen by interactome?

source(here::here("src/utils/utils_path.R"))

caInteractome <- interactome[interactome$bait%in%apmsDesign$groups[apmsDesign$type=="Active"], ]

cgc <- read.table(
  file.path(dataDir(), 'integration/cancer_consensus_genes/Census_10_09_21_2025.tsv'),
  header = T, sep = '\t'
)

tierOneNonFusion <- cgc[cgc$Tier == 1 & cgc$Role.in.Cancer != 'fusion', ]
# only select tier 1 and oncogenes or tumor supressors
clinicalRelevance <- caInteractome[caInteractome$prey %in% tierOneNonFusion$Gene.Symbol, ]

idxs <- match(clinicalRelevance$prey, tierOneNonFusion$Gene.Symbol )
clinicalRelevance$RoleInCancer <- tierOneNonFusion$Role.in.Cancer[idxs ]

preysOfTerm <- unique(clinicalRelevance$prey)

G_funct <- stringNetworkFromAPI(preysOfTerm, network_type = 'functional')
G_funct <- delete_edges(G_funct, which(E(G_funct)$score <0.5))
G_term <- stringNetworkFromAPI(preysOfTerm, network_type = 'physical')

singetons <- names(which(degree(G_term)<1) )
for (prot in singetons) {
  x <- names(neighbors(G_funct, prot) )
  x <- paste(x, collapse = ',')
  print(paste(prot, x))
}

G_term <- createPieChartNetork(G_term, G_term, clinicalRelevance, preysOfTerm)
#
allapms <- unique(unlist(apmsPreys) )
alltbid <- unique(unlist(tidPreys) )
#
V(G_term)$method <- ifelse(
  V(G_term)$name %in% intersect(allapms, alltbid),
  'both',
  ifelse(
    V(G_term)$name %in% allapms,
    'AP-MS',
     'TbID'
  )
)

idxs <- match(V(G_term)$name, tierOneNonFusion$Gene.Symbol )
V(G_term)$`Role in Cancer` <- tierOneNonFusion$Role.in.Cancer[idxs ]

write_graph(G_term,
            file = file.path(
              dataDir(),
                'processed/cytoscape/graphml/cgc_preys.graphml'
            ),
            format = "graphml")

idxs <- match(clinicalRelevance$bait, apmsDesign$groups)
clinicalRelevance$BaitProtein <- apmsDesign$paralog[idxs]
clinicalRelevance$BaitProtein <- gsub('ERK1', 'MAPK3', clinicalRelevance$BaitProtein  )
clinicalRelevance$BaitProtein <- gsub('ERK2', 'MAPK1', clinicalRelevance$BaitProtein  )
clinicalRelevance$BaitProtein <- gsub('MEK1', 'MAP2K1', clinicalRelevance$BaitProtein  )
clinicalRelevance$BaitProtein <- gsub('MEK2', 'MAP2K2', clinicalRelevance$BaitProtein  )

clinicalRelevance$known <- ''
for (protein in unique(clinicalRelevance$BaitProtein)) {
  known <- names(neighbors(ppi, protein ))
  idxs <- which(clinicalRelevance$BaitProtein==protein)
  
  isKnown <- ifelse(clinicalRelevance$prey[idxs] %in% known, 'known', 'novel')
  clinicalRelevance$known[idxs] <- isKnown
}

table(clinicalRelevance$prey, clinicalRelevance$known)


clinicalRelevance$method <- ifelse(
  clinicalRelevance$prey %in% intersect(allapms, alltbid),
  'both', 
  ifelse( clinicalRelevance$prey %in% allapms, 'AP-MS', 'TbID' )
)
table(clinicalRelevance$method)


clinicalRelevance$key <- paste0(clinicalRelevance$BaitProtein, '_', clinicalRelevance$prey)

tmp <- clinicalRelevance[, c('BaitProtein', 'bait', 'prey', 'key', 'method', 'known', 'RoleInCancer'), ]

tmp <- tmp[!duplicated(tmp$key), ]

tbidPct <- as.data.frame(table(tmp$known[tmp$method!='AP-MS'])/nrow(tmp[tmp$method!='AP-MS',]) * 100 )
apmsPct <- as.data.frame(table(tmp$known[tmp$method!='TbID'])/nrow(tmp[tmp$method!='TbID',]) * 100 )
bothPct <- as.data.frame(table(tmp$known[tmp$method=='both'])/nrow(tmp[tmp$method=='both',]) * 100 )
#
tbidPct$method <- 'TbID'
apmsPct$method <- 'AP-MS'
bothPct$method <- 'both'

df <- rbind(bothPct, apmsPct)
df <- rbind(df, tbidPct)

names(df) <- c('prey', 'percentage', 'method')


### % known/novel over combined CA data set:
pctCmb <- table(tmp$known)/nrow(tmp) * 100
df <- rbind(df,
            data.frame(
              prey=c('known', 'novel'),
              percentage=c(pctCmb['known'], pctCmb['novel']),
              method='Combined'
            )
            )

df$method <- factor(df$method, levels = c('AP-MS', 'TbID', 'both', 'Combined'))

### % known/novel by method or those in both
ggplot(df, aes(x=method, y=percentage, fill=prey)) +
  geom_bar(stat = 'identity') +
  scale_fill_brewer(palette = 'Paired') +
  xlab('') + ylab('%') +
  theme_cowplot()
ggsave(file.path(figureDir(), 'S5/5D_CGC_preys_pct_known_novel.pdf'),
       width = 4.5, height = 3.1)

write.table(df,
            file.path(figureDir(), 'S5/data/CGC_preys_pct_known_novel.tsv'),
       row.names = F, sep = '\t', quote = F)

tmp$bait <- NULL
write.table(tmp,
            file.path(figureDir(), 'S5/data/CGC_preys_by_method_and_known.tsv'),
            row.names = F, sep = '\t', quote = F
            )
