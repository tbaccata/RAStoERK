#' @author Sebastian Didusch
#' @date 25.09.15
#' ORA stratified by subcell. loc
#' 

source(here::here("src/utils/utils_path.R"))

stratByOrganelles <- read.table(
     file.path(figureDir(), 'Figure3/data/ora_stratified_by_subcell_organelle.tsv'),
     header = T, sep = '\t'
)

tmp <- stratByOrganelles$term_id[order(stratByOrganelles$subcell_localization, 
                                       stratByOrganelles$p_value, decreasing = F)]
tmp <- tmp[!duplicated(tmp)]

stratByOrganelles$term_id <- factor(stratByOrganelles$term_id, levels = tmp )

locs <- c("plasma_membrane",
          "endoLysoGolgi",
          "actin_cytoskeleton",
          "ER",
          # "14-3-3_scaffold",
          "stress_granule",
          # 'mitochondrion',
          #"proteasome",
          'nucleus', 
          'unclassified')
stratByOrganelles$subcell_localization <- factor(
  stratByOrganelles$subcell_localization,
  levels = locs
)

#######
#######
data <- read.delim(
  file.path(dataDir(), 'processed/enrichment/250919_stringdb_selected_ora_terms.tsv'),
  header = T, sep = '\t'
)

allLevels <- list()

out <- list()
prev <- ''
newLevels <- c()

for (loc in locs) {
  dbs <- c('GO Process', 'Reactome')
  for (db in dbs) {
    plotdata <- stratByOrganelles[stratByOrganelles$subcell_localization==loc, ]
    pattern <- paste0(unique(data$term_id[data$category==db]), collapse = '|')
    plotdata <- plotdata[grep(pattern, plotdata$term_id ), ]
    #
    x <- levels(plotdata$term_id)
    lvls <- x[x%in%plotdata$term_id ]
    for (trmid in lvls) {
      newLevels <- c(newLevels, trmid)
      
      subsetted <- plotdata[plotdata$term_id == trmid, ]
      curr <- unlist(strsplit(trmid, '_'))[1]
      
      out[[trmid]] <- subsetted
      
      
    }
    df <- data.frame(matrix(nrow = 1, ncol = ncol(plotdata ) ) )
    names(df) <- names(plotdata)
    df$term_id <- paste0(loc, '__',db, '_spacer')
    df$term_name <- ''
    df$source <- unique(plotdata$source)
    
    df$protein <- 'KRAS'
    newLevels <- c(newLevels, df$term_id)
    
    out[[ paste0(loc, '__',db, '_spacer')]] <- df
  }
  
}

plotdata <- do.call(rbind,out)
plotdata$term_id <- factor(
  plotdata$term_id,
  levels = rev(newLevels)
)

patterns <- c("Regulation of cytoskeletal remodeling and cell spreading by IPP complex components",
              'Recruitment of mitotic centrosome proteins and complexes',
              'transcriptional', 'transcriptional', 'and', 'mediated', 'negative', 'regulation', 'inhibition', 
              'endoplasmic reticulum', 'plasma membrane')
replacements <- c('Reg. of cytoskeletal remodel. by IPP complex components', 
                  'Recr. of mit. centros. prot. & complexes',
                  'transl.', 'transcr.', '&', 'med.', 'neg.', 'reg.', 'inhib.', 'ER', 'PM')

for (idx in seq_along(patterns)) {
  plotdata$term_name <- gsub(patterns[idx],
                             replacements[idx],
                             plotdata$term_name
  )
}

#toChange <- 
#plotdata$term_name[plotdata$term_name == toChange] <- 'Reg. of cyto. remodeling by IPP complex components'

ann <- plotdata %>%
  distinct(term_id, subcell_localization, term_name)

# keep the same factor ordering as in plotdata (important if you reordered y)
if (is.factor(plotdata$term_id)) {
  ann$term_id <- factor(as.character(ann$term_id), levels = levels(plotdata$term_id))
} else {
  ann$term_id <- factor(ann$term_id, levels = unique(plotdata$term_id))
}

#
plotdata$protein <- factor(plotdata$protein,
                           levels = c( "KRAS","NRAS", "HRAS",  "ARAF", "BRAF",
                                       "RAF1", "MEK1", "MEK2",  "ERK1", "ERK2" )
)
##

loc2color <- c(plasma_membrane='#708097ff',
               endoLysoGolgi='#26a9e0ff',
               actin_cytoskeleton='#c76893ff',
               ER='#603813ff',
               stress_granule='#72be44ff',
               nucleus='#fab346ff',
               unclassified='grey70')

plotdata$p_value[!is.na(plotdata$p_value) &
                   plotdata$p_value < 1e-10 ] <- 1e-10

ggplot(plotdata, aes(
  y = term_id,
  x = protein,
  size = intersection_size,
  color = -log10(p_value)
)) +
  geom_point() +
  scale_color_viridis_c(option = "viridis") +
  scale_size(range = c(1, 6)) +
  
  # allow a new fill scale for the annotation tiles
  ggnewscale::new_scale_fill() +
  
  # draw one tile per term; don't inherit the global aes (so size won't apply)
  geom_tile(
    data = ann,
    inherit.aes = FALSE,
    aes(x = 0, y = term_id, fill = subcell_localization),
    width = 0.5, height = 0.9,
    color = NA   # remove tile border so stroke/size doesn't show
  ) +
  scale_fill_manual(values = loc2color, name = "Localization") +
  
  # give space on the left for the annotation strip
  scale_x_discrete(expand = expansion(add = c(1, 0))) +
  
  # pretty y labels (use ann to match one label per term)
  scale_y_discrete(breaks = ann$term_id, labels = ann$term_name) +
  
  coord_cartesian(clip = "off") +    # if you want the strip to sit in the margin
  theme_cowplot() + 
  xlab('') + ylab('') +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid.major.x = element_line(color = "grey90", size = 0.2),
    panel.grid.minor.x = element_line(color = "grey90", size = 0.1),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    panel.grid.minor.y = element_line(color = "grey90", size = 0.1)
  )

ggsave(
  file.path(figureDir(), 'Figure4/4A_ora_stratified_by_subcell_organelle.pdf'),
  useDingbats = FALSE,
  width = 8.25, height = 12 )

### --- Cytoscape STRING networks ###
reacOfInterest <- c('Asparagine N-linked glycosylation', 
                    'Axon guidance',
                    'Membrane Trafficking', 
                    'Cell Cycle'
                    )

reacSelectedData <- stratByOrganelles[stratByOrganelles$source=='REAC',]

allapms <- unique(unlist(apmsPreys) )
alltbid <- unique(unlist(tidPreys) )

dir.create(  file.path(dataDir(), 'processed/cytoscape'), showWarnings = F)
dir.create(  file.path(dataDir(), 'processed/cytoscape/graphml'), showWarnings = F)

for (reacTerm in reacOfInterest) {
  preysOfTerm <- reacSelectedData$intersection[reacSelectedData$term_name == reacTerm] %>%
    strsplit(',') %>% unlist() %>% unique()
  
  termInteractome <- interactome[gsub(';.*', '', interactome$prey) %in%preysOfTerm, ]
  
  termFunctionalString <- stringNetworkFromAPI(preysOfTerm, network_type = 'functional')
  termPhysicalString <- stringNetworkFromAPI(preysOfTerm, network_type = 'physical')
  termString <- createPieChartNetork(termFunctionalString, termPhysicalString, termInteractome, preysOfTerm)
  
  V(termString)$method <- ifelse(
    V(termString)$name%in% intersect(allapms, alltbid),
    'both',
    ifelse(
      V(termString)$name%in% allapms, 'AP-MS', 'TbID'
    )
  )
  
  tmp <- termInteractome[!duplicated(termInteractome$prey), ]
  #
  V(termString)$localization <- tmp$subcell_localization[match(V(termString)$name, tmp$prey) ]
  
  df <- as_data_frame(termString, what = 'vertices')
  write.table(df, file.path(dataDir(), paste0('processed/cytoscape/graphml/',
                                              make.names(reacTerm), '_pies.tsv' ) ),
              row.names = F, sep = '\t', quote = F )
  
  write_graph(termString,
              file = file.path(dataDir(),
                               paste0('processed/cytoscape/graphml/',
                                 make.names(reacTerm), '.graphml') ) , format = "graphml")
}

