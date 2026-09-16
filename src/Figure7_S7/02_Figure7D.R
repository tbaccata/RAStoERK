source(here::here("src/utils/utils_path.R"))

wnt <- read.table(file.path(dataDir(), 
                            'integration/wnt_reac_hsa195721/string_reac_wnt_protein_annotations.tsv'),
                  sep = '\t', header = T)

wntInteractome <- interactome[gsub(';.*', '', interactome$prey) %in%wnt$node, ]

wntGenes <- unique(wntInteractome$prey)

wntGenes <- c(wntGenes,
              c('LZTS2', 'UBR5', 'USP9X', 'SIRT1', 'TCF3', 'USP47', 'GNA11')
)
wntPhsicalString <- stringNetworkFromAPI(wntGenes)

wntPhsicalString <- createPieChartNetork(wntPhsicalString, wntPhsicalString, wntInteractome, wntGenes)

allapms <- unique(unlist(apmsPreys) )
alltbid <- unique(unlist(tidPreys) )

V(wntPhsicalString)$method <- ifelse(
  V(wntPhsicalString)$name%in% intersect(allapms, alltbid),
  'both',
  ifelse(
    V(wntPhsicalString)$name%in% allapms, 'AP-MS', 'TbID'
  )
)

write_graph(wntPhsicalString, 
            file = file.path(figureDir(),
                             'processed/cytoscape/graphml/wnt_physical.graphml' ),
            format = "graphml")

write.table(wntInteractome,
            file = file.path(
              figureDir(),
              paste0(
                'Figure7/data/wnt_interactome_reactome.tsv'
              )
            ), sep = '\t', quote = F, row.names = F
)
