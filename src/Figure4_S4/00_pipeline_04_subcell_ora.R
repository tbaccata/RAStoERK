#' @author Sebastian Didusch
#' @date 25.09.15
#' ORA stratified by orgenelle IP subcell. localization
#' 
#' Results show paralog-specific interactors at the resolution
#' of compartments
#' 

source(here::here("src/utils/utils_path.R"))

endoLysoGolgi <- c('recycling_endosome', 'early_endosome', 'lysosome', 'trans-Golgi', 'Golgi' )

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


bait2prey2organelle <- data.frame(Bait='0', Prey='0', Organelle='0')
oraInputList <- list()
for (loc in locs) {
  organelle <- loc
  
  if (loc == 'endoLysoGolgi' ) {
    organelle <- endoLysoGolgi
  }
  
  orgInteractome <- interactome[interactome$subcell_localization%in%organelle,]
  proteins <- unique(apmsDesign$paralog[apmsDesign$type !='Control'])
  
  for ( protein in proteins) {
    orgProtInteractome <- orgInteractome[orgInteractome$bait%in%apmsDesign$groups[apmsDesign$paralog==protein], ]
    preys <- unique(orgProtInteractome$prey)
    
    
    
    if (length(preys) > 3) {
      oraInputList[[paste0(loc, '__', protein ) ]] <- preys
    }
    
    if (length(preys) > 0) {
      bait2prey2organelle <- rbind(
        bait2prey2organelle,
        data.frame(
          Bait=protein,
          Prey=preys,
          Organelle=loc
        )
      )
    }
  }
}
#

library(gprofiler2)
resGprofiler <- gost(oraInputList, 
                     sources = c('GO:BP', 'REAC'),
                     evcodes = T )

stratByOrganelles <- resGprofiler$result

stratByOrganelles$subcell_localization <- sapply(strsplit(stratByOrganelles$query, '__' ), 
                                                 function(x) x[1] )
stratByOrganelles$protein <- sapply(strsplit(stratByOrganelles$query, '__' ), 
                                    function(x) x[2] )

stratByOrganelles$term_id <- paste0(stratByOrganelles$subcell_localization,
                                    '_', 
                                    stratByOrganelles$term_id)

stratByOrganelles$parents <- NULL

write.table(
  stratByOrganelles,
  file.path(figureDir(), 'Figure4/data/ora_stratified_by_subcell_organelle.tsv'),
  row.names = F, sep = '\t', quote = F
)
