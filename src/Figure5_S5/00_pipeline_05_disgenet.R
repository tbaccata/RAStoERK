#' @author Sebastian Didusch
#' @date 25.09.02
#' 
#' Download disease - gene associations from disgenet

## define paths
source(here::here("src/utils/utils_path.R"))

library(biomaRt)
library(disgenet2r)

API_KEY = 'TOP-SECRET-API-KEY'
Sys.setenv(DISGENET_API_KEY= API_KEY)

### complete proteome
# Connect to Ensembl BioMart (human dataset)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Retrieve HGNC symbols for protein-coding genes
protein_coding_genes <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype"),
  filters = "biotype",
  values = "protein_coding",
  mart = mart
)

# Just the unique HGNC gene symbols
hgnc_symbols <- unique(protein_coding_genes$hgnc_symbol)
hgnc_symbols <- hgnc_symbols[hgnc_symbols != ""]

###
bin_size <- 99
start <- 1
end <- bin_size
intervals <- data.frame(start=start, end=end)

while (end < length(hgnc_symbols)) {
  start <- end + 1
  end <- end + bin_size
  
  if (end >= length(hgnc_symbols)) end <- length(hgnc_symbols)
  
  intervals <- rbind(intervals,
                     data.frame(
                       start=start, end=end  )
  )
}

outList <- list()
for (idx in 1:nrow(intervals)) {
  
  intervals[idx,]$start
  results <- gene2disease(
    gene     = hgnc_symbols[intervals[idx,]$start:intervals[idx,]$end],
    database = "CURATED",
    score =c(0.3, 1)
  )
  outList[[idx]] <- results
  print(paste(idx, '/', nrow(intervals) ) )
}

dis2proteome <- do.call(rbind,
                      lapply(outList, function(x) x@qresult) )
x <-sapply(dis2proteome, class)
write.table(dis2proteome[, names(x[x!='list'])],
            file.path(dataDir(), 'processed/proteome_gene2disease_curated_score_gt_0.3.tsv' ),
            row.names = F, sep = '\t', quote = F )

# hcDisGenes <- unique(dis2proteome$gene_symbol[dis2proteome$score > 0.66] )
###
activeOverlap <- read.table(file.path(figureDir(), 'S2/data/combined_active_interactome_binary_matrix.tsv'),
                            sep = '\t', header = T )

oraInput <- list(RAS=activeOverlap$Gene.names[activeOverlap$RAS==1 ], 
                 RAF=activeOverlap$Gene.names[activeOverlap$RAF==1 ], 
                 MEK=activeOverlap$Gene.names[activeOverlap$MEK==1], 
                 ERK=activeOverlap$Gene.names[activeOverlap$ERK==1] )

oraInput[['RASRAFMEK']] <-activeOverlap$Gene.names[activeOverlap$RAS==1 &
                                                     activeOverlap$RAF==1 &
                                                     activeOverlap$MEK==1 ]
oraInput[['RASRAF']] <-activeOverlap$Gene.names[activeOverlap$RAS==1 &
                                                  activeOverlap$RAF==1 ]
oraInput[['RAFMEK']] <-activeOverlap$Gene.names[activeOverlap$RAF==1 &
                                                  activeOverlap$MEK==1 ]
out <- list()
for (elem in names(oraInput)) {
  prots <- oraInput[[elem]]
  
  res <- disgenet2r::disease_enrichment(prots,
                                        common_entities = 3)
  qres <- res@qresult
  qres$query <- elem
  out[[elem]] <- qres
  
  print(paste(elem, nrow(qres )) )
}


df <- do.call(rbind, out )

signDf <- df[df$p_adjusted < 0.05, ]

dis2proteinList <- list()
for (elem in unique(signDf$diseaseName)) {
  x <- unique(unlist(strsplit(signDf$intersection[signDf$diseaseName==elem], 
                              ',') ) )
  dis2proteinList[[elem]] <- x
}

signDf <- signDf[ , !sapply(signDf, is.list)]
write.table(signDf,
            file.path(figureDir(), 'Figure5/data/activating_tier_overlap.tsv' ),
            row.names = F, sep = '\t', quote = F
)

#######################################
#### Non-Cancer disease enrichment
#######################################
###
enrGfpList <- list()
dfList <- list()
for (bait in names(apmsPreys) ) {
  
  print(bait)
  
  aPreys <- gsub(';.*', '', apmsPreys[[bait]] )
  tPreys <- gsub(';.*', '', tidPreys[[bait]] )
  res <- disgenet2r::disease_enrichment(unique(c(aPreys,tPreys)),
                                        common_entities = 3 )
  enrGfpList[[paste0('combined_',bait)]] <- res
  df <- res@qresult
  if (nrow(df) > 0 ) {
    dfList[[paste0('combined_',bait) ]] <- df
  }
}

expDesign <- apmsDesign[apmsDesign$batch==1, ]

gfpDisEnrichmentDf <- do.call(rbind, dfList)
gfpDisEnrichmentDf$query <- gsub('\\..*', '', rownames(gfpDisEnrichmentDf) )
gfpDisEnrichmentDf$bait <- gsub('.*_', '', gfpDisEnrichmentDf$query  )
gfpDisEnrichmentDf$method <- gsub('_.*', '', gfpDisEnrichmentDf$query  )

idxs <- match(gfpDisEnrichmentDf$bait, expDesign$groups)
gfpDisEnrichmentDf$type <- expDesign$type[idxs]

write.table(gfpDisEnrichmentDf,
            file.path(figureDir(), 
                      'S5/data/disease_enrichment_interactome.tsv' ),
            row.names = F, sep = '\t', quote = F)
