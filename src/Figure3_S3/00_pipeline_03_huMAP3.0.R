#' @author Sebastian Didusch
#' @date 25.08.24
#' 
#' integrate huMAP3.0
#' 1) upload custom gmt to gprofiler
#' 2) retrieve significant humap3.0 complexes using gost
#' 3) all sign. complexes saved in file.path(figureDir(), "Figure3/data/all_enriched_humap3.0_networks.pdf")
#' 4) selected complexes exported to graphml for further visualization using networkx 
#' (see src/notebooks) 


library(gprofiler2)

## define paths
source(here::here("src/utils/utils_path.R"))
source(file.path(srcDir(), "utils/enrichment_utils.R"))

######### Only need to do this once
######### Upload a custom gmt file to gprofiler
######### id can be saved and reused

# custom_id <- upload_GMT_file(file.path(dataDir(), 'integration/humap3.0/humap3.0.gmt'))

# Your custom annotations ID is gp__lCNf_VqBL_OwM.
# You can use this ID as an 'organism' name in all the related enrichment tests against this custom source.
# Just use: gost(my_genes, organism = 'gp__lCNf_VqBL_OwM')

uniqueDesign <- apmsDesign[apmsDesign$type != 'Control' & apmsDesign$batch == 1, ]
#
outDf <- list()
nodes <- c('RAS', 'RAF', 'MEK', 'ERK')
types <- unique(apmsDesign$type[apmsDesign$type != 'Control'] )

for (node in nodes) {
  for (type in types) {
    # (only one KRAS even though two KRAS mutants - converge on protein level)
    paras <- unique(apmsDesign$paralog[apmsDesign$type == type & apmsDesign$node == node])

    para2preys <- list()
    for (para in paras) {
      baits <- uniqueDesign$groups[uniqueDesign$paralog == para & uniqueDesign$type == type]
      preys <- unique(c(unlist(apmsPreys[baits] ),
                        unlist(tidPreys[baits]) ) )

      outDf[[paste0(type, '_', para, '_combined') ]] <- gsub(';.*', '', preys) # overlappingPreys
      #outDf[[paste0(type, '_', para, '_apms') ]] <- gsub(';.*', '',  unique(c(unlist(apmsPreys[baits]))) )
      #outDf[[paste0(type, '_', para, '_tid') ]] <- gsub(';.*', '', unique(c(unlist(tidPreys[baits])) ) )
    }
  }
}

humapGprofiler <- gost(outDf, organism = "gp__lCNf_VqBL_OwM", evcodes = T,
                      custom_bg = background, significant = T)
humapComplexes <- humapGprofiler$result
humapComplexes$parents <- NULL

write.table(humapComplexes,
            file.path(figureDir(), 'Figure3/data/humap3.0_enrichment.tsv'),
            row.names = F, sep = '\t', quote = F)

uniqueDesign <- apmsDesign[apmsDesign$type != 'Control' & apmsDesign$batch == 1, ]

### 
humap2preys <- list()
for (cid in unique(humapComplexes$term_id)) {
  preys <- humapComplexes$intersection[humapComplexes$term_id == cid] 
  preys <- unique(unlist(strsplit(preys, ',' ) ) )
  humap2preys[[cid]] <- preys
}

getRedundantTerms <- function(humap2preys) {
  #' function that calculates pairwise overlap coeff.
  #' and builds a directed graph of terms that have 100% overlap.
  #' source is always the bigger term, 
  #' the higher confident term or at random if terms have the same size.
  #' it returns all child nodes, 
  
  simDf <- calc_pairwise_overlaps(humap2preys)
  identicals <- simDf[simDf$overlap == 1, ]
  
  df <- data.frame(from='A', to='B')
  for (idx in 1:nrow(identicals )) {
    example <- identicals[idx, ]
    conf1 <- as.integer(gsub('.*_', '', example$sample1) )
    conf2 <- as.integer(gsub('.*_', '', example$sample2) )
    
    if (example$num_sample1 == example$num_sample2) {
      if (conf1 < conf2) {
        df <- rbind(
          df,
          data.frame(from=example$sample1, to=example$sample2 )
        )
      } else { 
        # if conf 2 < conf 1 take sample2 as parent,
        # otherwise take sample 2 as parent (if equal conf and size)
        df <- rbind(
          df,
          data.frame(from=example$sample2, to=example$sample1 )
        )
      } 
    } else if (example$num_sample1 > example$num_sample2) {
      df <- rbind(
        df,
        data.frame(from=example$sample1, to=example$sample2 )
      )
    } else {
      df <- rbind(
        df,
        data.frame(from=example$sample2, to=example$sample1 )
      )
    }
  }
  df <- df[df$from != 'A', ]
  redudancyGraph <- graph_from_data_frame(df, directed = T)
  #child_nodes <- which( degree(redudancyGraph, v = V(redudancyGraph), mode = "in" ) > 0)
  #child_nodes
  redudancyGraph
}

redundancyGraph <- getRedundantTerms(humap2preys)
child_nodes <- which( degree(redundancyGraph, v = V(redundancyGraph), mode = "in" ) > 0)

V(redundancyGraph)$color <- 'darkred'
V(redundancyGraph)$color[child_nodes] <- 'grey80'

# Visualize redudant huMAP3.0 complexes
# plot(redundancyGraph, vertex.size=3, edge.arrow.size=0.01, vertex.label.cex=0.003)
ccs <- igraph::components(redundancyGraph)

cc2mat <- list()
for (idx in unique(ccs$membership) ) {
  sign_terms <- names(ccs$membership[ccs$membership==idx] )
  tmpOra <- humapComplexes[humapComplexes$term_name %in% sign_terms, ]
  sign_baits <- tmpOra$query
  proteins <- (unlist(lapply(strsplit(sign_baits, '_'), function(x) x[2] ) ) )
  tiers <- unique(apmsDesign$node[apmsDesign$paralog %in% proteins] )
  
  trm2preys <- list()
  for (trm in sign_terms) {
    x <- unique(unlist(strsplit(tmpOra$intersection[tmpOra$term_name == trm], ',' )))
    trm2preys[[trm]] <- x
  }
  
  preys <- sort(unique(unlist(trm2preys)))
  mat <- t(+sapply(trm2preys, "%in%", x = preys))  ## matrix output
  colnames(mat) <- preys
  
  mat <- t(mat)
  mat <- as.data.frame(mat)
  
  signIn <- c()
  for (rname in rownames(mat) ) {
    baits <- interactome$bait[interactome$prey == rname]
    sign <- paste(unique(apmsDesign$paralog[apmsDesign$groups %in% baits]), 
                  collapse = ',')
    signIn <- c(signIn, sign)
  }
  mat$baits <- signIn
  #
  if (length(tiers) > 1)  cc2mat[[idx]] <- mat
}

saveRDS(cc2mat, 
        file.path(figureDir(), 
                          'Figure3/data/humap3.0_redudant_complexes.rds') )

#humapComplexes[humapComplexes$term_name %in% names(ccs$membership[ccs$membership==which.max(ccs$csize) ] ), ]
nonRedHumap <- humapComplexes[!humapComplexes$term_id %in% names(child_nodes), ]


### 
termList <- list()
info <- list()
for (termName in unique(nonRedHumap$term_name)) {
  tmpOra <- nonRedHumap[nonRedHumap$term_name == termName, ] 
  preys <- nonRedHumap$intersection[nonRedHumap$term_name == termName]
  
  #if (length(grep('^POL' , preys ) ) > 0 ) next
  
  preys <- preys %>% strsplit(',') %>% unlist() %>% unique()
  complexInteractome <- interactome[interactome$prey %in% preys, ]
  
  if ( (max(tmpOra$recall) < 0.75 && max(tmpOra$intersection_size) < 5 ) ||
       max(tmpOra$term_size) < 3 || max(tmpOra$intersection_size) < 2 ) next
  
  conf <- as.integer(gsub('.*_', '', tmpOra[1,]$term_id) )
  charOut <- ''
  if (conf == 1) {
    charOut <- 'Extremely High'
  } else if (conf == 2) {
    charOut <- 'Very High'
  } else if (conf == 3) {
    charOut <- 'High'
  } else if (conf == 4) {
    charOut <- 'Moderately High'
  } else if (conf == 5) {
    charOut <- 'Medium High'
  } else {
    charOut <- 'Medium'
  }
  
  
  signBaits <-strsplit(tmpOra$query, '_')
  signBaits <- unique(sapply(signBaits, function(x) x[2] ) )
  info[[termName]] <- paste('term size: ', max(tmpOra$term_size), 
                            ';max coverage: ',
                            round(max(tmpOra$recall), 2),
                            'significant in: ', 
                            paste(signBaits, collapse = ',') )
  
  outList <- list()
  for (type in c('Inactive', 'WT', 'Active')) {
    baits <- unique(apmsDesign$groups[apmsDesign$type == type] )
    complexTypeInteractome <- complexInteractome[complexInteractome$bait %in% baits, ]
    
    if (nrow(complexTypeInteractome) < 1) {
      outList[[type]] <- make_empty_graph(n = 0, directed = F)
      next
    }
    ppi <- getCorumPpiNetwork(complexTypeInteractome, include = 'all')
    G_complex <- graph_from_data_frame(ppi, directed = F)
    G_complex <- annotateCorumSignorNetwork(G_complex)
    
    outList[[type]] <- G_complex
  }
  termList[[termName]] <- outList
}

saveRDS(termList,
        file.path(figureDir(), 'Figure3/data/bait_humap3.0ora_termslist.rds'))

simGraph <- calc_pairwise_overlaps(humap2preys[names(termList)] )
g <- graph_from_data_frame(simGraph[simGraph$overlap > 0.3, c('sample1', 'sample2')],
                           directed = F
)
clst <- igraph::cluster_louvain(g)

newOrder <- c()
for (idx in unique(clst$membership) ) {
  terms <- clst$names[clst$membership == idx]
  newOrder <- c(newOrder, names(rev(sort(sapply(humap2preys[terms], length) ) ) ) )
  
  for (trm in terms) {
    info[[trm ]] <- paste0(info[[trm ]], ' cluster:', idx)
  }
}
newOrder <- c(newOrder, setdiff(names(termList), clst$names  ) )

###
termList <- termList[newOrder]
info <- info[newOrder]

dev.off()
pdf(file.path(figureDir(), "Figure3/data/all_enriched_humap3.0_networks.pdf"), width = 8, height = 8)  # open PDF device
for (k in seq_along(termList)) {
  
  nets <- termList[[k]]
  
  par(
    mfrow = c(1, 3),        # always 3 per row
    mar = c(1,1,1,1), 
    oma = c(4,0,4,0)        # space for title
  )
  
  for (i in 1:3) {
    g <- nets[[i]]
    if (is.null(g) || igraph::vcount(g) == 0) {
      plot.new()
    } else {
      E(G_complex)$width <- 2
      plot(
        g,
        layout = layout_in_circle,
        vertex.label.cex = 1, #V(G_kinases)$vertex.label.cex,
        vertex.label.dist = -2.5,
        edge.arrow.size = 0.001,
        #edge.width=E(G_kinases)$width,
        edge.color = E(g)$color,
        vertex.size = 15,
        vertex.label.color = 'black',
        vertex.color = V(g)$vertex.color
        #main = termName
      )
    }
  }
  
  term_name <- names(termList)[k]
  
  conf <- as.integer(gsub('.*_', '', term_name) )
  term_name <- gsub('.{2}$', '', term_name)
  
  charOut <- ''
  
  if (conf == 1) {
    charOut <- 'Extremely High'
  } else if (conf == 2) {
    charOut <- 'Very High'
  } else if (conf == 3) {
    charOut <- 'High'
  } else if (conf == 4) {
    charOut <- 'Moderately High'
  } else if (conf == 5) {
    charOut <- 'Medium High'
  } else {
    charOut <- 'Medium'
  }
  
  mtext(paste0(term_name, ', ', 
               charOut, ', | '), outer = TRUE, cex = 0.8, line = 2)
  mtext(info[k], outer = TRUE, side = 1, line = 2, cex = 0.8)
}
dev.off()  # close PDF device

################################
###
### bring complexes to Cytoscape
###
################################
p <- file.path(dataDir(), 'integration/humap3.0/humap3.0_ora_results_selection.dat')

dat <- read.delim(p, header = T, sep = '\t')

complexes <- c('PI3K', 'PKC', 'ERK signaling',
               'SAGA complex', 'RAF signaling', 'phosphatases',
               'Cell−extracellular matrix interactions')

outList <- list()
for (cplx in complexes) {
  ids <- dat$id[dat$annotation==cplx]
  
  termIds <- c()
  for (cid in ids) {
    tid <-  grep(cid, names(termList), value = T)
    termIds <- c(termIds, tid)
  }
  
  catList <- list()
  for (hid in termIds) {
    tmp <- termList[[hid]]
    dfInactive <- igraph::as_data_frame(tmp[['Inactive']])
    dfWt <- igraph::as_data_frame(tmp[['WT']])
    dfActive <- igraph::as_data_frame(tmp[['Active']])
    
    dfInactive$network <- 'Inactive'
    dfWt$network <- 'WT'
    dfActive$network <- 'Active'
    
    df <- rbind(
      rbind(dfWt, dfInactive), 
      dfActive)
    df$category <- cplx
    catList[[hid]] <- df
  }
  df <- do.call(rbind, catList)
  
  df$dupl <-  paste0(df$from,'_',df$to,'_',df$network, '_', df$type)
  df <- df[!duplicated(df$dupl),]
  
  outList[[cplx]] <- df
}
df <- do.call(rbind, outList)

df$strict_type <- df$type
df$strict_type[df$strict_type %in% c('core', 'shared')] <- 'both'

df$dupl <- NULL

df$from <- toupper(df$from)
df$to <- toupper(df$to)

write.table(df,
            file.path(dataDir(), 'processed/cytoscape/selected_humap3.0_complexes.tsv'),
            row.names = F, sep = '\t', quote = F
)
