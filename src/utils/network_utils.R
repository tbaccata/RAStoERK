library(RColorBrewer)
library(igraph)

ppi <- read.table(
  file.path(dataDir(), 'networks/intact.edgelist'),
  sep = '\t',
  header = F
)
ppi <- graph_from_data_frame(ppi, directed = F)

createGraphAttributes <- function(G) {
  if ('perturbed' %in% names(edge_attr(G)))  E(G)$color <- ifelse(E(G)$perturbed == 'enriched', '#ed2939', '#08519c')
  if('logFC' %in% names(edge_attr(G))) E(G)$width <- ifelse(abs(E(G)$logFC) > 3, 2, 0.2)
  
  if('significant' %in% names(edge_attr(G))) {
    E(G)$color <- ifelse(E(G)$significant == 'core', '#33A02C', 
                         ifelse(E(G)$significant == 'shared', '#B2DF8A',
                         ifelse(E(G)$significant == 'AP-MS', 'pink',
                         'skyblue'
                         )))
  } 
  
  V(G)$shape <- ifelse(V(G)$name %in% baitOrder,
                       "square", "circle")
  
  V(G)$vertex.size <- ifelse(V(G)$name %in% baitOrder,
                             6, 3)
  
  x <- rep(0, length(degree(G)))
  x[degree(G) < 2] <- 1.5
  x[degree(G) == 2] <- 1.75
  x[degree(G) == 3] <- 2
  x[degree(G) == 4] <- 2.25
  x[degree(G) == 5] <- 2.5
  x[degree(G) == 6] <- 2.75
  x[degree(G) == 7] <- 3
  x[degree(G) == 8] <- 3.25
  x[degree(G) == 9] <- 3.5
  x[degree(G) == 10] <- 3.75
  x[degree(G) == 11] <- 4
  x[degree(G) > 10] <- 5
  x[grep('kras|nras|hras|araf|braf|raf1|mek1|mek2|erk1|erk2', V(G)$name, ignore.case = T) ] <- 7
  
  V(G)$vertex.size <- x
  V(G)$vertex.color <- 'black'
  V(G)$vertex.color[grep('Kras', V(G)$name) ] <- '#a50f15'
  V(G)$vertex.color[grep('Nras', V(G)$name) ] <- '#de2d26'
  V(G)$vertex.color[grep('Hras', V(G)$name) ] <- '#fb6a4a'
  V(G)$vertex.color[grep('Araf', V(G)$name) ] <- '#08519c'
  V(G)$vertex.color[grep('Braf', V(G)$name) ] <- '#3182bd'
  V(G)$vertex.color[grep('Raf1', V(G)$name) ] <- '#6baed6'
  V(G)$vertex.color[grep('MEK1', V(G)$name) ] <- '#31a354'
  V(G)$vertex.color[grep('MEK2', V(G)$name) ] <- '#74c476'
  V(G)$vertex.color[grep('ERK1', V(G)$name) ] <- '#756bb1'
  V(G)$vertex.color[grep('ERK2', V(G)$name) ] <- '#9e9ac8'
  G
}

####
####
corumNetwork <- read.table(file.path(dataDir(), 'networks/corum.edgelist'),
                           header = F, sep = '\t')
signorNetwork <- read.table(file.path(dataDir(), 'networks/signor.edgelist'),
                            header = F, sep = '\t')
G_corum <- graph_from_data_frame(corumNetwork, directed = F)
G_signor <- graph_from_data_frame(signorNetwork, directed = F)

annotateCorumSignorNetwork <- function(G) {
  if ('perturbed' %in% names(edge_attr(G)))  E(G)$color <- ifelse(E(G)$perturbed == 'enriched', '#ed2939', '#08519c')
  if('logFC' %in% names(edge_attr(G))) E(G)$width <- ifelse(abs(E(G)$logFC) > 3, 2, 0.2)
  
  E(G)$width <- 1.5
  
  if('type' %in% names(edge_attr(G))) {
    E(G)$color <- ifelse(E(G)$type == 'core', 'grey80', 
                         ifelse(E(G)$type == 'shared', 'grey80',
                                ifelse(E(G)$type == 'AP-MS', 'skyblue',
                                       ifelse(E(G)$type == 'TbID', 'pink',
                                              ifelse(E(G)$type == 'CORUM', 'black',
                                                     ifelse(E(G)$type == 'SIGNOR', '#33a02c', 'grey60')
                                              )))))
  } 
  
  V(G)$shape <- ifelse(V(G)$name %in% baitOrder,
                       "square", "circle")
  
  V(G)$vertex.color <- 'black'
  V(G)$vertex.color[V(G)$name %in% baitOrder ] <- 'grey80'
  V(G)$label.family <- 'Helvetica'
  G
}


getCorumPpiNetwork <- function(interactomeSubset, include = 'CORUM') {
  ppi <- interactomeSubset[, c('bait', 'prey', 'significant')]
  names(ppi) <- c('from', 'to', 'type')
  
  corumEdges <- igraph::as_data_frame(
    simplify(induced_subgraph(G_corum, which(V(G_corum)$name %in% unique(ppi$to) ) ) )
  )
  if (nrow(corumEdges) > 0)  corumEdges$type <- 'CORUM'
  #
  signorEdges <- igraph::as_data_frame(
    simplify(induced_subgraph(G_signor, which(V(G_signor)$name %in% unique(ppi$to) ) ) )
  )
  if (nrow(signorEdges) > 0)  signorEdges$type <- 'SIGNOR'
  
  out <- ppi
  if (include == 'CORUM') {
    out <- rbind(out, corumEdges)
  } else if (include == 'SIGNOR') {
    out <- rbind(out, signorEdges)
  } else {
    out <- rbind(out, corumEdges)
    out <- rbind(out, signorEdges)
  }
  out
}


createPieChartNetork <- function(baseNetwork, physicalNetwork, subInteractome, preysOfTerm) {
  #' @param baseNetworkical igraph network with phys. and funct. edges
  #' @param physicalNetwork igraph network with phys. edges only
  #' @param subInteractome interactome
  #' @param preysOfTerm vector of node names to be searched in baseNetwork
  #' 
  #' return annotated network for CytoScape visualization
  allTiers <- c('RAS', 'RAF', 'MEK', 'ERK')
  idxs <- which(V(baseNetwork)$name %in% preysOfTerm)
  
  G_term <- induced_subgraph(baseNetwork, vids = idxs)
  #
  # V(G_term)$organelle <- ifelse(V(G_term)$name %in% tmp$prey, tmp$LeonettiLocalization[tmp$prey %in%
  #                                                                                        preysOfTerm], '')
  #
  natOrder <- names(annotation_colors$paralog)
  
  # Initialize attributes
  V(G_term)$signProts <- NA
  V(G_term)$nTiers    <- NA
  V(G_term)$tiers    <- NA
  V(G_term)$specific  <- NA
  V(G_term)$RAS  <- NA
  V(G_term)$RAF  <- NA
  V(G_term)$MEK  <- NA
  V(G_term)$ERK  <- NA
  #
  V(G_term)$KRAS  <- NA
  V(G_term)$NRAS  <- NA
  V(G_term)$HRAS  <- NA
  #
  V(G_term)$ARAF  <- NA
  V(G_term)$BRAF  <- NA
  V(G_term)$RAF1  <- NA
  #
  V(G_term)$MEK1  <- NA
  V(G_term)$MEK2  <- NA
  #
  V(G_term)$ERK1  <- NA
  V(G_term)$ERK2  <- NA
  
  V(G_term)$method  <- NA
  
  cras <- craf <- cmek <- cerk <- c()
  for (prey in unique(subInteractome$prey)) {
    prots <- subInteractome$bait[subInteractome$prey == prey]
    
    tiers <- unique(apmsDesign$node[apmsDesign$groups %in% prots])
    prots <- unique(apmsDesign$paralog[apmsDesign$groups %in% prots])
    
    V(G_term)$RAS[V(G_term)$name == prey] <- ifelse('RAS' %in% tiers, 1, 0)
    V(G_term)$RAF[V(G_term)$name == prey] <- ifelse('RAF' %in% tiers, 1, 0)
    V(G_term)$MEK[V(G_term)$name == prey] <- ifelse('MEK' %in% tiers, 1, 0)
    V(G_term)$ERK[V(G_term)$name == prey] <- ifelse('ERK' %in% tiers, 1, 0)
    
    #
    V(G_term)$KRAS[V(G_term)$name == prey] <- ifelse('KRAS' %in% prots, 1, 0)
    V(G_term)$NRAS[V(G_term)$name == prey] <- ifelse('NRAS' %in% prots, 1, 0)
    V(G_term)$HRAS[V(G_term)$name == prey] <- ifelse('HRAS' %in% prots, 1, 0)
    
    V(G_term)$ARAF[V(G_term)$name == prey] <- ifelse('ARAF' %in% prots, 1, 0)
    V(G_term)$BRAF[V(G_term)$name == prey] <- ifelse('BRAF' %in% prots, 1, 0)
    V(G_term)$RAF1[V(G_term)$name == prey] <- ifelse('RAF1' %in% prots, 1, 0)
    
    V(G_term)$MEK1[V(G_term)$name == prey] <- ifelse('MEK1' %in% prots, 1, 0)
    V(G_term)$MEK2[V(G_term)$name == prey] <- ifelse('MEK2' %in% prots, 1, 0)
    
    V(G_term)$ERK1[V(G_term)$name == prey] <- ifelse('ERK1' %in% prots, 1, 0)
    V(G_term)$ERK2[V(G_term)$name == prey] <- ifelse('ERK2' %in% prots, 1, 0)
    
    signProts <- ''
    if (length(prots) > 1) {
      signProts <- paste(natOrder[natOrder %in% prots], collapse = '_')
      V(G_term)$specific[V(G_term)$name == prey] <- 'shared'
    } else {
      signProts <- prots
      V(G_term)$specific[V(G_term)$name == prey] <- prots
    }
    V(G_term)$tiers[V(G_term)$name == prey] <- paste(allTiers[allTiers %in% tiers], collapse = '_')
    
    nTiers <- length(prots)
    # Assign attributes to the matching node
    V(G_term)$signProts[V(G_term)$name == prey] <- signProts
    V(G_term)$nTiers[V(G_term)$name == prey] <- nTiers
  }
  
  cntDf <- data.frame(V(G_term)$KRAS, V(G_term)$NRAS, V(G_term)$HRAS,
             V(G_term)$ARAF, V(G_term)$BRAF, V(G_term)$RAF1,
             V(G_term)$MEK1, V(G_term)$MEK2,
             V(G_term)$ERK1, V(G_term)$ERK2
             )
  # for cytoscape.js - % of "signal" between 0-100
    V(G_term)$pie1 <- V(G_term)$KRAS/apply(cntDf, 1, sum) * 100
    V(G_term)$pie2 <- V(G_term)$NRAS/apply(cntDf, 1, sum) * 100
    V(G_term)$pie3 <- V(G_term)$HRAS/apply(cntDf, 1, sum) * 100
    
    V(G_term)$pie4 <- V(G_term)$ARAF/apply(cntDf, 1, sum) * 100
    V(G_term)$pie5 <- V(G_term)$BRAF/apply(cntDf, 1, sum) * 100
    V(G_term)$pie6 <- V(G_term)$RAF1/apply(cntDf, 1, sum) * 100
    
    V(G_term)$pie7 <- V(G_term)$MEK1/apply(cntDf, 1, sum) * 100
    V(G_term)$pie8 <- V(G_term)$MEK2/apply(cntDf, 1, sum) * 100
    
    V(G_term)$pie9 <- V(G_term)$ERK1/apply(cntDf, 1, sum) * 100
    V(G_term)$pie10 <- V(G_term)$ERK2/apply(cntDf, 1, sum) * 100
  
  E(G_term)$physical <- "no"  # default
  for (e in E(G_term)) {
    v1 <- ends(G_term, e)[1]
    v2 <- ends(G_term, e)[2]
    if (are_adjacent(physicalNetwork, v1, v2)) {
      E(G_term)[e]$physical <- "yes"
    }
  }
  G_term
}


if (!require(httr)) install.packages("httr")
if (!require(jsonlite)) install.packages("jsonlite")

stringNetworkFromAPI <- function(proteins, network_type = 'physical') {
  #'@param proteins list of proteins to query (HGNC)
  #'@param network_type either "physical" or "functional"
  #'
  #'return igraph object of STRING network
  
  if (network_type != "physical" && network_type !=  "functional") {
    stop('network_type must be either "physical" or "functional"' )
  }
  species <- 9606  # Human
  min_score <- 400
  string_version <- "12.0"
  base_url <- paste0("https://www.string-db.org/api/")
  network_url <- paste0(base_url, "json/network")
  
  version <- 'v12.0'
  ids <- paste(proteins, collapse = "%0d")
  
  network_resp <- GET(network_url,
                      query = list(
                        identifiers = ids,
                        species = species,
                        required_score = min_score,
                        network_type = network_type,
                        version=version
                      ))
  
  network <- fromJSON(content(network_resp, "text", encoding = "UTF-8"))
  
  edges <- data.frame(
    protein1 = network$preferredName_A,
    protein2 = network$preferredName_B,
    score = network$score,
    stringsAsFactors = FALSE
  )
  G <- graph_from_data_frame(edges, directed = F)
  
  missing <- setdiff(proteins, V(G)$name)
  
  if ( length(missing) > 0 ) {
    G <- add_vertices(G, nv = length(missing), attr = list(name = missing ))
  }
  G
}
