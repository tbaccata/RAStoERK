#' @author Sebastian Didusch
#' @date 24.10.05
#' 
#' Script to provide:
#' Global variables,
#' Functions (pipeline, plotting, data analysis),
#' Colors

library(limma)
library(vsn)
library(DEqMS)
library(imputeLCMD)
library(igraph)

library(pheatmap)
# library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(reshape2)
library(GGally)


### COLORS ######################
##### BEGIN Colors
paraCols=c(KRAS='#f26924',
  NRAS='#f68b22',
  HRAS='#faa61e',
  ARAF='#174a92',
  BRAF='#4f88c7',
  RAF1='#a1bce2',
  MEK1='#009147',
  MEK2='#88be40',
  ERK1='#a62877',
  ERK2='#dfa0c2')
# names(paraCols) <- c('KRAS', 'NRAS', 'HRAS', 'ARAF', 'BRAF', 'RAF1', 'MEK1', 'MEK2', 'ERK1', 'ERK2')

groups <- c('#0077c0',
            '#00a2e4',
            '#7bc7ef',
            '#bcdef5')

groups <- c(
  '#f6891e',
  '#5c86af',
  '#44a843',
  '#c2649b'
)

names(groups) <- c('RAS', 'RAF', 'MEK', 'ERK')

annotation_colors <- list(
  type=c(WT='#a7a9ab', Inactive='#f0f1f1', Active='#fec111' ),
  paralog=paraCols,
  node = groups
)
##### END Colors
#################################

source("https://raw.githubusercontent.com/tbaccata/amica/refs/heads/master/amica/R/utils.R")
source("https://raw.githubusercontent.com/tbaccata/amica/refs/heads/master/amica/R/ProteomicsData.R")
# source('~/Documents/dev/shiny/amica/amica/R/utils.R')
# source('~/Documents/dev/shiny/amica/amica/R/ProteomicsData.R')

rawPrefix="^Intensity."
intensityPrefix="LFQ.intensity."
abundancePrefix="iBAQ"
razorUniqueCount="razorUniqueCount"
razorUniqueCountPrefix="razorUniqueCount."
spectraCount="spectraCount"
proteinId="Majority.protein.IDs"
geneName="Gene.names"
contaminantCol="Potential.contaminant"
pvalPrefix="P.Value_"
padjPrefix="adj.P.Val_"
logfcPrefix="logFC_"
avgExprPrefix="AveExpr_"
filterVal="+"

baitOrder <- c(
  'KrasWT', "KrasS17N",'KrasG12V','KrasG12D',
  "NrasWT", "NrasS17N", "NrasQ61K", "NrasQ61R",
  "HrasWT", "HrasS17N", "HrasG13R", "HrasQ61R",
  "ArafWT", "ArafD447A", "ArafS214P",
  "BrafWT", "BrafD594N", "BrafG469A", "BrafV600E",
  "Raf1WT", "Raf1D486A", "Raf1S257L", "Raf1S259F",
  "MEK1WT", "MEK1AA", "MEK1DD",
  "MEK2WT", "MEK2AA", "MEK2DD",
  "ERK1WT", "ERK1AA", "ERK1DD", 
  "ERK2WT", "ERK2AA", "ERK2DD"
)

################################################################################
imputeIntensities <- function(df_int,
                              method = "normal",
                              downshift = 1.8,
                              width = 0.3, seed = 12345) {
  set.seed(seed)
  
  if (method == "normal") {
    df_int[] <-
      lapply(df_int, function(x)
        replace(x, is.na(x), rnorm(
          mean = median(x, na.rm = T) - downshift * sd(x, na.rm = T),
          sd = sd(x, na.rm = T) * width,
          n = sum(is.na(x)) )
        )
      )
  }
  if (method == "min") {
    const <- min(df_int, na.rm = T)
    df_int[is.na(df_int)] <- const
  }
  if (method == "global") {
    medians <- apply(df_int, MARGIN = 2, FUN=median, na.rm=TRUE)
    median_global <- median(medians)
    sd_global <- sd(as.numeric(as.matrix(df_int)), na.rm=TRUE)
    
    mu_imputed <- median_global - downshift*sd_global
    sd_imputed <- width*sd_global
    
    df_int[] <-
      lapply(df_int, function(x)
        replace(x, is.na(x), rnorm(
          mean = mu_imputed,
          sd = sd_imputed,
          n = sum(is.na(x)) )
        )
      )
  }
  return(df_int)
}


renormalizeIntensities <- function(data, method = "None", proteinsOfInterest=NULL) {
  if (method == "Quantile") {
    data <- normalizeBetweenArrays(data, method = "quantile")
  } else if (method == "VSN") {
    data <- normalizeVSN(2 ^ data)
  } else if (method == "Median.centering") {
    data[] <- lapply(data, function(x) {
      gMedian = median(x, na.rm = TRUE)
      x - gMedian
    })
  } else if (method == 'Median.sweeping') {
    validRows <- rowSums(data, na.rm = T) > 0
    validDf <- data[validRows, ]
    
    medianDf <- apply(validDf, MARGIN = 2, FUN=median, na.rm=TRUE)
    gMedian <- median(medianDf)
    dMedians <- gMedian - medianDf
    
    data <- sweep(data, FUN="+", MARGIN = 2, STATS = dMedians)
  } else if (method == "renormalize2Proteins") {
    
    if (is.null(proteinsOfInterest)) {
      stop('Please provide valid protein ids of proteins you want to renormalize.')
    }
    
    rawPOIDf <- data[proteinsOfInterest,]
    validRows <- rowSums(rawPOIDf, na.rm = T) > 0
    
    validPOIDf <- rawPOIDf[validRows, ]
    
    medianPOIDf <- apply(validPOIDf, MARGIN = 2, FUN=median, na.rm=TRUE)
    gPOIMedian <- median(medianPOIDf)
    
    dMedians <- gPOIMedian - medianPOIDf
    
    data <- sweep(data, FUN="+", MARGIN = 2, STATS = dMedians)
    
  }
  return(as.data.frame(data))
}



readInSpectronaut <- function(protData, design, intensityThreshold = 1) {
  
  protData[[contaminantCol]] <- ""
  
  if (!("PG.ProteinAccessions" %in% names(protData))) {
    stop(paste0("Error: ", 
                "PG.ProteinAccessions", " is not a column name of the uploaded input data."))
  }
  if (!("PG.Genes" %in% names(protData))) {
    stop(paste0("Error: ", 
                "PG.Genes", " is not a column name of the uploaded input data."))
  }
  
  names(protData)[which(names(protData)=="PG.ProteinAccessions")] <- "Majority.protein.IDs"
  names(protData)[which(names(protData)=='PG.Genes')] <- "Gene.names"
  
  if ('PG.RunEvidenceCount' %in% names(protData)) {
    names(protData)[which(names(protData)=='PG.RunEvidenceCount')] <- "razorUniqueCount"
  }
  
  msmsIdx <- grep('PG.NrOfPrecursorsIdentified.*Experiment.wide.', names(protData))
  razIdx <- grep('PG.NrOfStrippedSequencesIdentified.*Experiment.wide.', names(protData))
  
  # keep birA
  protData$Majority.protein.IDs[protData$Gene.names=='birA'] <- 'P06709'
  protData$Majority.protein.IDs[protData$Gene.names=='GFP'] <- 'P42212'
  row.names(protData) <- protData$Majority.protein.IDs
  
  protData[grep('^contam', protData$Majority.protein.IDs), contaminantCol] <- '+'
  protData[grep('^contam', protData$Gene.names), contaminantCol] <- '+'
  
  if ( length( msmsIdx ) > 0 ) {
    names(protData)[msmsIdx] <- "spectraCount"
  }
  if ( length( razIdx ) > 0 ) {
    names(protData)[razIdx] <- "razorUniqueCount"
  }
  
  samples <- make.names(design$samples)
  
  ### razor/unique peptide count
  if (length(grep('NrOfStrippedSequencesUsedForQuantification',names(protData))) > 0 ) {
    idxs <- grep('NrOfStrippedSequencesUsedForQuantification',names(protData))
    
    if (length(idxs) != length(samples)) {
      stop(paste0("Error: Not all samples in experimental design could be matched in PG report."))
    }
    
    names(protData)[idxs] <- paste0('razorUniqueCount.', samples)
  }
  
  ### spectra counts
  if (length(grep('NrOfPrecursorsUsedForQuantification',names(protData))) > 0 ) {
    idxs <- grep('NrOfPrecursorsUsedForQuantification',names(protData))
    
    if (length(idxs) != length(samples)) {
      stop(paste0("Error: Not all samples in experimental design could be matched in PG report."))
    }
    
    names(protData)[idxs] <- paste0('spectraCount.', samples)
  }
  
  protData$Gene.names <- ifelse(
    protData$Gene.names == "",
    protData$Majority.protein.IDs,
    protData$Gene.names
  )
  
  sampleColumns <- grep('PG.Quantity', names(protData) )
  if (length(sampleColumns) != length(samples)) {
    stop(paste0("Error: Not all samples in experimental design could be matched in PG report."))
  }
  
  assayList <- list()
  intensities <- protData[, sampleColumns]
  names(intensities) <- samples
  intensities[intensities==0] <- NA
  intensities[!is.na(intensities) & intensities < intensityThreshold] <- NA
  intensities <- log2(intensities)
  assayList[['LFQIntensity']] <- intensities
  
  dropIdx <- sampleColumns
  
  if (length(assayList) < 1) {
    stop(paste0("Error: ", 
                "There are intensities in the uploaded Spectronaut file."))
  }
  
  se <- ProteomicsData(
    assays = assayList,
    rowData = protData[,-dropIdx],
    colData = design
  )
  
  se
}


# groupComparisons <-
#   function(imp_df_int, comparisons, expDesign, pep.count.table = NULL) {
#     df <- data.frame()
#     idx <- 1
#     
#     
#     for (i in 1:nrow(comparisons)) {
#       group1 <- comparisons[i, 1]
#       group2 <- comparisons[i, 2]
#       
#       group1names <- expDesign[expDesign$groups==group1, "samples"]
#       group2names <- expDesign[expDesign$groups==group2, "samples"]
#       
#       group1_idx <- grep(paste0("^",group1names,"$", collapse = "|"), colnames(imp_df_int))
#       group2_idx <- grep(paste0("^",group2names,"$", collapse = "|"), colnames(imp_df_int))
#       
#       rel_idxs <- c(group2_idx, group1_idx)
#       relevant_group_names <- c(rep(group2, length(group2_idx)), rep(group1, length(group1_idx)) )
#       comparison <- subset(imp_df_int, select = rel_idxs)
#       
#       limmaResults <- NULL #data.frame(row.names = rownames(imp_df_int) )
#       
#       if (length(group1_idx) < 2 | length(group2_idx) < 2) {
#         logFC <- subset(imp_df_int, select = group1_idx) - subset(imp_df_int, select = group2_idx)
#         
#         AveExpr <-
#           (subset(imp_df_int, select = group1_idx) + 
#              subset(imp_df_int, select = group2_idx)) / 2
#         
#         tmp <- data.frame(logFC, AveExpr)
#         
#         colnames(tmp) <- c("logFC", "AveExpr")
#         limmaResults <- rbind(limmaResults, tmp)
#       } else {
#         
#         class = factor(relevant_group_names, levels=c(group2, group1))
#         
#         
#         design = model.matrix(~0+class) # fitting without intercept
#         colnames(design) <- c(group2, group1)
#         
#         fit1 = lmFit(comparison, design = design)
#         contrastNames <- c(paste0(group1,"-",group2))
#         cont <- makeContrasts(contrasts=contrastNames, levels = colnames(design))
#         #cont <- makeContrasts(contrasts=eval(paste0(group1,"-",group2)), levels = colnames(design))
#         fit2 = contrasts.fit(fit1,contrasts = cont)
#         fit3 <- eBayes(fit2)
#         
#         limmaResults <- NULL
#         
#         if (!is.null(pep.count.table) ){
#           fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]
#           fit4 = spectraCounteBayes(fit3)
#           limmaResults = outputResult(fit4,coef_col = 1)
#           #limmaResults <- limmaResults[order(as.numeric( rownames(limmaResults) )), ]
#           limmaResults <- limmaResults[rownames(imp_df_int),]
#           limmaResults$P.Value <- limmaResults$sca.P.Value
#           limmaResults$adj.P.Val <- limmaResults$sca.adj.pval
#           limmaResults$sca.P.Value <- NULL
#           limmaResults$sca.adj.pval <- NULL
#           
#         } else {
#           limmaResults <- topTable(
#             fit3,
#             coef = 1,
#             number = Inf,
#             adjust = "BH",
#             sort.by = "none"
#           )
#         }
#       }
#       
#       colnames(limmaResults) <-
#         paste0(colnames(limmaResults), "_", group1, "__vs__", group2)
#       
#       if (idx < 2) {
#         df <- limmaResults
#       } else {
#         df <- cbind(df, limmaResults)
#       }
#       idx <- idx + 1
#     }
#     
#     rownames(df) <- rownames(imp_df_int)
#     return(df)
#   }


groupComparisonsBatchCorrected <-
  function(impDf, comparisons, design, pep.count.table = NULL) {
    df <- data.frame()
    idx <- 1
    
    
    for (i in 1:nrow(comparisons)) {
      group1 <- comparisons[i, 1]
      group2 <- comparisons[i, 2]
      
      #### limma
      group1names <- design[design$groups==group1, "samples"]
      group2names <- design[design$groups==group2, "samples"]
      
      group1_idx <- grep(paste0("^",group1names,"$", collapse = "|"), colnames(impDf))
      group2_idx <- grep(paste0("^",group2names,"$", collapse = "|"), colnames(impDf))
      
      rel_idxs <- c(group2_idx, group1_idx)
      relevant_group_names <- c(rep(group2, length(group2_idx)), rep(group1, length(group1_idx)) )
      comparison <- subset(impDf, select = rel_idxs)
      
      batches <- c()
      for (cname in colnames(comparison)) {
        batches <- c(batches, design$batch[design$samples == cname])
      }
      batches <- factor(batches)
      
      limmaResults <-  NULL #data.frame(row.names = rownames(impDf) )
      
      class = factor(relevant_group_names, levels=c(group2, group1))
      modelDesign <- model.matrix(~class)
      dupcor <- duplicateCorrelation(comparison,modelDesign,block=batches)
      
      fit1 <- lmFit(comparison, modelDesign, block= batches, correlation = dupcor$consensus)
      
      fit_ebayes <- eBayes(fit1)
      ####
      
      #### DEqMS
      fit_ebayes$count = pep.count.table[rownames(fit_ebayes$coefficients),"count"]
      fit4 = spectraCounteBayes(fit_ebayes)
      limmaResults = outputResult(fit4,coef_col = 2)
      #limmaResults <- limmaResults[order(as.numeric( rownames(limmaResults) )), ]
      limmaResults <- limmaResults[rownames(impDf),]
      limmaResults$P.Value <- limmaResults$sca.P.Value
      limmaResults$adj.P.Val <- limmaResults$sca.adj.pval
      limmaResults$sca.P.Value <- NULL
      limmaResults$sca.adj.pval <- NULL
      limmaResults$sca.t <- NULL
      limmaResults$gene <- NULL
      ####
      
      colnames(limmaResults) <-
        paste0(colnames(limmaResults), "_", group1, "__vs__", group2)
      
      if (idx < 2) {
        df <- limmaResults
      } else {
        df <- cbind(df, limmaResults)
      }
      idx <- idx + 1
    }
    rownames(df) <- rownames(impDf)
    return(df)
  }


preProcess <- function(protData, design, intensityThreshold, baitGroups=NULL, minMSMS=3, minRazor=2, minValueGroup=3) {
  se <- readInSpectronaut(protData, design, intensityThreshold=intensityThreshold )
  
  impDf <- assay(se, 'LFQIntensity')
  
  se <- setAssay(x = se,
                 assay = impDf,
                 assayName = 'LFQIntensity')
  
  rnames <- filterOnMinValuesRnames(
    y = se,
    minMSMS = minMSMS,
    minRazor = minRazor
  )
  
  ##
  tmp <- filterOnValidValues(
    data = impDf[rnames, ],
    mappings = colData(se),
    groupsToConsider = baitGroups,
    minValue = minValueGroup,
    method = "in_one_group")
  
  
  rowsDf <- rowData(se)
  rowsDf$quantified <- ""
  rowsDf[tmp, 'quantified'] <- "+"
  
  se <-
    setRowData(se, rowsDf)
  
  se
}

#' Compute pca plot
#' 
#' @param plotData data.frame of intensities in wide format.
#' @param myGroupColors vector of colors.
#' @param boxplot_base int describing the base plot font size.
#' @param boxplot_legend int describing the legend size.
#' @param assayNames char 
#' @param groupFactors vector of groups to be ordered.
#' @param groupInputs vector of groups to be plotted.
#' @return list of df_out containing coordinates and percentages char vector of % variance.
#' @examples
computePCA <-
  function(plotData,
           expDesign,
           groupFactors = NULL,
           groupInputs = NULL) {
    # filter only by selected groups
    if (is.null(groupInputs))
      groupInputs <- unique(expDesign$groups)
    
    if (!is.null(groupFactors) &&
        all(groupFactors %in% unique(expDesign$groups))) {
      expDesign$groups <-
        factor(expDesign$groups, levels = c(groupFactors,
                                            setdiff(groupInputs, groupFactors)))
    }
    
    plotData <-
      plotData[, expDesign$samples[expDesign$groups %in% groupInputs]]
    plotData <- plotData[complete.cases(plotData),]
    
    # compute PCA
    pca <- prcomp(as.data.frame(t(plotData)))
    df_out <- as.data.frame(pca$x)
    
    df_out$group <-
      expDesign$groups[expDesign$samples %in% row.names(df_out)]
    df_out$sample <- rownames(df_out)
    df_out$key <- df_out$sample
    df_out$show_id <- FALSE
    
    tmp <- summary(pca)
    percentage <-
      round(100 * tmp$importance['Proportion of Variance',], 2)
    percentage <-
      paste(colnames(df_out), "(", paste(as.character(percentage), "%", ")", sep =
                                           ""))
    return(list(df_out = df_out, percentage = percentage))
  }


################################################################################

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

calc_pairwise_overlaps <- function(sets) {
  # Ensure that all sets are unique character vectors
  sets_are_vectors <- vapply(sets, is.vector, logical(1))
  if (any(!sets_are_vectors)) {
    stop("Sets must be vectors")
  }
  sets_are_atomic <- vapply(sets, is.atomic, logical(1))
  if (any(!sets_are_atomic)) {
    stop("Sets must be atomic vectors, i.e. not lists")
  }
  sets <- lapply(sets, as.character)
  is_unique <- function(x) length(unique(x)) == length(x)
  sets_are_unique <- vapply(sets, is_unique, logical(1))
  if (any(!sets_are_unique)) {
    stop("Sets must be unique, i.e. no duplicated elements")
  }
  
  n_sets <- length(sets)
  set_names <- names(sets)
  n_overlaps <- choose(n = n_sets, k = 2)
  
  vec_name1 <- character(length = n_overlaps)
  vec_name2 <- character(length = n_overlaps)
  vec_num_sample1 <- integer(length = n_overlaps)
  vec_num_sample2 <- integer(length = n_overlaps)
  vec_num_shared <- integer(length = n_overlaps)
  vec_overlap <- numeric(length = n_overlaps)
  vec_jaccard <- numeric(length = n_overlaps)
  overlaps_index <- 1
  
  for (i in seq_len(n_sets - 1)) {
    name1 <- set_names[i]
    set1 <- sets[[i]]
    for (j in seq(i + 1, n_sets)) {
      name2 <- set_names[j]
      set2 <- sets[[j]]
      
      set_intersect <- set1[match(set2, set1, 0L)]
      set_union <- .Internal(unique(c(set1, set2), incomparables = FALSE,
                                    fromLast = FALSE, nmax = NA))
      num_shared <- length(set_intersect)
      overlap <- num_shared / min(length(set1), length(set2))
      jaccard <- num_shared / length(set_union)
      
      vec_name1[overlaps_index] <- name1
      vec_name2[overlaps_index] <- name2
      vec_num_sample1[overlaps_index] <- length(set1)
      vec_num_sample2[overlaps_index] <- length(set2)
      vec_num_shared[overlaps_index] <- num_shared
      vec_overlap[overlaps_index] <- overlap
      vec_jaccard[overlaps_index] <- jaccard
      
      overlaps_index <- overlaps_index + 1
    }
  }
  
  result <- data.frame(sample1 = vec_name1,
                       sample2 = vec_name2,
                       num_sample1 = vec_num_sample1,
                       num_sample2 = vec_num_sample2,
                       num_shared = vec_num_shared,
                       overlap = round(vec_overlap, 3),
                       jaccard = round(vec_jaccard, 3)
  )
  return(result)
}
################################################################################
