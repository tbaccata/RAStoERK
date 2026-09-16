#' @author Sebastian Didusch
#' @date 25.10.03
#' 
#' write interactome to file after
#' Figure 5 has been processed
#' (in Fig. 5 disease genes are queried from DISGENET)
#' 


library(here)

source(here::here("src/utils/utils_path.R"))

p <- file.path(dataDir(), 'processed/proteome_gene2disease_curated_score_gt_0.3.tsv' )
if (!file.exists(p)) {
  stop(paste0('DISGENET has not been queried yet!\n',
              'Please execute src/pipeline/02_disgenet.R first.'
  ))
}
dis2proteome <- read.delim(p, header = T, sep = '\t')
interactome <- read.table(file.path(interimDataDir(), 'erkpathmap_interactome.tsv'),
                          header = T, sep = '\t')

p <- file.path(dataDir(), 'integration/subcell_localization/organelle_ip_umap.tsv')
if (!file.exists(p)) {
  stop(paste0('organelle IP has not been processed yet!\n',
              'Please execute src/Figure4_S4/00_pipeline_04_subcell_ora.R first.'
  ))
}
ips <- read.table(file.path(dataDir(), 'integration/subcell_localization/organelle_ip_umap.tsv'),
                  header = T, sep = '\t')

# only highlight canonical swissprot id for core pathway
ips$Gene_name_canonical[ips$Majority.protein.IDs=='A0A024RAV5;A0A3G1LBH1;Q14015;Q14014'] <- 'KRAS_trembl'
ips$Gene_name_canonical[ips$Majority.protein.IDs%in%c('H7C560;A0A2R8Y8E0', 'D7RF68;A0A2R8Y492')] <- 'BRAF_trembl'
ips$Gene_name_canonical[ips$Majority.protein.IDs%in%c('A0A384P5S9;A0A024R889;Q9UG16', 'A0A0D9SGF6;B4DGT1')] <- 'SPTAN1_trembl'

tmpIps <- ips[, c('Gene_name_canonical', 'graph_localization_annotation') ]
tmpIps <- tmpIps[!duplicated(tmpIps$Gene_name_canonical), ]

hcDisGenes <- unique(dis2proteome$gene_symbol[dis2proteome$score > 0.66] )

interactome$HCDP <- ifelse(
  gsub(';.*', '', interactome$prey)%in%hcDisGenes,
  'yes', 'no'
) 

interactome$key <- gsub(';.*', '', interactome$prey)

interactome <- merge(interactome, 
                     tmpIps[, c("Gene_name_canonical","graph_localization_annotation")],
                     by.x = 'key', by.y = 'Gene_name_canonical', all.x = T)

tmp <- interactome[, c('bait', 'prey', 'significant',
                       "HCDP", "graph_localization_annotation",
                       "logFC_APMS", "padj_APMS", "logFC_TbID", "padj_TbID" )]
names(tmp)[which(names(tmp) == "graph_localization_annotation")] <- "subcell_localization"

allApmsPreys <- unique(unlist(apmsPreys))
allTidPreys <- unique(unlist(tidPreys))

tmp$method <- ifelse(
  tmp$prey %in% intersect(allApmsPreys, allTidPreys),
  'Both',
  ifelse(tmp$prey %in% allApmsPreys, 'AP-MS', 'TbID')
)

write.table(tmp, file.path(dataDir(), 'processed/erkpathmap_interactome.tsv'),
            row.names = F,sep = '\t',quote = F)
