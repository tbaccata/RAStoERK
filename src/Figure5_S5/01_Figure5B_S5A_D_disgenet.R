#' @author Sebastian Didusch
#' @date 25.10.03
#' 
#' plot enrichment from DisGeNET

source(here::here("src/utils/utils_path.R"))

dis2proteome <- read.delim(
  file.path(dataDir(), 'processed/proteome_gene2disease_curated_score_gt_0.3.tsv' ),
  header = T, sep = '\t'
)

tierOverlap <- read.table(
  file.path(figureDir(), 'Figure5/data/activating_tier_overlap.tsv' ),
  header = T, sep = '\t'
)

gfpDisEnrichmentDf <- read.delim(
  file.path(figureDir(), 
            'S5/data/disease_enrichment_interactome.tsv' ),
  header = T, sep = '\t'
)
gfpDisEnrichmentDf$bait <- factor(gfpDisEnrichmentDf$bait, levels = baitOrder )

###
########### Tier overlap
signDf <- tierOverlap[order(tierOverlap$p_adjusted), ]

signDf$diseaseName <- factor(
  signDf$diseaseName,
  levels = rev(signDf$diseaseName[!duplicated(signDf$diseaseName)])
)
signDf$query <- factor(signDf$query,
                       levels = c('RASRAFMEK', 
                                  'RASRAF', 
                                  'RAFMEK',
                                  'RAS', 'RAF', 
                                  'MEK', 'ERK') )

redundantTerms <- c('GIANT PIGMENTED HAIRY NEVUS', 'Pilomyxoid astrocytoma',
                    'Noonan syndrome and Noonan-related syndrome',
                    'Noonan Syndrome 1', 'Epithelioma'
)

tierPlotDf <- signDf[!signDf$diseaseName%in%redundantTerms,]

#### no cancer
##### Combined
cmbDis <- gfpDisEnrichmentDf[gfpDisEnrichmentDf$method=='combined' & 
                               gfpDisEnrichmentDf$type=='Active' &
                               gfpDisEnrichmentDf$diseaseName %in%
                               gfpDisEnrichmentDf$diseaseName[gfpDisEnrichmentDf$p_adjusted < 0.05],]

cancerPattern <- 'epithelioma|tumor|neoplasm|cancer|astrocyt|sarcoma|carcinoma|glioma|blastoma'
cmbDis <- cmbDis[grep(cancerPattern, cmbDis$diseaseName, ignore.case = T, invert = T ),]
cmbDis <- cmbDis[order(cmbDis$p_adjusted), ]

cmbDis$diseaseName <- factor(
  cmbDis$diseaseName, 
  levels = rev(cmbDis$diseaseName[!duplicated(cmbDis$diseaseName)] )
)
cmbDis <- cmbDis[order(cmbDis$p_adjusted), ]
cmbDis$diseaseName <- factor(
  cmbDis$diseaseName, 
  levels = rev(cmbDis$diseaseName[!duplicated(cmbDis$diseaseName)] )
)

rasopathies <- c('leopard', 'noonan', 'rasopath', 
                 'costello', 'cardio-facio-cutaneous', 'nevus',
                 'hydrops fetalis', 'vascular anom',
                 'Congenital arteriovenous malformation')



toFilter <- c('Noonan syndrome and Noonan-related syndrome', 
              "Noonan syndrome-like disorder with loose anagen hair",
              "NOONAN SYNDROME 3",
              'Noonan Syndrome 1', 'GIANT PIGMENTED HAIRY NEVUS',
              'Nevus sebaceous', 'Nevus Sebaceus of Jadassohn', 'Histiocytosis, Langerhans-Cell',
              'Congenital malformation syndromes associated with short stature', 
              'NEVUS, KERATINOCYTIC, NONEPIDERMOLYTIC', 
              'Diffuse Large B-Cell Lymphoma',
              #'Pauciarticular juvenile rheumatoid arthritis', 
              #'Acute polyarticular juvenile rheumatoid arthritis',
              #'Monoarticular juvenile rheumatoid arthritis', 
              'Unipolar Depression', 'Major depression, single episode (disorder)',
              'Trigeminal Neuralgia', 'Sepsis', 'Cardiomyopathy, Familial Hypertrophic, 1', 
              'Papillomatosis', 'Miller-McKusick-Malvaux-Syndrome (3M Syndrome)',
              'Diffuse Large B-Cell Lymphoma',
              'CAPILLARY MALFORMATION-ARTERIOVENOUS MALFORMATION 1',
              'Congenital Lipomatous Overgrowth, Vascular Malformations, and Epidermal Nevi'
)

plotdf <- cmbDis[cmbDis$p_adjusted<0.05,]
paralogPlotDf <- plotdf[!plotdf$diseaseName%in%toFilter,]

newLevls <- c(
  grep(paste(rasopathies, collapse = '|'), levels(paralogPlotDf$diseaseName), value = T, ignore.case = T, invert = T),
  grep(paste(rasopathies, collapse = '|'), levels(paralogPlotDf$diseaseName), value = T, ignore.case = T)
)
paralogPlotDf$diseaseName <- factor(
  paralogPlotDf$diseaseName, levels = newLevls
)

color_range <- range(-log10(c(tierPlotDf$p_adjusted, paralogPlotDf$p_adjusted)) )
size_range <- range(c(tierPlotDf$intersectionSize, paralogPlotDf$intersectionSize)) 

newLevls <- c(
  grep(paste(rasopathies, collapse = '|'), levels(tierPlotDf$diseaseName), value = T, ignore.case = T, invert = T),
  grep(paste(rasopathies, collapse = '|'), levels(tierPlotDf$diseaseName), value = T, ignore.case = T)
)
tierPlotDf$diseaseName <- factor(
  tierPlotDf$diseaseName, levels = newLevls
)

ggplot(tierPlotDf,
       aes(y=diseaseName, 
           x=query,
           size=intersectionSize, 
           color=-log10(p_adjusted))) + 
  geom_point() + 
  scale_size_continuous(limits = size_range) + 
  #scale_color_gradient(limits = color_range, low = "blue", high = "red")
  scale_color_viridis_c(limits = color_range, option = 'viridis') +
  # scale_color_viridis_c(option = 'viridis') +
  ylab('') + xlab('') +
  theme_cowplot() +
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

ggsave(file.path(figureDir(), 'Figure5/5B_activating_tier_overlap_pruned.pdf' ),
       width = 8, height = 5.5
)


ggplot(paralogPlotDf,
       aes(y=diseaseName, 
           x=bait,
           size=intersectionSize, 
           color=-log10(p_adjusted))) + 
  geom_point() + 
  scale_size_continuous(limits = size_range) + 
  #scale_color_gradient(limits = color_range, low = "blue", high = "red")
  scale_color_viridis_c(limits = color_range, option = 'viridis') +
  # scale_color_viridis_c(option = 'viridis') +
  ylab('') + xlab('') +
  theme_cowplot() +
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
ggsave(file.path(figureDir(), 'S5/S5D_disease_enrichment_interactome_noncancer.pdf' ),
       width = 8.5, height = 5.5
)


#### for cytoscape
inflammation <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Inflammation'])
arthritis <- unique(dis2proteome$gene_symbol[grep('arthritis', dis2proteome$disease_name) ])
hypertensive_disease <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Hypertensive disease'])
ischemia <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Myocardial Ischemia'])
rasopathies <- unique(dis2proteome$gene_symbol[grep('noonan|costello|rasopath|Cardio-facio-cutaneous', 
                                                    dis2proteome$disease_name, ignore.case = T) ])
hivInfections <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='HIV Infections'])
hydropsFetalis <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Hydrops Fetalis, Non-Immune'])
hypertrophCardiomyo <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Hypertrophic Cardiomyopathy'])
vascAnom <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Vascular anomaly'])
endometriosis <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Endometriosis'])
nevus <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Epidermal Nevus'])
cholestasis <- unique(dis2proteome$gene_symbol[dis2proteome$disease_name=='Cholestasis, Extrahepatic'])

nrasq61rpathies <- unique(dis2proteome$gene_symbol[grep('Klippel-Trenaunay-Weber|Megalencephaly cutis|Familial cerebral cavernous malformation', 
                                                        dis2proteome$disease_name, ignore.case = T) ])



prey2Seldis <- list(inflammation=inflammation,
                    arthritis=arthritis,
                    hypertensive_disease=hypertensive_disease,
                    ischemia=ischemia,
                    rasopathies=rasopathies,
                    NRASQ61R=nrasq61rpathies,
                    'Hydrops Fetalis' = hydropsFetalis,
                    'Hypertrophic Cardiomyopathy'=hypertrophCardiomyo,
                    'HIV Infections' = hivInfections,
                    "Vascular anomaly" = vascAnom,
                    'Epidermal nevus' = nevus,
                    cholestasis = cholestasis,
                    endometriosis = endometriosis
)
preys <- sort(unique(unlist(prey2Seldis)))

mat <- t(+sapply(prey2Seldis, "%in%", x = preys))  ## matrix output
colnames(mat) <- preys

mat <- t(mat)
mat <- as.data.frame(mat)
mat$Gene.names <- rownames(mat)

### import this table to DISGENET CA network in CytoScape
write.table(mat,
            file.path(figureDir(), 'Figure5/data/5C_selected_diseases.tsv' ),
            row.names = F, sep = '\t', quote = F  )
