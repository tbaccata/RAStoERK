source(here::here("src/utils/utils_path.R"))

mergedQuant <- read.table(file.path(interimDataDir(), 'all_interactions_apms_turboid.tsv'),
                          header = T, sep = '\t')


ttpQuant <- mergedQuant[mergedQuant$Gene.names=='ZFP36', ]
ttpQuant <- ttpQuant[grep('Araf|Braf|Raf1|MEK', ttpQuant$bait), ]
#
idxs <- match(baitOrder[baitOrder%in%ttpQuant$bait], ttpQuant$bait )
ttpQuant <- ttpQuant[idxs, ]
rownames(ttpQuant) <- ttpQuant$bait

#colPalette <- viridis(100) # 
colPalette <- colorRampPalette( brewer.pal(8, 'Blues'))(100)

plotdf <- ttpQuant[, c('logFC_APMS', 'logFC_TbID')]
names(plotdf) <- c('AP-MS', 'TbID')

pdf(file.path(figureDir(), 'Figure7/7A_ttp_log2fc_heatmap_blue.pdf'),
    width = 4.5, height = 6
)
ComplexHeatmap::pheatmap(
  plotdf,
  name = 'log2FC',
  cluster_rows = F,
  cluster_cols = F,
  cellwidth = 15,
  cellheight = 15,
  border_color = 'white',
  color = colPalette
)
dev.off()

write.table(
  ttpQuant[, grep("combined", names(ttpQuant), invert = T )],
  file.path(figureDir(), 'Figure7/data/ttp_log2fc_heatmap.tsv'),
  row.names = F, sep = '\t', quote = F
)

#### qPCR
p <- file.path(dataDir(), 'raw/qpcrs')
resFiles <- list.files(p )

for (fname in resFiles) {
  dat <- read.table(file.path(p, fname), sep = '\t', header = T)
  rownames(dat) <- dat$Gene
  dat <- dat[c('TNF', 'VEGF', 'ZFP36', 'FOS', 'CDKN1A'),]
  #
  prot <- gsub('.tsv', '', fname)
  pdf(file.path(figureDir(), paste0('Figure7/7C_', prot, '_qcpr_heatmap.pdf') ),
      width = 3.5, height = 3.5
  )
  print(ComplexHeatmap::pheatmap(
    main = prot,
    scale = 'row',
    dat[, 2:ncol(dat)],
    #scale = 'row',
    name = 'z-score',
    cluster_rows = F,
    cluster_cols = F,
    cellwidth = 15,
    cellheight = 15,
    border_color = 'white',
    labels_row = dat[, 1],
    color = brewer.pal(8, 'Greys')
  )
  )
  dev.off()
}
