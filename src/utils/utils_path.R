library(here)

paths <- list(
  base = here(),
  src = here("src"),
  data = here("data"),
  data_raw = here("data", "raw"),
  data_interim = here("data", "interim"),
  data_processed = here("data", "processed"),
  results = here("results")
)

processedDataDir   <- function(...) here("data", "processed", ...)
figureDir   <- function(...) here("Figures", ...)

rpath   <- function(...) here("results", ...)
dataDir   <- function(...) here("data", ...)
srcDir      <- function(...) here("src", ...)

interimDataDir <- function(...) here("data", "interim", ...)

apmsPreyFile <- file.path(interimDataDir(), 'apms_sn19_multiple_ctrl_preys.rds' )
tidPreyFile <- file.path(interimDataDir(), 'turboid_sn19_multiple_ctrl_preys.rds' )

if (!file.exists(apmsPreyFile)) {
  stop(paste0('Thresholding interactome has not been done yet!\n',
              'Please execute src/pipeline/threshold_interactome.R first.'
              ))
}

source(file.path(srcDir(), "utils/pipeline_utils.R"))
source(file.path(srcDir(), "utils/network_utils.R"))

### data almost always needed
apmsPreys <- readRDS(apmsPreyFile)
apmsDesign <- read.table(file.path(processedDataDir(), 'apms/design_apms.txt' ),
                                      header = T, sep = '\t')
tidPreys <- readRDS(tidPreyFile)
tidDesign <- read.table(file.path(processedDataDir(), 'turboid/design_turboid.txt' ),
                         header = T, sep = '\t')
tidDesign$samples <- as.character(tidDesign$samples)

#
p <- file.path(processedDataDir(), 'erkpathmap_interactome.tsv')

if (!file.exists(p)) {
  p <- file.path(interimDataDir(), 'erkpathmap_interactome.tsv')
} 

interactome <- read.table(
  file.path(processedDataDir(), 'erkpathmap_interactome.tsv'),
  sep = '\t', header = T)

# create output dir
for (idx in seq_along(1:7) ) {
  dir.create(file.path(figureDir(), paste0('Figure', idx )), showWarnings = FALSE)
  dir.create(file.path(figureDir(), paste0('Figure', idx, '/data' )), showWarnings = FALSE )
  dir.create(file.path(figureDir(), paste0('S', idx )), showWarnings = FALSE)
  dir.create(file.path(figureDir(), paste0('S', idx, '/data' )), showWarnings = FALSE )
}


