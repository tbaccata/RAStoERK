#' @author Sebastian Didusch
#' @date 25.06.18
#' 

source(here::here("src/utils/utils_path.R"))

apmsSe <- readRDS(file.path(processedDataDir(), 'apms/se_apms.rds'))
apms_MPLID25011_se <- readRDS(file.path(processedDataDir(), 'apms/se_apms_MPLID25011.rds'))
tidSe <- readRDS(file.path(processedDataDir(), 'turboid/se_turboid.rds'))
tid_MPLID25011_se <- readRDS(file.path(processedDataDir(), 'turboid/se_turboid_MPLID25011.rds'))

background <- unique(
  c(
    gsub(';.*', '', rowData(apmsSe)[isQuantRnames(apmsSe), 'Gene.names'] ),
    gsub(';.*', '', rowData(apmsSe)[isQuantRnames(apms_MPLID25011_se), 'Gene.names'] ),
    gsub(';.*', '', rowData(apmsSe)[isQuantRnames(tidSe), 'Gene.names'] ),
    gsub(';.*', '', rowData(apmsSe)[isQuantRnames(tid_MPLID25011_se), 'Gene.names'] )
  )
)
background <- background[background!='']