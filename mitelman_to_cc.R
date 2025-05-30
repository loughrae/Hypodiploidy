setwd("~/PhD/Hypodiploidy")
library(tidyverse)
library(data.table)

select <- dplyr::select

args <- commandArgs(trailingOnly = TRUE)

mit_download <- args[1]
cc_out <- args[2]
cc_err <- args[3]

mit_for_cc <- fread(mit_download, sep = ',') %>% 
  mutate(Inv = paste(Refno, CaseNo, InvNo, sep = '-')) %>% #sample ID
  mutate(Karyotype = ifelse(KaryLong == '', KaryShort, KaryLong)) %>%
#  filter(Sex %in% c('F', 'M')) %>%
  select(Inv, Karyotype) 

setwd('./CytoConverter')
source('cytoscript_vinput.R')
res = CytoConverter(as.matrix(data.frame(mit_for_cc$Inv, mit_for_cc$Karyotype)))
res$Result %>% write.table(cc_out, col.names = T, row.names = F, quote = F, sep = '\t')
res$Error_log %>% as.data.frame() %>% write.table(cc_err, col.names = T, row.names = F, quote = F, sep = '\t')




