library(data.table)
library(tidyverse)

select <- dplyr::select

mit_for_cc <- fread('bq-results-20230925-210156-1695675751704/bq-results-20230925-210156-1695675751704.csv', sep = ',') %>%
  mutate(Inv = paste(RefNo, CaseNo, InvNo, sep = '-')) %>% #sample ID
  mutate(Karyotype = ifelse(KaryLong == '', KaryShort, KaryLong)) %>%
  filter(Sex %in% c('F', 'M')) %>% 
  select(Inv, Karyotype)

setwd('./CytoConverter')
source('cytoscript_vinput.R')
res = CytoConverter(as.matrix(data.frame(mit_for_cc$Inv, mit_for_cc$Karyotype)))
res$Result %>% write.table('../cytoconverter_result.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
res$Error_log %>% as.data.frame() %>% write.table('../cytoconverter_errorlog.tsv', col.names = T, row.names = F, quote = F, sep = '\t')







