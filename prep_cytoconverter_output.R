source('setup.R')

print('Filtering CytoConverter output')

#process Mitelman database download
mit_all <- fread('bq-results-20230925-210156-1695675751704/bq-results-20230925-210156-1695675751704.csv', sep = ',')  %>%
  mutate(Inv = paste(RefNo, CaseNo, InvNo, sep = '-')) %>% #sample ID
  mutate(Karyotype = ifelse(KaryLong == '', KaryShort, KaryLong)) %>%
  select(RefNo:InvNo, Inv, Karyotype, Clones, Sex, Age) %>% #tidy up
  filter(!is.na(Sex)) %>%
  separate_rows(Karyotype, sep = '/') %>% #separate samples into clones
  group_by(Inv) %>%
  mutate(clone = row_number()) %>%
  mutate(id = paste(Inv, clone, sep = '_')) %>% #clone ID
  separate(Karyotype, into = c('First', NA), extra = 'merge', sep = ',', remove = F) %>% #extract chromosome count
  mutate(ranges = str_count(First, '-')) %>%
  ungroup()

## identify samples with Y chromosomes mentioned or with chromosome count ranges
Y_inv <- mit_all %>% mutate(Ymentions = str_count(Karyotype, 'Y')) %>% filter(Ymentions != 0) %>% pull(Inv)
range_samples <- mit_all %>% filter(ranges > 0) %>% pull(Inv)

## make first metadata table
mit_meta <- mit_all %>% distinct(id, .keep_all = T) %>% select(Inv, id, Sex, Age) %>% mutate(Y = Inv %in% Y_inv)
mit_meta %>% fwrite('mit_meta.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

## filter cytoconverter output to remove clones with errors, zero-length segments, and samples with ranges, and change Gain to +1 and Loss to -1
cc_errors <- read.table('cytoconverter_errorlog.tsv', header = T, sep = '\t')
all_cc <- read.table('cytoconverter_result.tsv', sep = '\t', header = T) %>% 
  clean_names() %>%
  separate(sample_id, into = c('Ref', 'Case', 'Sample'), sep = '-', remove = F) %>% 
  separate(Sample, into = c('Sam', 'Subclone'), sep = '_') %>% 
  mutate(Inv = paste(Ref, Case, Sam, sep = '-')) %>%
  filter(!sample_id %in% cc_errors$Sample.ID) %>% #so including ones where one of the clones had an issue? hmmm...
  filter(!Inv %in% range_samples) %>%
  mutate(change = ifelse(type == 'Gain', 1, -1)) %>%
  filter(start != end) 

all_cc %>%
  select(chr, start, end, sample_id, change) %>%
  arrange(chr, start, end) %>%
  write.table('all_cytoconverted.bed', sep = '\t', col.names = F, row.names = F, quote = F)

mit_dips <- mit_all %>% mutate(pluses = str_count(Karyotype, '\\+')) %>% filter(!id %in% all_cc$sample_id, !id %in% cc_errors$Sample.ID, First == 46, pluses == 0) %>% pull(id)
mit_dips %>% write.table(file = 'mit_diploids.tsv', sep = '\t', quote = F, row.names = F, col.names = F)
