source('setup.R')

mit_dips <- fread('mit_diploids.tsv', header = F) %>% pull(V1)
mit_meta <- fread('mit_meta.tsv') 

mitcn <- fread('all_mitelman_CNs.bed', header = F) %>%
  mutate(id = V6) %>%
  rename(chr = V1, start = V2, end = V3, netCN = V4, collated = V5) %>%
  filter(!(id %in% mit_dips & chr == 'chr1' & start == 0 & end == 1)) %>% #because empty.bed is chr1-0-1!
  mutate(start = ifelse(id %in% mit_dips & start == 2, 0, start)) %>% #i could say if chr == chr1 but shouldn't have to
  separate(id, into = c('Ref', 'Case', 'Inv_Clone'), sep = '-', remove = F) %>%
  separate(Inv_Clone, into = c('Sample', 'Subclone')) %>%
  mutate(patient_id = paste(Ref, Case, sep = '-')) %>%
  mutate(case_id = paste(Ref, Case, Sample, sep = '-')) %>%
  left_join(mit_meta, by = c('id')) %>%
  filter(!chr %in% c('chrX', 'chrY')) %>%
  mutate(CN = ifelse(netCN == -999, 2, 2 + netCN)) %>%
  group_by(id) %>%
  mutate(neg_CN = sum(CN < 0)) %>%
  filter(neg_CN == 0) %>%
  ungroup()

#check for cases where the same segment is both gained and lost
contrasting_alterations <- mitcn %>%
  separate_rows(collated, sep = ',') %>%
  group_by(chr, start, end, id) %>%
  mutate(opposites = n_distinct(collated)) %>%
  group_by(id) %>%
  summarize(uhoh = sum(opposites != 1)) %>%
  filter(uhoh > 0) %>%
  pull(id) %>%
  unique()


mitcn_prep <- mitcn %>%
  filter(!id %in% contrasting_alterations) %>% 
  group_by(patient_id) %>%
  mutate(samples_per_patient = n_distinct(Inv)) %>%
  filter(samples_per_patient == 1) %>%
  ungroup() %>% 
  mutate(len = end - start) %>%
  group_by(id, case_id) %>%
  mutate(ploidy = sum(CN*len)/sum(len)) %>%
  group_by(id, case_id, chr) %>% 
  mutate(chr_len = sum(len)) %>%
  group_by(id, case_id, chr) %>%
  mutate(loss_perc = sum(len[CN < 2])/sum(len)) %>%
  mutate(loss = ifelse(loss_perc > 0.9, 'Lost', 'Retained')) %>%
  ungroup()
  
#chromosome level summary
mitcn_chr_summary <- mitcn_prep %>%
  group_by(id, case_id, chr, chr_len, ploidy, Subclone, Sex, loss, loss_perc) %>%
  summarize(chr_ploidy = sum(CN*len)/sum(len)) %>% 
  mutate(chr_somy = round(chr_ploidy, 0)) %>% 
  group_by(id, case_id) %>%
  mutate(n_chr_nosex = sum(chr_somy), n_nullisomic = sum(chr_somy == 0)) %>% 
  mutate(n_trisomic = sum(chr_somy == 3), n_tetrasomic = sum(chr_somy == 4), n_disomic = sum(chr_somy == 2)) %>%
  filter(n_nullisomic == 0, n_chr_nosex >= 22) %>% 
  mutate(group = ifelse(n_chr_nosex < 28, 'Near-Haploid', ifelse(n_chr_nosex < hypo_threshold, 'Low-Hypodiploid', ifelse(n_chr_nosex >= 49 & n_chr_nosex <= 65, 'Hyperdiploid', 'Other')))) %>%
  mutate(dipl = id %in% mit_dips) %>%
  group_by(case_id) %>%
  mutate(any_lows = n_distinct(id[group %in% c('Near-Haploid', 'Low-Hypodiploid')])) %>%
  mutate(any_hyper = n_distinct(id[group == 'Hyperdiploid'])) %>%
  mutate(n_clones = n_distinct(id)) 

mitcn_chr_summary %>% fwrite('mitcn_chr_summary.tsv', quote = F, sep = '\t', col.names = T, row.names = F)

mitcn_meta <- mitcn_chr_summary %>%
  distinct(id, .keep_all = T) %>%
  select(id, case_id, Subclone, n_chr_nosex, group, any_lows, any_hyper, n_clones, chr_len, Sex, ploidy, dipl, n_nullisomic, n_disomic, n_trisomic, n_tetrasomic) %>%
  ungroup()

mitcn_meta %>% fwrite('mitcn_meta.tsv', quote = F, sep = '\t', col.names = T, row.names = F)

mitcn_prep <- mitcn_prep %>% filter(id %in% mitcn_meta$id)

mitcn_prep %>% fwrite('mitcn_prep.tsv', quote = F, sep = '\t', col.names = T, row.names = F)

mitcn_prep %>% 
  select(chr, start, end, id, CN) %>% #no need to do start - 1
  write.table('mitcn_all_for_arms.bed', sep = '\t', quote = F, col.names = F, row.names = F)

mitcn_prep %>%
  left_join(mitcn_meta, by = c('case_id', 'id')) %>%
  filter(n_clones > 1, any_lows > 0) %>% 
  select(chr, start, end, id, CN, case_id) %>%
  arrange(chr, start, end, id) %>%
  write.table('mitALL_preMEDICC.bed', row.names = F, col.names = F, sep = '\t', quote = F)