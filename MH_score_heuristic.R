source('setup.R')
fasc <- fread('filtered_ascat.tsv')
mitcn_prep <- fread('mitcn_prep.tsv')
tcga_classes <- fread('TCGA_ploidy_classes.tsv')

print('MH score heuristic')

wgd_mits <- read.table("wgd_statuses.tsv", sep = "\t") %>%
  mutate(V1 = sub("_summary.tsv:wgd_status", "", V1)) %>%
  mutate(V1 = sub("MITMED/", "", V1)) %>%
  rename(sample = V1, WGD = V2)

nsamps <- read.table("nsamples.tsv", sep = "\t") %>%
  mutate(V1 = sub("_summary.tsv:nsamples", "", V1)) %>%
  mutate(V1 = sub("MITMED/", "", V1)) %>%
  rename(sample = V1, nsampl = V2)

med <- wgd_mits %>% left_join(nsamps, by = c("sample"))

wgd_cases <- med %>%
  filter(WGD != "no WGD") %>%
  pull(sample) 

med_single_ids <- med %>%
  filter(nsampl == 3) %>% #this actually means two clones, because of the diploid pseudoroot
  filter(WGD != "no WGD") %>%
  mutate(WGD = sub("WGD on branch ", "", WGD)) %>%
  pull(WGD)

## Calculate MH score for mit
scores_mit <- mitcn_prep %>%
  group_by(id, case_id, group, chr, chr_len, CN, patient_id, Ref, Case) %>%
  summarize(len_cn = sum(len)) %>%
  mutate(frac_cn = len_cn / chr_len) %>%
  group_by(id, case_id, group) %>%
  filter(frac_cn >= 0.9) %>%
  mutate(valid_chromosomes = n_distinct(chr)) %>%
  filter(valid_chromosomes >= 10) %>%
  mutate(CN = paste0("n_", CN), CN = ifelse(CN %in% c("n_1", "n_2", "n_3", "n_4"), CN, "Other")) %>%
  count(id, case_id, group, CN) %>%
  pivot_wider(names_from = CN, values_from = n, values_fill = 0) %>%
  mutate(diff = n_4 - n_3) %>%
  ungroup()


# true positives:
true_pos <- mitcn_prep %>%
  select(id, patient_id, case_id, n_clones, any_lows, any_hyper, group, n_chr_nosex) %>%
  filter(n_clones > 1, any_lows > 0, n_clones > any_lows, n_chr_nosex < 78) %>% # include hyperdiploids and near-triploids
  distinct(id, .keep_all = T) %>%
  filter(id %in% med_single_ids) %>%
  left_join(scores_mit, by = c("id", "case_id", "group")) %>%
  mutate(label = "True-Pos") %>%
  mutate(chrs = n_chr_nosex)

# true negatives: hyperdiploid cases without a hypodiploid history (determined by LOH)
heh <- read.xlsx("Supplementary_data_1_Woodward et al.xlsx", startRow = 2) %>%
  filter(!is.na(chr1)) %>%
  pivot_longer(cols = chr1:chrY, names_to = "Chromosome", values_to = "CN") %>%
  select(Chromosome, Modal.no.including.subclonality, Chromosome, CN, Case.no.) %>%
  filter(!Chromosome %in% c("chrX", "chrY")) %>%
  filter(CN %in% c("Tri", "Tetra_2_2", "Tetra_3_1")) %>%
  count(Case.no., Modal.no.including.subclonality, CN) %>%
  pivot_wider(names_from = CN, values_from = n, values_fill = 0) %>%
  rename(n_3 = Tri) %>%
  mutate(n_4 = Tetra_2_2 + Tetra_3_1) %>%
  mutate(diff = n_4 - n_3) %>%
  mutate(label = "True-Neg")

heh_reformat <- heh %>%
  rename(chrs = Modal.no.including.subclonality) %>%
  select(label, n_3, n_4, diff, chrs)

heh %>% count(n_3 < 3) # 575 correct, 2 wrong - changes with the nosex sigh.
heh %>% count(diff > 0) # 576 correct, 1 wrong

rescues <- scores_mit %>%
  left_join(mitcn_prep, by = c("id", "case_id", "group")) %>%
  distinct(id, .keep_all = T) %>%
  filter(diff > 0, n_chr_nosex < 78, any_lows == 0, n_chr_nosex >= 49) %>%
  pull(case_id)

hyper_rescues <- scores_mit %>%
  left_join(mitcn_prep, by = c("id", "case_id", "group")) %>%
  distinct(id, .keep_all = T) %>%
  filter(diff > 0, n_chr_nosex < 78, group == "Hyperdiploid", any_lows == 0, n_chr_nosex >= 49) %>%
  pull(case_id)

tcga_scores <- fasc %>% 
  mutate(len = (End - Start) + 1) %>%
  filter(GDC_Aliquot %in% tcga_classes$GDC_Aliquot) %>% 
  group_by(GDC_Aliquot, Chromosome) %>% 
  mutate(chr_len = sum(len)) %>% 
  group_by(GDC_Aliquot, Chromosome, chr_len, Copy_Number) %>% 
  summarize(len_cn = sum(len)) %>%
  mutate(frac_cn = len_cn/chr_len) %>% 
  group_by(GDC_Aliquot) %>% 
  filter(frac_cn >= 0.9) %>% 
  mutate(valid_chromosomes = n_distinct(Chromosome)) %>% 
  filter(valid_chromosomes >= 10) %>%
  mutate(CN = paste0('n_', Copy_Number), CN = ifelse(CN %in% c('n_1', 'n_2', 'n_3', 'n_4'), CN, 'Other')) %>% 
  count(GDC_Aliquot, CN) %>% 
  pivot_wider(names_from = CN, values_from = n, values_fill = 0) %>% 
  mutate(diff = n_4 - n_3) %>%
  ungroup()

