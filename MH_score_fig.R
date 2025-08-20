source('setup.R')
source('MH_score_heuristic.R')

tcga_classes <- fread('TCGA_ploidy_classes.tsv')
mitcn_meta <- fread('mitcn_meta.tsv')
mitcn_chr_summary <- fread('mitcn_chr_summary.tsv')


##fig1c: masked hypodiploidy visualisation
viz_ids <- mitcn_meta %>%
  filter(n_clones == 2) %>%
  filter(any_lows > 0, any_hyper > 0) %>%
  arrange(case_id, ploidy) %>%
  dplyr::slice(1:8) %>%
  pull(id)

viz <- mitcn_chr_summary %>%
  filter(id %in% viz_ids) %>%
  mutate(group = ifelse(group == 'Low-Hypodiploid', 'LH', ifelse(group == 'Near-Haploid', 'NH', ifelse(group == 'Hyperdiploid', 'HeH', 'what')))) %>%
  ggplot(aes(x = factor(sub("chr", "", chr), levels = 1:22), y = group, fill = as.factor(chr_somy))) +
  geom_tile() +
  facet_wrap(~case_id, nrow = 4, scales = "free_y") +
  labs(fill = "Somy", x = "Chromosome", y = "") +
  scale_fill_viridis_d(option = "mako") +
  theme_large() +
  theme(legend.position = "right", strip.background = element_blank(), strip.text = element_blank()) +
  labs(subtitle = "(a)")

##fig1d: trisomies with autosome count
heh <- read.xlsx('Supplementary_data_1_Woodward et al.xlsx', startRow = 2) %>%
  filter(!is.na(chr1)) %>%
  pivot_longer(cols = chr1:chrY, names_to = 'Chromosome', values_to = 'CN') %>% 
  select(Chromosome, Modal.no.including.subclonality, Chromosome, CN, Case.no.) %>% 
  filter(!Chromosome %in% c('chrX', 'chrY')) %>%
  filter(CN %in% c('Tri', 'Tetra_2_2', 'Tetra_3_1')) %>% 
  count(Case.no., Modal.no.including.subclonality, CN) %>% 
  pivot_wider(names_from = CN, values_from = n, values_fill = 0) %>% 
  rename(n_3 = Tri) %>%
  mutate(n_4 = Tetra_2_2 + Tetra_3_1) %>% 
  mutate(diff = n_4 - n_3) %>%
  mutate(label = 'True-Neg')

heh_aut <- read.xlsx("Supplementary_data_1_Woodward et al.xlsx", startRow = 2) %>%
  filter(!is.na(chr1)) %>%
  filter(Subclonality.of.whole.chromosomes == "No") %>%
  pivot_longer(cols = chr1:chrY, names_to = "Chromosome", values_to = "CN") %>%
  group_by(Case.no.) %>%
  filter(!Chromosome %in% c("chrX", "chrY")) %>%
  mutate(n_sub = sum(str_count(CN, "Sub"))) %>% #there were still some subclonal chromosomes left
  filter(n_sub == 0) %>%
  ungroup() %>%
  mutate(reCN = case_when(CN %in% c("Di", "UPID") ~ 2, CN %in% c("Tri", "UPT") ~ 3, CN == "Mono" ~ 1, CN == "Penta" ~ 5, CN %in% c("Tetra", "Tetra_3_1", "Tetra_2_2") ~ 4, TRUE ~ 999)) %>%
  group_by(Case.no.) %>%
  summarize(n_autosome = sum(reCN)) %>%
  left_join(heh, by = c("Case.no.")) %>%
  select(label, n_3, n_4, diff, n_autosome) %>%
  rename(chrs = n_autosome)

scatter <- true_pos %>%
  rename(n_3 = n_trisomic, n_4 = n_tetrasomic) %>%
  select(label, n_3, n_4, diff, chrs) %>%
  rbind(heh_aut) %>%
  mutate(label = ifelse(label == "True-Pos", "Masked Hypodiploid", ifelse(label == "True-Neg", "Hyperdiploid", "What"))) %>%
  ggplot(aes(x = chrs, y = n_3, colour = label, group = label)) +
  geom_jitter() +
  stat_cor() +
  ylab("Trisomies") +
  xlab("Autosome Count") +
  labs(subtitle = '(b)', x = 'Autosome Count', y = 'Trisomies') +
  theme_large() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("darkgreen", "darkblue"))

##fig1e: MH score overlap
true_pos_reformat <- true_pos %>%
  select(label, diff, chrs) %>%
  rbind(heh_reformat[, c('label', 'diff', 'chrs')]) %>%
  mutate(label = ifelse(label == "True-Pos", "Masked Hypodiploid", ifelse(label == "True-Neg", "Hyperdiploid", "What")))

overlap <- true_pos_reformat %>%
  ggplot(aes(x = diff, fill = label)) +
  geom_bar(position = "stack") +
  theme_large() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("darkgreen", "darkblue")) + 
  labs(subtitle = "(c)", fill = "", x = "Tetrasomies - Trisomies", y = "Number of Cases")

## fig1f: MH score confusion matrix
confusion <- true_pos_reformat %>% 
  count(label, score = diff > 0) %>%
  group_by(label) %>%
  mutate(tot = sum(n)) %>%
  mutate(perc = n / tot) %>%
  ggplot(aes(x = label, y = score, fill = perc)) +
  geom_tile() +
  geom_label(size = 7, aes(x = label, y = score, colour = perc > 0.5, label = paste0(round(perc * 100, 2), "%\n", n))) +
  xlab("Label") +
  ylab("Tetrasomies > Trisomies") +
  theme_large_classic() +
  labs(subtitle = "(d)") +
  labs(fill = "") +
  theme(legend.position = "none") +
  xlab("") +
  scale_color_manual(values = c("FALSE" = "white", "TRUE" = "black"), guide = "none")


target_query <- GDCquery(project = 'TARGET-ALL-P2', data.category = 'Copy Number Variation', access = 'open', workflow.type = 'ASCAT2', data.type = 'Allele-specific Copy Number Segment')

target <- GDCprepare(target_query) #note I did GDCdownload(query) before this

dedup_aliquots <- target %>% #deal with one duplicated patient
  separate(Sample, into = c("A", "B"), sep = ";") %>%
  mutate(pt = word(A, 1, 3, sep = "-"), collapse = "-") %>%
  distinct(GDC_Aliquot, .keep_all = T) %>%
  distinct(pt, .keep_all = T) %>%
  pull(GDC_Aliquot)

target <- target %>% filter(GDC_Aliquot %in% dedup_aliquots)

rare_tetra <- target %>%
  mutate(loh_ind = Minor_Copy_Number == 0) %>%
  mutate(len = End - Start, prod = Copy_Number*len) %>%
  group_by(GDC_Aliquot) %>%
  summarize(loh = sum(loh_ind*len)/sum(len), ploidy = sum(prod)/sum(len)) %>%
  ggplot(aes(x = ploidy, y = loh)) + geom_point() + stat_cor() +
  theme_large() + labs(x = 'Ploidy', y = 'Loss of heterozygosity') + ggtitle('(a)')


tcga_heuristic <- tcga_scores %>%
  left_join(tcga_classes, by = c("GDC_Aliquot")) %>%
  mutate(interest = ifelse(group != "Other" & wgd == "WGD", "Masked Hypo", ifelse(group == "Other" & n_chr_nosex >= 49 & n_chr_nosex <= 65, "Hyperdiploid", "Other"))) %>%
  filter(interest != "Other") 

tcga_overlap <- tcga_heuristic %>%
  select(interest, n_3, n_4, proj) %>%
  mutate(diff = n_4 - n_3) %>%
  ggplot(aes(x = diff, fill = interest)) +
  geom_bar(position = "dodge") +
  theme_large() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("darkgreen", "darkblue")) +
  xlab("Tetrasomies - Trisomies") +
  labs(fill = "") +
  ylab("Number of cases") + ggtitle('(b)')

enough_masked_hypo <- tcga_classes %>% filter(group != 'Other', wgd == 'WGD') %>% count(proj) %>% filter(n >= 15) %>% pull(proj) 

tcga_heuristic_accuracy <- tcga_heuristic %>%
  filter(proj %in% enough_masked_hypo) %>%
  group_by(proj, interest) %>%
  summarize(pos_rate = mean(diff > 0)) 

tcga_heuristic_accuracy_plot <- tcga_heuristic_accuracy %>%
  ggplot(aes(x = proj, fill = interest, y = pos_rate)) +
  geom_col(position = "dodge") +
  theme_large() +
  theme(legend.position = "bottom") +
  labs(fill = "", x = "Cancer Type", y = "Proportion of Cases with MH Score > 0") +
  coord_flip() +
  scale_fill_manual(values = c("darkgreen", "darkblue")) + ggtitle('(c)')

false_neg_ids <- true_pos %>% filter(diff <= 0) %>% pull(case_id)

false_neg_viz <- mitcn_chr_summary %>%
  filter(case_id %in% false_neg_ids) %>%
  mutate(group = factor(group, levels = c('Near-Haploid', 'Low-Hypodiploid', 'Hyperdiploid', 'Other'))) %>%
  ggplot(aes(x = factor(sub("chr", "", chr), levels = 1:22), y = group, fill = as.factor(chr_somy))) +
  geom_tile() +
  facet_wrap(~case_id, ncol = 1, scales = "free_y") +
  labs(fill = "Somy", x = "Chromosome", y = "") +
  scale_fill_viridis_d(option = "mako") +
  theme_large() +
  theme(legend.position = "right", strip.background = element_blank(), strip.text = element_blank()) +
  ggtitle('(c)')

heh_spec <- heh_aut %>% #no subclonality because I can't count autosomes when that's there
  group_by(chrs) %>%
  mutate(id = row_number()) %>%
  ggplot(aes(x = chrs, y = id, colour = diff > 0)) +
  geom_point() +
  theme_large() +
  ylab("Case") +
  labs(colour = "MH score > 0", subtitle = "Woodward et al (high-hyperdiploid cases, no subclonality, N = 414)", x = "Autosomes") +
  geom_vline(xintercept = 45.5, colour = "darkred") +
  geom_vline(xintercept = 76, colour = "darkred") +
  scale_colour_manual(values = c(`FALSE` = "#CC79A7", `TRUE` = "darkblue")) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(45, 89)) +
  ggtitle('(f)')

mh <- true_pos %>%
  group_by(chrs) %>%
  mutate(id = row_number()) %>%
  ggplot(aes(x = chrs, y = id, colour = diff > 0)) +
  geom_point(size = 2) +
  theme_large() +
  ylab("Case") +
  labs(colour = "MH score > 0", subtitle = "Subclonal masked hypodiploids, N = 126", x = "Autosomes") +
  geom_vline(xintercept = 45.5, colour = "darkred") +
  geom_vline(xintercept = 76, colour = "darkred") +
  scale_colour_manual(values = c(`FALSE` = "#CC79A7", `TRUE` = "darkblue")) +
  theme(legend.position = "bottom") +
  scale_x_continuous(limits = c(45, 89)) +
  ggtitle('(e)')

fully_masked_plot <- scores_mit %>%
  left_join(mitcn_meta, by = c("id", "case_id", "group")) %>%
  filter(diff > 0, n_chr_nosex < 77, any_lows == 0, n_chr_nosex >= 46) %>%
  ggplot(aes(x = n_chr_nosex)) +
  geom_bar(fill = "black") +
  theme_large() +
  xlab("Autosomes") +
  ylab("Inferred Fully-Masked Cases") +
  labs(subtitle = "Inferred Fully-Masked Cases") +
  ggtitle('(d)')

target_classes <- target %>%
  mutate(len = (End - Start) + 1) %>%
  rename(CN = Copy_Number, mCN = Minor_Copy_Number, chr = Chromosome, start = Start, end = End) %>%
  filter(!chr %in% c("chrX", "chrY")) %>%
  mutate(loh_indicator = ifelse(mCN == 0, 1, 0)) %>%
  mutate(prod = CN * len, loh_len = loh_indicator * len) %>%
  group_by(GDC_Aliquot) %>%
  mutate(ploidy = sum(prod) / sum(len)) %>% # calculate mean copy number per sample
  group_by(GDC_Aliquot, chr, ploidy) %>%
  summarize(loh = sum(loh_len) / sum(len), chr_somy = round(sum(prod) / sum(len), 0)) %>% # calculate proportion of LOH per chromosome per sample
  mutate(loss = ifelse(loh > 0.9, "Lost", "Retained")) %>% # this should also cover nullisomy
  mutate(min_somy = ifelse(loss == "Retained", 2, pmin(1, chr_somy))) %>% # base_somy = 2 because excluding the sex chromosomes
  ungroup() %>%
  group_by(GDC_Aliquot, ploidy) %>%
  summarize(n_loh_nosex = sum(loss == "Lost"), n_chr_nosex = sum(chr_somy), min_chr_nosex = sum(min_somy), n_disomic_nosex = sum(chr_somy == 2), n_trisomic = sum(chr_somy == 3), n_tetrasomic = sum(chr_somy == 4), n_nullisomic_nosex = sum(chr_somy == 0)) %>%
  ungroup()

target_spec <- target_classes %>%
  filter(min_chr_nosex >= hypo_threshold) %>%
  mutate(diff = n_tetrasomic - n_trisomic) %>%
  filter(n_chr_nosex > 44) %>%
  group_by(n_chr_nosex) %>%
  mutate(id = row_number()) %>%
  ggplot(aes(x = n_chr_nosex, colour = diff > 0, y = id)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 45.5, colour = "darkred") +
  geom_vline(xintercept = 76, colour = "darkred") +
  theme_large() +
  labs(colour = "MH score > 0", subtitle = "TARGET (non-hypodiploid cases with > 44 autosomes, N = 66)") +
  scale_colour_manual(values = c(`FALSE` = "#CC79A7", `TRUE` = "darkblue")) +
  ylab("Cases") +
  xlab("Autosomes") +
  theme(legend.position = "bottom")

fig2 <- (viz) / (scatter | overlap | confusion) / (mh | (heh_spec / target_spec))
ggsave(plot = fig2, file =  'paper/MH_Score_fig.png', width = 28, height = 30)
ggsave(plot = fig2, file = 'paper/hypo_fig2_MH_Score.pdf', width = 28, height = 30)

figsMH <- (rare_tetra | tcga_overlap) / (false_neg_viz | fully_masked_plot) 
ggsave(plot = figsMH, file =  'paper/Supp_MH_score.png', width = 25, height = 20)
ggsave(plot = figsMH, file = 'paper/hypo_supp_fig1_MH_score.pdf', width = 25, height = 20)

# number of cases analysed by MEDICC
nrow(med) #148 --> 150
nrow(true_pos) #125 --> 126
nrow(heh) #577
true_pos %>% summarise(median(n_trisomic)) #0
heh %>% summarize(median(n_3)) #6

#false-negs
mh_ids <- true_pos %>% pull(case_id) %>% unique()
mitcn_meta %>% filter(case_id %in% mh_ids) %>% count(group) #LHs make up 36% of hypos in the test set
mitcn_meta %>% filter(case_id %in% false_neg_ids) %>% count(group) #but 66% of false negs
true_pos %>% filter(diff <= 0) %>% summarize(median(diff))
true_pos %>% filter(diff <= 0) %>% summarize(median(n_trisomic))
#sensitivity in the hyperdiploid range
true_pos %>% filter(group == 'Hyperdiploid') %>% count(diff > 0)
#rescues
NROW(rescues)
NROW(hyper_rescues)
mitcn_meta %>% filter(id %in% hyper_rescues) %>% filter(n_clones > 1) #one subclonal
coexisting <- mitcn_meta %>% distinct(case_id, .keep_all = T) %>% filter(any_lows > 0, n_clones > any_lows) %>% mutate(doubled = case_id %in% wgd_cases) # coexisting: 141 --> 142
exclusive_hypo <- mitcn_meta %>% distinct(case_id, .keep_all = T) %>% filter(any_lows == n_clones)  #153 --> 159
mitcn_meta %>% filter(id %in% rescues) %>% distinct(case_id) #17
#total hypodiploids: 
nrow(exclusive_hypo) + nrow(coexisting) + NROW(rescues)


#number of hyperdiploids with a subclonal hypo
mitcn_meta %>% filter(group == 'Hyperdiploid', any_lows > 0) %>% distinct(case_id) %>% nrow() #69
#total cases with a hyperdiploid clone --> 989
mitcn_meta %>% filter(group == 'Hyperdiploid') %>% distinct(case_id) %>% nrow()



#the sus rescues 
mitcn_meta %>% filter(id %in% rescues, !id %in% hyper_rescues) 
target_classes %>% filter(n_chr_nosex == 46) %>% count(n_tetrasomic <= n_trisomic)

# clinical info about the viz ID patients
mit_meta <- fread('mit_meta.tsv')
mit_meta %>% filter(id %in% viz_ids)
