source('setup.R')
source('MH_score_heuristic.R')

tcga_classes <- fread('TCGA_ploidy_classes.tsv')
mitcn_meta <- fread('mitcn_meta.tsv')
mitcn_chr_summary <- fread('mitcn_chr_summary.tsv')

print('Making Figure 1')

## fig1a: autosome count distributions for Mitelman and TCGA datasets
autosome_count_mit <- mitcn_meta %>%
  filter(n_chr_nosex < 44) %>%
  ggplot(aes(x = n_chr_nosex)) +
  geom_bar(fill = "black") +
  labs(subtitle = "(a) Acute lymphoblastic leukemia", x = "Autosome Count", y = "Cases") +
  geom_vline(xintercept = 27.5, colour = 'red', linetype = 'dashed') + 
  geom_vline(xintercept = hypo_threshold - 0.5, colour = 'blue', linetype = 'dashed') +
  theme_large_classic()

autosomes_graph <- function(df, title) {
  df %>%
    filter(min_chr_nosex < 44) %>%
    mutate(wgd = factor(wgd, levels = c('WGD', 'No WGD'))) %>%
    ggplot(aes(x = min_chr_nosex, fill = wgd)) +
    geom_bar() +
    scale_x_continuous(limits = c(22, 44), breaks = seq(25, 40, by = 5)) +
    labs(subtitle = title, x = "Inferred Autosome Count", y = "Cases", fill = "") +
    theme_large() +
    geom_vline(xintercept = 27.5, colour = 'darkblue', linetype = 'dashed') +
    geom_vline(xintercept = hypo_threshold - 0.5, colour = 'darkblue', linetype = 'dashed') +
    scale_fill_manual(values = c(`No WGD` = 'black', `WGD` = "lightsteelblue")) + 
    theme_large_classic() + theme(legend.position = 'bottom') 
}

cts_tcga <- tcga_classes %>% autosomes_graph('Pan-TCGA') + theme(legend.position = 'none')
cts_brca <- tcga_classes %>% filter(proj == 'BRCA') %>% autosomes_graph('Breast')
cts_kich <- tcga_classes %>% filter(proj == 'KICH') %>% autosomes_graph('Kidney Chromophobe')
cts_acc <- tcga_classes %>% filter(proj == 'ACC') %>% autosomes_graph('Adrenocortical')
cts_sarc <- tcga_classes %>% filter(proj == 'SARC') %>% autosomes_graph('Sarcoma')

## fig1b: hypodiploid proportions
hypo_props <- tcga_classes %>%
  group_by(proj) %>%
  mutate(n_total = n()) %>%  
  ungroup() %>%
  filter(group != 'Other') %>%
  group_by(proj, group) %>%
  summarise(n = n(), n_total = first(n_total), .groups = 'drop') %>%
  mutate(prop = n / n_total) %>%
  mutate(group = factor(group, levels = c('Low-Hypodiploid', 'Near-Haploid'))) %>%
  ggplot(aes(x = reorder(proj, prop, FUN = sum), y = prop, fill = group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c('darkgreen', 'darkblue')) +
  labs(fill = '', x = 'Cancer Type', y = 'Proportion of All Cases', subtitle = '(b)') +
  theme_large_classic() +
  theme(legend.position = 'bottom') 

class_props <- tcga_classes %>%
  group_by(proj) %>%
  mutate(prop_hypo = mean(group != "Other")) %>%
  mutate(Class = fct_rev(factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Diploid", "Aneuploid", "Polyploid")))) %>%
  ggplot(aes(x = reorder(proj, prop_hypo), fill = Class)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = class_palette) +
  theme_large_classic() +
  theme(axis.text.y = element_text(size = 12)) +
  coord_flip() +
  theme(legend.position = "bottom") +
  labs(x = "", y = "Proportion of Cases", subtitle = '(b)') +
  guides(fill = guide_legend(reverse = TRUE))

##fig1c: masked hypodiploidy visualisation
viz_ids <- mitcn_meta %>%
  filter(n_clones == 2) %>%
  filter(any_lows > 0, any_hyper > 0) %>%
  arrange(case_id, ploidy) %>%
  dplyr::slice(1:8) %>%
  pull(id)

viz <- mitcn_chr_summary %>%
  filter(id %in% viz_ids) %>%
  ggplot(aes(x = factor(sub("chr", "", chr), levels = 1:22), y = group, fill = as.factor(chr_somy))) +
  geom_tile() +
  facet_wrap(~case_id, nrow = 4, scales = "free_y") +
  labs(fill = "Somy", x = "Chromosome", y = "") +
  scale_fill_viridis_d(option = "mako") +
  theme_large() +
  theme(legend.position = "right", strip.background = element_blank(), strip.text = element_blank()) +
  labs(subtitle = "(c)")

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
  labs(subtitle = '(d)', x = 'Autosome Count', y = 'Trisomies') +
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
  labs(subtitle = "(e)", fill = "", x = "Tetrasomies - Trisomies", y = "Number of Cases")

## fig1f: MH score confusion matrix
confusion <- true_pos_reformat %>% 
  count(label, score = diff > 0) %>%
  group_by(label) %>%
  mutate(tot = sum(n)) %>%
  mutate(perc = n / tot) %>%
  ggplot(aes(x = label, y = score, fill = perc)) +
  geom_tile() +
  geom_label(colour = "white", aes(x = label, y = score, label = paste0(round(perc * 100, 2), "%\n", n))) +
  xlab("Label") +
  ylab("Tetrasomies > Trisomies") +
  theme_large_classic() +
  labs(subtitle = "(f)") +
  labs(fill = "") +
  theme(legend.position = "none") +
  xlab("")

fig1a <- (autosome_count_mit / cts_tcga / (cts_brca | cts_kich) / (cts_acc | cts_sarc))

fig1a / (class_props | viz ) / ( scatter | overlap | confusion) + plot_layout(heights = c(1, 1, 1, 1, 2, 2))
ggsave('paper/paper_fig1_2025.png', width = 25, height = 32)
ggsave('paper/paper_fig1_2025.pdf', width = 24, height = 32)


## Supplementary Figure 1

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
  ggtitle('(d)')

fully_masked_plot <- scores_mit %>%
  left_join(mitcn_meta, by = c("id", "case_id", "group")) %>%
  filter(diff > 0, n_chr_nosex < 77, any_lows == 0, n_chr_nosex >= 46) %>%
  ggplot(aes(x = n_chr_nosex)) +
  geom_bar(fill = "black") +
  theme_large() +
  xlab("Autosomes") +
  ylab("Inferred Fully-Masked Cases") +
  labs(subtitle = "Inferred Fully-Masked Cases") +
  ggtitle('(e)')

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

#figs1 <- (rare_tetra | tcga_overlap) / (tcga_heuristic_accuracy_plot | false_neg_viz)

figs1 <- (rare_tetra | tcga_overlap) / (false_neg_viz | mh) / (fully_masked_plot | (heh_spec / target_spec))
ggsave(plot = figs1, file =  'paper/paper_suppfig1_2025.png', width = 25, height = 30)
ggsave(plot = figs1, file = 'paper/paper_suppfig1_2025.pdf', width = 25, height = 30)

## Stats in text ##

# number of ALL patients
mitcn_meta %>% distinct(case_id) %>% nrow() #6907
# number of hypodiploid clones
mitcn_meta %>% filter(group %in% c('Low-Hypodiploid', 'Near-Haploid')) %>% nrow() #305 --> 313
# number of patients with hypodiploid clones
mitcn_meta %>% filter(group %in% c('Low-Hypodiploid', 'Near-Haploid')) %>% distinct(case_id) %>% nrow() #294 --> 301
# Mit LH v NH split
mitcn_meta %>% count(group) #163 NH, 142 LH --> 165 NH, 148 LH
# proportion of hypodiploids between the two modes (Mit)
mitcn_meta %>% filter(group %in% c('Low-Hypodiploid', 'Near-Haploid')) %>% count(n_chr_nosex >= 28 & n_chr_nosex <= 31) %>% mutate(prop = n/sum(n)) #0.052 --> 0.0479
# lowest chr count 
mitcn_meta %>% summarize(min(n_chr_nosex))
mitcn_chr_summary %>% filter(n_chr_nosex == 23) %>% filter(chr_somy != 1) %>% count(chr)
# number of TCGA current v former hypodiploids (note: some Aneuploid cases have n_chr_nosex < hypo_threshold but not enough LOH...)
tcga_classes %>% filter(group != 'Other') %>% count(n_chr_nosex < hypo_threshold)
# LH v NH ratio (TCGA)
tcga_classes %>% filter(group != 'Other') %>% count(Class) %>% mutate(prop = n/sum(n))
# LH frequency by cancer type
tcga_classes %>% group_by(proj) %>% summarize(tot = n(), n_lh = sum(Class == 'Low-Hypodiploid'), prop_lh = n_lh/tot) %>% arrange(desc(prop_lh))
tcga_classes %>% group_by(proj) %>% summarize(tot = n(), n_lh = sum(Class == 'Low-Hypodiploid'), prop_lh = n_lh/tot) %>% arrange(desc(prop_lh)) %>% filter(n_lh >= 5)
# NH frequency by cancer type 
tcga_classes %>% group_by(proj) %>% summarize(tot = n(), n_nh = sum(Class == 'Near-Haploid'), prop_nh = n_nh/tot) %>% arrange(desc(prop_nh)) %>% filter(n_nh > 0)
tcga_classes %>% filter(Class == 'Near-Haploid') %>% count(wgd)
tcga_classes %>% filter(Class == 'Near-Haploid') %>% pull(n_chr_nosex) %>% summary()

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
mitcn_meta %>% filter(id %in% rescues) %>% filter(n_clones > 1) #one subclonal
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

#hypo rates in TCGA blood cancers
tcga_classes %>% group_by(proj) %>% summarize(n_hypos = sum(group != 'Other'), hypo_rate = n_hypos/n()) %>% arrange(hypo_rate) %>% filter(proj %in% c('DLBC', 'LAML'))

# clinical info about the viz ID patients
mit_meta <- fread('mit_meta.tsv')
mit_meta %>% filter(id %in% viz_ids)