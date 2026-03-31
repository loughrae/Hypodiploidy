source('setup.R')
source('MH_score_heuristic.R')

print('Making Figure 3')

tcga_classes <- fread('TCGA_ploidy_classes.tsv')

fasc <- fread("filtered_ascat.tsv") %>%
  filter(GDC_Aliquot %in% tcga_classes$GDC_Aliquot) %>%
  left_join(tcga_classes, by = c("GDC_Aliquot", "proj")) 

mitcn_meta <- fread('mitcn_meta.tsv') 

cnh <- read.csv("cnh_vanDijk2021.csv") %>% rowwise() %>%
  mutate(Patient = paste(str_split(Samplename, pattern = "-")[[1]][1:3], collapse = "-"))


thrip <- read.csv("Rasnic2021_chromothripsis_TCGA.csv", skip = 1) %>% clean_names()

enough_lh <- tcga_classes %>% filter(group == 'Low-Hypodiploid') %>% count(proj) %>% filter(n >= 15) %>% pull(proj) 

pcawg_cn <- fread('~/Downloads/20170119_final_consensus_copynumber_donor')
pcawg_clin <- fread('~/Downloads/pcawg_donor_clinical_August2016_v9.tsv')
pcawg_ploidy <- fread('~/Downloads/consensus.20170217.purity.ploidy_donor')
pcawg_exclusions <- fread('~/Downloads/donor_wgs_exclusion_white_gray')



pcawg_losses <- pcawg_cn %>%
  filter(!is.na(total_cn), !is.na(minor_cn)) %>%
  mutate(len = (end - start) + 1) %>%
  rename(CN = total_cn, mCN = minor_cn) %>%
  mutate(loh_indicator = ifelse(mCN == 0, 1, 0)) %>%
  mutate(prod = CN * len, loh_len = loh_indicator * len) %>%
  group_by(sampleID) %>%
  mutate(ploidy = sum(prod) / sum(len)) %>% # calculate mean copy number per sample
  group_by(sampleID, chr, ploidy) %>%
  summarize(loh = sum(loh_len) / sum(len), chr_somy = round(sum(prod) / sum(len), 0)) %>% # calculate proportion of LOH per chromosome per sample
  mutate(loss = ifelse(loh > 0.9, "Lost", "Retained")) %>% # this should also cover nullisomy
  mutate(min_somy = ifelse(loss == "Retained", 2, pmin(1, chr_somy))) %>% # base_somy = 2 because excluding the sex chromosomes
  ungroup()


pcawg_classes <- pcawg_losses %>%
  filter(chr %in% 1:22) %>%
  left_join(pcawg_ploidy, by = c("sampleID" = "samplename")) %>%
  filter(wgd_uncertain == FALSE) %>%
  left_join(pcawg_exclusions, by = c("sampleID" = "icgc_donor_id")) %>%
  filter(donor_wgs_exclusion_white_gray == "Whitelist") %>%
  group_by(sampleID, ploidy.x, wgd_status) %>%
  summarize(n_loh_nosex = sum(loss == "Lost"), n_chr_nosex = sum(chr_somy), min_chr_nosex = sum(min_somy), n_disomic_nosex = sum(chr_somy == 2), n_nullisomic_nosex = sum(chr_somy == 0)) %>%
  mutate(group = ifelse(min_chr_nosex < 28, "Near-Haploid", ifelse(min_chr_nosex < hypo_threshold, "Low-Hypodiploid", "Other"))) %>%
  mutate(Class = ifelse(group != "Other", group, ifelse(wgd_status == "wgd", "Polyploid",
    ifelse(n_disomic_nosex == 22, "Diploid", "Aneuploid")
  ))) %>% # create 5 ploidy classes
  ungroup() %>%
  filter(n_chr_nosex >= 22, min_chr_nosex >= 22, n_nullisomic_nosex == 0)

pcawg_forsurv <- pcawg_classes %>%
  left_join(pcawg_clin, by = c("sampleID" = "icgc_donor_id")) %>%
  mutate(TCGA = startsWith(submitted_donor_id, "TCGA")) %>%
  filter(TCGA == F) %>%
  filter(donor_vital_status %in% c('deceased', 'alive')) %>%
  mutate(code = as.numeric(donor_vital_status == 'deceased'))


pairwise_comparisons <- list(
  c("Diploid", "Polyploid"),
  c("Diploid", "Low-Hypodiploid"),
  c("Polyploid", "Low-Hypodiploid"),
  c('Low-Hypodiploid', 'Near-Haploid'),
  c('Low-Hypodiploid', 'Aneuploid')
)

## fig3a: WGD rate by hypodiploidy status
mit_wgd_rate <- mitcn_meta %>%
  mutate(group_rescued = ifelse(id %in% rescues & (n_chr_nosex / 2) < 28, 'Near-Haploid', ifelse(id %in% rescues, 'Low-Hypodiploid', group))) %>%
  mutate(high_wgd_clone = ifelse(ploidy > 2.7, "WGD", "Not WGD")) %>%
  group_by(case_id) %>%
  mutate(any_wgd = sum(high_wgd_clone == "WGD" | case_id %in% wgd_cases | id %in% rescues)) %>%
  mutate(any_group = ifelse(sum(group_rescued == "Near-Haploid") > 0, "Near-Haploid", ifelse(sum(group_rescued == "Low-Hypodiploid") > 0, "Low-Hypodiploid", "Other"))) %>%
  distinct(case_id, .keep_all = T) %>%
  mutate(wgd = ifelse(any_wgd > 0, "WGD", "No WGD")) %>%
  mutate(any_group = ifelse(any_group == "Near-Haploid", "NH", ifelse(any_group == "Low-Hypodiploid", "LH", any_group))) %>%
  ungroup() %>%
  mutate(any_group = factor(any_group, levels = c('Other', 'LH', 'NH'))) %>%
  ggplot(aes(x = any_group, fill = wgd)) +
  geom_bar(position = "fill") +
  labs(fill = "", x = "", y = "Proportion of Cases with WGD") +
  theme_large_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(`No WGD` = "black", `WGD` = "lightsteelblue")) +
  ggtitle("(a) ALL") +
  geom_text(data = . %>% count(any_group), aes(x = any_group, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

mit_wgd_rate_hyper_rescues <- mitcn_meta %>%
  mutate(group_rescued = ifelse(id %in% hyper_rescues & (n_chr_nosex / 2) < 28, 'Near-Haploid', ifelse(id %in% hyper_rescues, 'Low-Hypodiploid', group))) %>%
  mutate(high_wgd_clone = ifelse(ploidy > 2.7, "WGD", "Not WGD")) %>%
  group_by(case_id) %>%
  mutate(any_wgd = sum(high_wgd_clone == "WGD" | case_id %in% wgd_cases | id %in% hyper_rescues)) %>%
  mutate(any_group = ifelse(sum(group_rescued == "Near-Haploid") > 0, "Near-Haploid", ifelse(sum(group_rescued == "Low-Hypodiploid") > 0, "Low-Hypodiploid", "Other"))) %>%
  distinct(case_id, .keep_all = T) %>%
  mutate(wgd = ifelse(any_wgd > 0, "WGD", "No WGD")) %>%
  ungroup() %>%
  mutate(any_group = ifelse(any_group == "Near-Haploid", "NH", ifelse(any_group == "Low-Hypodiploid", "LH", any_group))) %>%
  ggplot(aes(x = factor(any_group, levels = c("Other", "LH", "NH")), fill = wgd)) +
  geom_bar(position = "fill") +
  labs(fill = "", x = "", y = "Proportion of Cases with WGD") +
  theme_large_classic() +
  scale_fill_manual(values = c(`No WGD` = "black", `WGD` = "lightsteelblue")) +
  ggtitle("(a)") +
  theme(legend.position = 'bottom') +
  geom_text(data = . %>% count(any_group), aes(x = any_group, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

tcga_wgd_rate <- tcga_classes %>%
  mutate(group = ifelse(group == "Near-Haploid", "NH", ifelse(group == "Low-Hypodiploid", "LH", group))) %>%
  ggplot(aes(x = factor(group, levels = c("Other", "LH", "NH")), fill = wgd)) +
  geom_bar(position = "fill") +
  ggtitle("TCGA") +
  labs(fill = "", x = "", y = "Proportion of Cases with WGD") +
  scale_fill_manual(values = c(`No WGD` = "black", `WGD` = "lightsteelblue")) +
  theme_large_classic() +
  theme(legend.position = "right") +
  geom_text(data = . %>% count(group), aes(x = group, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

## fig3b: WGD rate vs hypodiploidy rate (WGD rate calculated based on non-hypodiploid cases)
wgd_v_hypo <- tcga_classes %>%
  filter(!proj %in% c("ACC", "KICH")) %>% #filter out ACC and KICH (outliers)
  group_by(proj) %>%
  summarise(wgd_rate = sum(wgd == "WGD") / n(), hypo_rate = mean(group == "Low-Hypodiploid"), wgd_rate_nohypo = sum(wgd == "WGD" & group == "Other") / sum(group == "Other")) %>%
  ggplot(aes(x = hypo_rate, y = wgd_rate_nohypo)) +
  geom_point() +
  geom_smooth(size = 1, method = "lm", colour = "black", fill = "lightsteelblue") +
  stat_cor() +
  geom_label_repel(aes(label = proj)) +
  theme_large(base_size = 18) +
  ggtitle('(b)') +
  labs(x = "Low-Hypodiploidy Rate", y = "WGD Rate (non-hypodiploid)") 

tcga_cnh <- tcga_classes %>%
  left_join(cnh, by = c("Patient")) %>%
  filter(!is.na(CNH)) 


short_palette <- c(
  "NH" = "#1b9e77", 
  "LH"    = "#66c2a5", 
  "Polyploid"       = "#3288bd", 
  "Aneuploid"           = "#a6bddb" ,  
  "Diploid" = 'lightblue'
)

## fig3e: CNH by ploidy class (violins)
cnh_violins <- tcga_cnh %>%
  mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>%
  mutate(Class = factor(Class, levels = c("NH", "LH", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ggplot(aes(x = Class, y = CNH)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons_short, method = "wilcox.test",  tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  ggtitle('(e)') +
  labs(x = "", y = "Copy Number Heterogeneity") +
  scale_fill_manual(values = short_palette) +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


tcga_cnh %>% filter(group != 'Other') %>% do(tidy(lm(CNH ~ group + proj, data = .))) #NH remains significantly lower CNH after controlling for proj

#Compute exact p-values for CNH
cnh_data <- tcga_cnh %>%
  mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', 
                        ifelse(Class == 'Near-Haploid', 'NH', Class))) %>%
  mutate(Class = factor(Class, levels = c("NH", "LH", "Diploid", "Aneuploid", "Polyploid")))

cnh_pval_table <- cnh_data %>%
  rstatix::pairwise_wilcox_test(
    CNH ~ Class,
    comparisons = pairwise_comparisons_short,
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


cnh_pval_table_plot <- ggtexttable(
  cnh_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(cnh_pval_table_plot, file = 'paper/tables/cnh_pval_table.pdf', width = 4, height = 2, units = 'in')


cnh_violins_nondoubledhypos <- tcga_cnh %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>%
  mutate(Class = factor(Class, levels = c("NH", "LH", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ggplot(aes(x = Class, y = CNH)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons_short, method = "wilcox.test",  tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = "", y = "Copy Number Heterogeneity") +
  ggtitle('(a)') +
  scale_fill_manual(values = short_palette) + 
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

#Compute exact p-values for CNH non-doubled hypos
cnh_nondoubled_data <- tcga_cnh %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>%
  mutate(Class = factor(Class, levels = c("NH", "LH", "Diploid", "Aneuploid", "Polyploid")))

cnh_nondoubled_pval_table <- cnh_nondoubled_data %>%
  rstatix::pairwise_wilcox_test(
    CNH ~ Class,
    comparisons = pairwise_comparisons_short,
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


cnh_nondoubled_pval_table_plot <- ggtexttable(
  cnh_nondoubled_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(cnh_nondoubled_pval_table_plot, file = 'paper/tables/cnh_nondoubled_pval_table.pdf', width = 4, height = 2, units = 'in')
## fig3f: CNH vs hypo rate
hypo_cnh_CI <- tcga_cnh %>%
  filter(!proj %in% c("KICH", "ACC")) %>%
  group_by(proj) %>%
  summarise(wgd_rate = sum(wgd == "WGD") / n(), hypo_rate = mean(group == "Low-Hypodiploid"), median_cnh = median(CNH[Class %in% c("Diploid", "Aneuploid")])) %>%
  ggplot(aes(x = hypo_rate, y = median_cnh)) +
  geom_point() +
  geom_smooth(size = 1, method = "lm", colour = "black", fill = "lightsteelblue") +
  stat_cor() +
  geom_label_repel(aes(label = proj)) +
  xlab("Low-Hypodiploidy Rate") +
  ylab("Median copy number heterogeneity") +
  theme_large(base_size = 18) +
  ggtitle('(f)')

## stereotyped CNH 
stereotyped <- tcga_cnh %>%
  filter(Class == "Low-Hypodiploid", proj %in% enough_lh) %>%
  mutate(highlight = proj %in% c("ACC", "KICH")) %>%
  ggplot(aes(x = highlight, y = CNH, fill = highlight)) +
  geom_violin() +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(colour = "black", alpha = 0.2, height = 0) +
  stat_compare_means(hjust = -1, method = 'wilcox.test', aes(label = after_stat(paste0('p = ', format(p, scientific = TRUE, digits = 4))))) +
  scale_fill_manual(values = c(`FALSE` = "#D9D9D9", `TRUE` = "#8DA0CB")) +
  theme_large() +
  theme(legend.position = "none") +
  xlab("Stereotyped Hypodiploid") +
  ggtitle('(g)') +
  geom_text(data = . %>% count(highlight), aes(x = highlight, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


## number of copy number segments 
segs_summary <- fasc %>%
  mutate(len = (as.numeric(End) + 1) - as.numeric(Start)) %>%
  group_by(proj, GDC_Aliquot, Class, wgd, group) %>%
  summarize(num_segments = n()) %>%
  mutate(Class = factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ungroup() %>%
  mutate(log10_nsegs = log10(num_segments))

segs_violins <- segs_summary %>%
  ggplot(aes(x = Class, y = log10(num_segments))) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  #scale_y_continuous(trans = "log10") +
  labs(x = "", y = "Number of copy number segments (log10)") + 
  scale_fill_manual(values = class_palette) +
  ggtitle('(c)') +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


#get exact pvals
segnum_pval_table <- segs_summary %>%
  rstatix::pairwise_wilcox_test(
    log10_nsegs ~ Class,
    comparisons = pairwise_comparisons,
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


segnum_pval_table_plot <- ggtexttable(
  segnum_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(segnum_pval_table_plot, file = 'paper/tables/segnum_pval_table.pdf', width = 4, height = 2, units = 'in')


segs_summary %>% filter(group != 'Other') %>% do(tidy(lm(log10(num_segments) ~ group + proj, data = .))) #comparing NH to LH

segs_violins_nondoubledhypos <- segs_summary %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  ggplot(aes(x = Class, y = log10(num_segments))) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test",  tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  #scale_y_continuous(trans = "log10") +
  labs(x = "", y = "Number of copy number segments (log10)") + 
  ggtitle('(d)') +
  scale_fill_manual(values = class_palette) +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

segnum_nondoubled_pval_table <- segs_summary %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  rstatix::pairwise_wilcox_test(
    log10_nsegs ~ Class,
    comparisons = pairwise_comparisons,
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


segnum_nondoubled_pval_table_plot <- ggtexttable(
  segnum_nondoubled_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(segnum_nondoubled_pval_table_plot, file = 'paper/tables/segnum_nondoubled_pval_table.pdf', width = 4, height = 2, units = 'in')


hypo_segs_CI <- fasc %>%
  mutate(len = (as.numeric(End) + 1) - as.numeric(Start)) %>%
  group_by(proj, GDC_Aliquot, Class, wgd, group) %>%
  summarize(num_segments = n()) %>%
  filter(!proj %in% c("KICH", "ACC")) %>%
  group_by(proj) %>%
  summarize(wgd_rate = sum(wgd == "WGD") / n(), hypo_rate = mean(group == "Low-Hypodiploid"), median_nseg_sample_normalploidy = median(num_segments[Class %in% c("Diploid", "Aneuploid")])) %>%
  ggplot(aes(x = hypo_rate, y = median_nseg_sample_normalploidy)) +
  geom_point() +
  geom_smooth(size = 1, method = "lm", colour = "black", fill = 'lightsteelblue') +
  stat_cor() +
  geom_label_repel(aes(label = proj)) +
  xlab("Hypodiploidy Rate") +
  ylab("Median Segment Count") +
  theme_large(base_size = 18) +
  ggtitle('(d)') 

stereotyped_segs <- segs_summary %>% 
  filter(Class == "Low-Hypodiploid", proj %in% enough_lh) %>%
  mutate(highlight = proj %in% c("ACC", "KICH")) %>%
  ggplot(aes(x = highlight, y = num_segments, fill = highlight)) +
  geom_violin() +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(colour = "black", alpha = 0.2, height = 0) +
  stat_compare_means(method = 'wilcox.test', aes(label = after_stat(paste0('p = ', format(p, scientific = TRUE, digits = 4))))) +
  scale_fill_manual(values = c(`FALSE` = "#D9D9D9", `TRUE` = "#8DA0CB")) +
  theme_large() +
  theme(legend.position = "none") +
  ggtitle('(c)') +
  xlab("Stereotyped Hypodiploid") + ylab('Number of CN segments') +
  geom_text(data = . %>% count(highlight), aes(x = highlight, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


## max chromosome proportion
chr_prop_summary <- fasc %>%
  mutate(len = (as.numeric(End) + 1) - as.numeric(Start)) %>%
  group_by(proj, GDC_Aliquot, Chromosome, Class, wgd, group) %>%
  summarize(prop_max = max(len) / sum(len)) %>%
  mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>%
  mutate(Class = factor(Class, levels = c("NH", "LH", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ungroup() 

chr_prop_pval_table <- chr_prop_summary %>%
  rstatix::pairwise_wilcox_test(
    prop_max ~ Class,
    comparisons = pairwise_comparisons_short,
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


chr_prop_pval_table_plot <- ggtexttable(
  chr_prop_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(chr_prop_pval_table_plot, file = 'paper/tables/chr_prop_pval_table.pdf', width = 4, height = 2, units = 'in')

chr_prop_nondoubled_pval_table <- chr_prop_summary %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  rstatix::pairwise_wilcox_test(
    prop_max ~ Class,
    comparisons = pairwise_comparisons_short,
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


chr_prop_nondoubled_pval_table_plot <- ggtexttable(
  chr_prop_nondoubled_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(chr_prop_nondoubled_pval_table_plot, file = 'paper/tables/chr_prop_nondoubled_pval_table.pdf', width = 4, height = 2, units = 'in')

maxprop_violins_func <- function(df, sub) {
  df %>%
    ggplot(aes(x = Class, y = prop_max)) +
    geom_violin(aes(fill = Class)) +
    geom_boxplot(alpha = 0.2) +
    stat_compare_means(comparisons = pairwise_comparisons_short, method = "wilcox.test",  tip.length = 0.02) +
    theme_large(base_size = 18) +
    theme(legend.position = "none") +
    labs(x = '', y = "Proportion of chromosome in longest CN segment") +
    ggtitle(sub) +
    scale_fill_manual(values = short_palette) +
    geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 
  
  
}
maxprop_violins <- chr_prop_summary %>% maxprop_violins_func(sub = '(c)')
  
maxprop_violins_nodoubledhypos <- chr_prop_summary %>%
  filter(group == "Other" | wgd == "No WGD") %>%
 maxprop_violins_func(sub = '(e)')


# Chromothripsis bar 
chromoanagenesis <- tcga_classes %>%
  left_join(thrip, by = c("Patient" = "tcga_id")) %>%
  filter(!is.na(is_chromoanagenesis)) %>%
  mutate(chromothripsis = is_chromoanagenesis == 1) %>%
  mutate(Class = factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ggplot(aes(x = Class, fill = chromothripsis)) +
  geom_bar(position = "fill") +
  labs(fill = "Chromoanagenesis", x = "", y = "Proportion of Cases") +
  theme_large_classic() +
  scale_fill_manual(values = c(`FALSE` = "black", `TRUE` = "lightsteelblue")) +
  theme(legend.position = "bottom") +
  ggtitle('(f)') +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 
  

## Survival Analysis 

clin  <- fread('~/Downloads/TCGA-CDR-SupplementalTableS1(2).csv') #Liu paper, from the PanCanAtlas page (https://gdc.cancer.gov/about-data/publications/pancanatlas)

forsurv <- tcga_classes %>%
  left_join(clin, by = c('Patient' = 'bcr_patient_barcode')) %>%
  mutate(age_at_index = age_at_initial_pathologic_diagnosis) %>% 
  filter(!is.na(PFI.time), !is.na(OS.time)) %>%
  mutate(Class = factor(Class, levels = c('Diploid', 'Near-Haploid', 'Low-Hypodiploid', 'Aneuploid', 'Polyploid'))) 


hypo_effect <- survfit(Surv(OS.time, OS) ~ Class, data = forsurv %>% mutate(Class = as.character(Class)) %>% mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>% mutate(Class = factor(Class, levels = c('Diploid', 'NH', 'LH', 'Aneuploid', 'Polyploid'))))
labels_wrapped <- sub("Class=", "", names(hypo_effect$strata))
surv_ploidy_classes <- ggsurvplot(hypo_effect, data = forsurv, ggtheme = theme_large() + theme(legend.text = element_text(size = 16)), legend.title = '', legend.labs = labels_wrapped)

wgd_effect_hypos <- survfit(Surv(OS.time, OS) ~ wgd, data = forsurv %>% filter(Class == 'Low-Hypodiploid'))
surv_wgd_hypo <- ggsurvplot(wgd_effect_hypos, data = forsurv %>% filter(Class == 'Low-Hypodiploid'), ggtheme = theme_large() + theme(legend.text = element_text(size = 16)), pval = T, palette = c("#E41A1C", "#377EB8"), conf.int = TRUE, legend.title = '', legend.labs = c('No WGD', 'WGD'))

hypo_effect_pfs <- survfit(Surv(PFI.time, PFI) ~ Class, data = forsurv %>% mutate(Class = as.character(Class)) %>% mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>% mutate(Class = factor(Class, levels = c('Diploid', 'NH', 'LH', 'Aneuploid', 'Polyploid'))))
labels_wrapped <- sub("Class=", "", names(hypo_effect_pfs$strata))
surv_ploidy_classes_pfs <- ggsurvplot(hypo_effect_pfs, data = forsurv, ggtheme = theme_large() + theme(legend.text = element_text(size = 16)), legend.title = '', legend.labs = labels_wrapped) + ggtitle('(a)')

# surv: stable vs CIN

forsurv_filtered <- forsurv %>%
  mutate(Class = as.character(Class)) %>%
  filter(Class %in% c("Low-Hypodiploid"), proj %in% enough_lh) %>%
  mutate(typ = ifelse(proj %in% c("KICH", "ACC"), "Stable", "CIN")) 

km_fit <- survfit(Surv(OS.time, OS) ~ typ, data = forsurv_filtered)

surv_stereotyped <- ggsurvplot(
  km_fit,
  data = forsurv_filtered,
  risk.table = FALSE,         
  pval = TRUE,               
  conf.int = TRUE,           
  legend.labs = c("CIN", "Stable"),  
  legend.title = "",
  xlab = "Time",
  ylab = "Survival probability",
  palette = c("#E41A1C", "#377EB8"),
  ggtheme = theme_large() + theme(legend.text = element_text(size = 16))
)


pcawg_hypo_effect <- survfit(Surv(donor_survival_time, code) ~ Class, data = pcawg_forsurv %>% mutate(Class = as.character(Class)) %>% mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>% mutate(Class = factor(Class, levels = c('Diploid', 'NH', 'LH', 'Aneuploid', 'Polyploid'))))
surv_pcawg_classes <- ggsurvplot(pcawg_hypo_effect, data = pcawg_forsurv %>% mutate(Class = as.character(Class)) %>% mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>% mutate(Class = factor(Class, levels = c('Diploid', 'NH', 'LH', 'Aneuploid', 'Polyploid')))) + ggtitle('(k)') 

pcawg_hypo_wgd <- survfit(Surv(donor_survival_time, code) ~ wgd_status, data = pcawg_forsurv %>% filter(Class == 'Low-Hypodiploid'))
pcawg_hypo_wgd_plot <- ggsurvplot(pcawg_hypo_wgd, data = pcawg_forsurv %>% filter(Class == 'Low-Hypodiploid'), pval = T) + ggtitle('(l)') 

coxph(Surv(donor_survival_time, code) ~ wgd_status + donor_age_at_diagnosis + donor_sex + project, data = pcawg_forsurv %>% filter(Class == 'Low-Hypodiploid') %>% separate(project_code, into = c('project', 'country'), sep = '-', remove = F)) %>% tidy()

pcawg_lh_cancer_types <- pcawg_forsurv %>%
  filter(Class == "Low-Hypodiploid") %>%
  separate(project_code, into = c("project", "country"), sep = "-", remove = F) %>%
  ggplot(aes(x = project)) +
  geom_bar() +
  theme_large() +
  ylab("Number of cases") +
  xlab("Cancer Type") +
  labs(subtitle = "PCAWG Low-Hypodiploid Cases") + 
  ggtitle('(j)') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## MH score sensitivity in TCGA

# < 76 autosomes 
mh_tcga_76 <- tcga_scores %>%
  left_join(tcga_classes, by = c("GDC_Aliquot")) %>%
  filter(group != "Other", wgd == "WGD", n_chr_nosex <= 76) %>%
  group_by(proj) %>%
  summarize(n_dh = n(), sens = mean(diff > 0)) %>%
  filter(n_dh >= 15) %>%
  ggplot(aes(x = reorder(proj, sens), y = sens)) +
  geom_col(fill = "black", colour = "black") +
  theme_large() +
  labs(x = "", y = "MH score > 0", subtitle = "Up to 76 autosomes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle("(d)") +
  geom_text(
    aes(y = -0.02, label = paste0("n = ", n_dh)),  # position below x-axis
    vjust = 1, 
    size = 5
  )

# all -- it's more the specificity that changes in tetrasomy-dominant tetraploids, so not relevant to point
mh_tcga_all <- tcga_scores %>% left_join(tcga_classes, by = c('GDC_Aliquot')) %>% filter(group != 'Other', wgd == 'WGD') %>% group_by(proj) %>% summarize(n_dh = n(), sens = mean(diff > 0)) %>% filter(n_dh >= 15) %>% ggplot(aes(x = reorder(proj, sens), y = sens)) + geom_col(fill = 'black', colour = 'black') + theme_large() + labs(x = '', y = 'MH score > 0', subtitle = 'All')  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

remove_axis <- theme(
  axis.title = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank()
)

surv_gg <- as.ggplot(surv_ploidy_classes$plot) + ggtitle('(h)') + theme_large() + xlab('') + theme(panel.border = element_blank()) + remove_axis + guides(color = guide_legend(nrow = 2))
surv_wgd_gg <- as.ggplot(surv_wgd_hypo$plot) + ggtitle('(j)') + theme_large() + xlab('') + theme(panel.border = element_blank()) + remove_axis
surv_stereotyped_gg <- as.ggplot(surv_stereotyped$plot) + ggtitle('(i)') + theme_large() + xlab('') + theme(panel.border = element_blank()) + remove_axis
fig3_hypo_paper <- (((mit_wgd_rate | tcga_wgd_rate) | wgd_v_hypo) + plot_layout(widths = c(1,1,2))) / (segs_violins | hypo_segs_CI) / (cnh_violins | hypo_cnh_CI | stereotyped) / (surv_gg | surv_stereotyped_gg | surv_wgd_gg)
ggsave(plot = fig3_hypo_paper, file = 'paper/fig4_hypo_paper.png', width = 25, height = 30)
ggsave(plot = fig3_hypo_paper, file = 'paper/fig4_CIN_revised.pdf', width = 25, height = 30)

# supp_fig3 <- (mit_wgd_rate_hyper_rescues | maxprop_violins) / (segs_violins_nondoubledhypos | maxprop_violins_nodoubledhypos) / (chromothripsis | cnh_violins_nondoubledhypos) / (stereotyped_segs | (mh_tcga_76))
# ggsave(plot = supp_fig3, file = 'paper/supp_fig3_CIN.png', width = 25, height = 30)
# ggsave(plot = supp_fig3, file = 'paper/hypo_supp_fig3_CIN.pdf', width = 25, height = 30)

### Stats in text ###

# Mit WGD rate: with all rescues
mitcn_meta %>%
  mutate(group_rescued = ifelse(id %in% rescues & (n_chr_nosex / 2) < 28, 'Near-Haploid', ifelse(id %in% rescues, 'Low-Hypodiploid', group))) %>%
  mutate(high_wgd_clone = ifelse(ploidy > 2.7, "WGD", "Not WGD")) %>%
  group_by(case_id) %>%
  mutate(any_wgd = sum(high_wgd_clone == "WGD" | case_id %in% wgd_cases | id %in% rescues)) %>%
  mutate(any_group = ifelse(sum(group_rescued == "Near-Haploid") > 0, "Near-Haploid", ifelse(sum(group_rescued == "Low-Hypodiploid") > 0, "Low-Hypodiploid", "Other"))) %>%
  distinct(case_id, .keep_all = T) %>%
  mutate(wgd = ifelse(any_wgd > 0, "WGD", "No WGD")) %>%
  mutate(any_group = ifelse(any_group == "Near-Haploid", "NH", ifelse(any_group == "Low-Hypodiploid", "LH", any_group))) %>% group_by(any_group) %>% summarize(mean(wgd == 'WGD'))

# Mit WGD rate: with only rescues in the hyperdiploid range
mitcn_meta %>%
  mutate(group_rescued = ifelse(id %in% hyper_rescues & (n_chr_nosex / 2) < 28, 'Near-Haploid', ifelse(id %in% hyper_rescues, 'Low-Hypodiploid', group))) %>%
  mutate(high_wgd_clone = ifelse(ploidy > 2.7, "WGD", "Not WGD")) %>%
  group_by(case_id) %>%
  mutate(any_wgd = sum(high_wgd_clone == "WGD" | case_id %in% wgd_cases | id %in% hyper_rescues)) %>%
  mutate(any_group = ifelse(sum(group_rescued == "Near-Haploid") > 0, "Near-Haploid", ifelse(sum(group_rescued == "Low-Hypodiploid") > 0, "Low-Hypodiploid", "Other"))) %>%
  distinct(case_id, .keep_all = T) %>%
  mutate(wgd = ifelse(any_wgd > 0, "WGD", "No WGD")) %>%
  mutate(any_group = ifelse(any_group == "Near-Haploid", "NH", ifelse(any_group == "Low-Hypodiploid", "LH", any_group))) %>% group_by(any_group) %>% summarize(mean(wgd == 'WGD'))
#note Other WGD rate increases slightly because of two IDs that were in rescues but not hyper_rescues and had ploidy > 2.7 so qualified on that ground

#TCGA WGD rate
tcga_classes %>% group_by(group) %>% summarize(prop_wgd = mean(wgd == 'WGD'))
#chromothripsis rates
tcga_classes %>% left_join(thrip, by = c('Patient' = 'tcga_id')) %>% filter(!is.na(is_chromoanagenesis)) %>%  group_by(Class) %>% summarize(mean(is_chromoanagenesis))
## NH less likely to have chromothripsis than LH
tcga_classes %>% left_join(thrip, by = c('Patient' = 'tcga_id')) %>% filter(!is.na(is_chromoanagenesis)) %>%  filter(group != 'Other') %>% count(Class, is_chromoanagenesis) %>% pivot_wider(values_from = n, names_from = is_chromoanagenesis) %>% column_to_rownames(var = 'Class') %>% as.matrix() %>% fisher.test()
## LH more likely than Dip
tcga_classes %>% left_join(thrip, by = c('Patient' = 'tcga_id')) %>% filter(!is.na(is_chromoanagenesis)) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid')) %>% count(Class, is_chromoanagenesis) %>% pivot_wider(values_from = n, names_from = is_chromoanagenesis) %>% column_to_rownames(var = 'Class') %>% as.matrix() %>% fisher.test()
## LH more likely than aneuploid
tcga_classes %>% left_join(thrip, by = c('Patient' = 'tcga_id')) %>% filter(!is.na(is_chromoanagenesis)) %>% filter(Class %in% c('Low-Hypodiploid', 'Aneuploid')) %>% count(Class, is_chromoanagenesis) %>% pivot_wider(values_from = n, names_from = is_chromoanagenesis) %>% column_to_rownames(var = 'Class') %>% as.matrix() %>% fisher.test()

## stable hypodiploids: num segments
segs_summary %>% filter(Class == 'Low-Hypodiploid') %>% filter(proj %in% enough_lh) %>% group_by(proj %in% c('ACC', 'KICH')) %>% summarize(median(num_segments))
## stable hypodiploids: prop max
chr_prop_summary %>% filter(Class == 'Low-Hypodiploid') %>% filter(proj %in% enough_lh) %>% group_by(proj %in% c('ACC', 'KICH')) %>% summarize(median(prop_max))
## stable hypodiploids: CNH
tcga_cnh %>% filter(Class == 'Low-Hypodiploid') %>% filter(proj %in% enough_lh) %>% group_by(proj %in% c('ACC', 'KICH')) %>% summarize(median(CNH))
## tissue CNH 
tcga_cnh %>% ggplot(aes(x = reorder(proj, CNH, FUN = median), y = CNH)) + geom_boxplot() + theme_large() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
# CNH v Class, controlling for tissue: first against diploid, then against aneuploid
tcga_cnh %>% mutate(Class = factor(Class, levels = c('Diploid', 'Near-Haploid', 'Low-Hypodiploid', 'Aneuploid', 'Polyploid'))) %>% do(tidy(glm(CNH ~ Class + proj, data = .)))
tcga_cnh %>% mutate(Class = factor(Class, levels = c('Aneuploid', 'Near-Haploid', 'Low-Hypodiploid', 'Diploid', 'Polyploid'))) %>% do(tidy(glm(CNH ~ Class + proj, data = .)))
# num segs v Class, controlling for tissue: first against dip, then against aneuploid
segs_summary %>% mutate(Class = factor(Class, levels = c('Diploid', 'Near-Haploid', 'Low-Hypodiploid', 'Aneuploid', 'Polyploid'))) %>% do(tidy(glm(num_segments ~ Class + proj, data = .)))
segs_summary %>% mutate(Class = factor(Class, levels = c('Aneuploid', 'Near-Haploid', 'Low-Hypodiploid', 'Diploid', 'Polyploid'))) %>% do(tidy(glm(num_segments ~ Class + proj, data = .)))



#MH score stability measure in genome doubled tumours: stereotyped vs CIN (the >= 15 is important to make sure it's tumours that were determined to be CIN by chromosome loss patterns)
tcga_scores %>%
  left_join(tcga_classes, by = c("GDC_Aliquot")) %>%
  filter(group != "Other", wgd == "WGD", n_chr_nosex <= 76) %>%
  group_by(proj) %>%
  mutate(n_dh = n()) %>%
  filter(n_dh >= 15) %>%
  ungroup() %>%
  mutate(typ = ifelse(proj %in% c("ACC", "KICH"), proj, "Other")) %>%
  group_by(typ) %>%
  summarize(sens = mean(diff > 0))

## Survival Analysis ##
# LH v Diploid - p = 5.32e-17
coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid'))) %>% tidy()
# LH v Diploid - controlling for cancer type - p = 0.0237
coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid'))) %>% tidy()
# NH v Diploid - p = 0.0251
coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Near-Haploid', 'Diploid'))) %>% tidy()
# LH v Aneuploid - p = 0.00819
coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Aneuploid'))) %>% tidy()
# LH v Aneuploid - controlling for cancer type - p = 0.0179
coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Aneuploid'))) %>% tidy()
# n.s. with polyploid
coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Polyploid'))) %>% tidy()
# doubled v non-doubled LH - n.s. (0.192)
coxph(Surv(OS.time, OS) ~ wgd + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class == 'Low-Hypodiploid')) %>% tidy()
# doubled v non-doubled LH by proj - n.s. (0.588)
coxph(Surv(OS.time, OS) ~ wgd + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class == 'Low-Hypodiploid')) %>% tidy()
# stable v CIN hypodiploids
coxph(Surv(OS.time, OS) ~ typ + age_at_index + gender + race, data = forsurv_filtered) %>% tidy()

### Reviewer Comments ###

# Chromothripsis

pairwise_comparisons_short <- list(
  c("Diploid", "Polyploid"),
  c("Diploid", "LH"),
  c("Polyploid", "LH"),
  c('LH', 'Aneuploid'),
  c('LH', 'NH')
)

chromo_cc <- fread("Chromothripsis_Cortes_Ciriano.csv") %>%
  filter(startsWith(submitted_donor_id, "TCGA")) %>%
  group_by(submitted_donor_id) %>%
  summarize(n_highish = sum(chromo_label %in% c("High confidence", "Linked to high confidence")), n_lowish = sum(chromo_label %in% c("Low confidence", "Linked to low confidence")), n_chromo = sum(chromo_label != "No")) %>%
  left_join(tcga_classes, by = c("submitted_donor_id" = "Patient")) %>%
  filter(!is.na(Class)) %>%
  mutate(Class = ifelse(Class == "Low-Hypodiploid", "LH", ifelse(Class == "Near-Haploid", "NH", Class))) %>%
  mutate(Class = factor(Class, levels = c("NH", "LH", "Diploid", "Aneuploid", "Polyploid")))

chromothripsis_cc_number_plot <- chromo_cc %>%
  ggplot(aes(x = Class, y = n_chromo)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons_short_noNH, method = "wilcox.test",  tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  ggtitle('(d)') +
  labs(x = "", y = "Number of chromothriptic chromosomes") +
  scale_fill_manual(values = short_palette) +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


chromothripsis_cc_number_plot_highconf <- chromo_cc %>%
  ggplot(aes(x = Class, y = n_highish)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons_short_noNH, method = "wilcox.test",  tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  ggtitle('(g)') +
  labs(x = "", y = "Number of chromothriptic chromosomes") +
  scale_fill_manual(values = short_palette) +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

chromo_cc %>%
  ggplot(aes(x = Class, fill = n_highish > 0)) +
  geom_bar(position = 'fill') +
  theme_large(base_size = 18) +
  labs(x = "", y = ">= 1 chromothriptic chromosome") + theme(legend.position = 'bottom')

#presence - bar charts
chromo_cc %>%  mutate(Class = factor(Class, levels = c("NH", "LH", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ggplot(aes(x = Class, fill = n_chromo > 0)) +
  geom_bar(position = 'fill') +
  theme_large(base_size = 18) +
  labs(x = "", y = ">= 1 chromothriptic chromosome") + theme(legend.position = 'bottom')

chromo_cc %>% filter(!is.na(Class)) %>% group_by(Class) %>% count(n_chromo > 0) %>% mutate(prop = n/sum(n))
#proportions


## Controlling LH-CIN associations for tissue and purity
# WGD vs hypo_status, controlling for tissue
tcga_classes %>%
  mutate(big_group = ifelse(group != "Other", "Hypo", group)) %>%
  mutate(WGD_binary = wgd == "WGD") %>%
  mutate(big_group = factor(big_group, levels = c("Other", "Hypo"))) %>%
  do(tidy(glm(WGD_binary ~ big_group + proj, data = .)))

purity <- fread("~/Downloads/Aran2015_purity.csv", skip = 3) %>%
  clean_names() %>%
  mutate(Patient = str_extract(sample_id, "^[^-]+-[^-]+-[^-]+"))

tcga_purity <- tcga_classes %>% left_join(purity, by = c('Patient')) %>% add_count(Patient, name = 'purities_per_patient') %>% filter(purities_per_patient == 1) %>% filter(!is.na(cpe))

# WGD vs hypo_status, controlling for tissue and purity (CPE)
tcga_purity %>%
  mutate(big_group = ifelse(group != "Other", "Hypo", group)) %>%
  mutate(WGD_binary = wgd == "WGD") %>%
  mutate(big_group = factor(big_group, levels = c("Other", "Hypo"))) %>%
  do(tidy(glm(WGD_binary ~ big_group + proj + cpe, data = .)))

# CNH v Class, controlling for tissue and purity: first against diploid, then against aneuploid
tcga_purity %>%
  left_join(cnh, by = c('Patient')) %>%
  filter(!is.na(CNH)) %>%
  mutate(Class = factor(Class, levels = c("Diploid", "Near-Haploid", "Low-Hypodiploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(glm(CNH ~ Class + proj + cpe, data = .)))

tcga_purity %>%
  left_join(cnh, by = c('Patient')) %>%
  filter(!is.na(CNH)) %>%
  mutate(Class = factor(Class, levels = c("Aneuploid", "Near-Haploid", "Low-Hypodiploid", "Diploid", "Polyploid"))) %>%
  do(tidy(glm(CNH ~ Class + proj + cpe, data = .)))

#against LH -- comparing LH to NH
tcga_purity %>%
  left_join(cnh, by = c('Patient')) %>%
  filter(!is.na(CNH)) %>%
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(glm(CNH ~ Class + proj + cpe, data = .)))

# num segs v Class, controlling for tissue and purity: first against dip, then against aneuploid
segs_summary %>%
  left_join(tcga_purity[, c("GDC_Aliquot", "Patient", "cpe")], by = c("GDC_Aliquot")) %>%
  filter(!is.na(Patient)) %>%
  mutate(Class = factor(Class, levels = c("Diploid", "Near-Haploid", "Low-Hypodiploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(glm(num_segments ~ Class + proj + cpe, data = .)))

segs_summary %>%
  left_join(tcga_purity[, c("GDC_Aliquot", "Patient", "cpe")], by = c("GDC_Aliquot")) %>%
  filter(!is.na(Patient)) %>%
  mutate(Class = factor(Class, levels = c("Aneuploid", "Near-Haploid", "Low-Hypodiploid", "Diploid", "Polyploid"))) %>%
  do(tidy(glm(num_segments ~ Class + proj + cpe, data = .)))

#against LH -- comparing LH to NH -- but note there are only 15 NHs in tcga_purity. 
segs_summary %>%
  left_join(tcga_purity[, c("GDC_Aliquot", "Patient", "cpe")], by = c("GDC_Aliquot")) %>%
  filter(!is.na(Patient)) %>%
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(glm(num_segments ~ Class + proj + cpe, data = .)))

#just cancer type
tcga_cnh %>%
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(glm(CNH ~ Class + proj, data = .)))

segs_summary %>%
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(glm(log10(num_segments) ~ Class + proj, data = .)))




# now for non-doubled hypos: CNH

tcga_purity %>%
  left_join(cnh, by = c('Patient')) %>%
  filter(!is.na(CNH)) %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  mutate(Class = factor(Class, levels = c("Diploid", "Near-Haploid", "Low-Hypodiploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(glm(CNH ~ Class + proj + cpe, data = .)))

tcga_purity %>%
  left_join(cnh, by = c('Patient')) %>%
  filter(!is.na(CNH)) %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  mutate(Class = factor(Class, levels = c("Aneuploid", "Near-Haploid", "Low-Hypodiploid", "Diploid", "Polyploid"))) %>%
  do(tidy(glm(CNH ~ Class + proj + cpe, data = .)))

# for non-doubled hypos: number of segments

segs_summary %>%
  left_join(tcga_purity[, c("GDC_Aliquot", "Patient", "cpe")], by = c("GDC_Aliquot")) %>%
  filter(!is.na(Patient)) %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  mutate(Class = factor(Class, levels = c("Diploid", "Near-Haploid", "Low-Hypodiploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(glm(num_segments ~ Class + proj + cpe, data = .)))

segs_summary %>%
  left_join(tcga_purity[, c("GDC_Aliquot", "Patient", "cpe")], by = c("GDC_Aliquot")) %>%
  filter(!is.na(Patient)) %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  mutate(Class = factor(Class, levels = c("Aneuploid", "Near-Haploid", "Low-Hypodiploid", "Diploid", "Polyploid"))) %>%
  do(tidy(glm(num_segments ~ Class + proj + cpe, data = .)))

#stereotyped LHs have lower CNH: controlled for purity
tcga_purity %>%
  left_join(cnh, by = c('Patient')) %>%
  filter(!is.na(CNH)) %>% filter(Class == 'Low-Hypodiploid') %>% filter(proj %in% enough_lh) %>% mutate(highlight = proj %in% c("ACC", "KICH")) %>% do(tidy(glm(CNH ~ highlight + cpe, data = .)))

#stereotyped LHs have lower segment counts: controlled for purity
segs_summary %>%
  left_join(tcga_purity[, c("GDC_Aliquot", "Patient", "cpe")], by = c("GDC_Aliquot")) %>%
  filter(!is.na(Patient)) %>% filter(Class == 'Low-Hypodiploid') %>% filter(proj %in% enough_lh) %>% mutate(highlight = proj %in% c("ACC", "KICH")) %>% do(tidy(glm(num_segments ~ highlight + cpe, data = .)))


## Survival Analysis: PFI ##

coxph(Surv(PFI.time, PFI) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid'))) %>% tidy()
# LH v Diploid - controlling for cancer type 
coxph(Surv(PFI.time, PFI) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid'))) %>% tidy()
# NH v Diploid 
coxph(Surv(PFI.time, PFI) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Near-Haploid', 'Diploid'))) %>% tidy()
# LH v Aneuploid 
coxph(Surv(PFI.time, PFI) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Aneuploid'))) %>% tidy()
# LH v Aneuploid 
coxph(Surv(PFI.time, PFI) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Aneuploid'))) %>% tidy()
# n.s. with polyploid
coxph(Surv(PFI.time, PFI) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Polyploid'))) %>% tidy()
# doubled v non-doubled LH 
coxph(Surv(PFI.time, PFI) ~ wgd + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class == 'Low-Hypodiploid')) %>% tidy()
# doubled v non-doubled LH by proj 
coxph(Surv(PFI.time, PFI) ~ wgd + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class == 'Low-Hypodiploid')) %>% tidy()
# stable v CIN hypodiploids
coxph(Surv(PFI.time, PFI) ~ typ + age_at_index + gender + race, data = forsurv_filtered) %>% tidy()

## Survival Analysis New Plots: new Fig. S4

forsurv2 <- forsurv %>% 
  mutate(stage = ifelse(startsWith(ajcc_pathologic_tumor_stage, 'Stage IV'), 'IV', ifelse(startsWith(ajcc_pathologic_tumor_stage, 'Stage III'), 'III', ifelse(startsWith(ajcc_pathologic_tumor_stage, 'Stage II'), 'II', ifelse(startsWith(ajcc_pathologic_tumor_stage, 'Stage I'), 'I', 'Other'))))) %>%
  filter(stage %in% c('I', 'II', 'III', 'IV')) %>% 
  mutate(Class = as.character(Class)) %>% 
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Diploid", "Aneuploid", "Polyploid")))
                                                 
all_classes_coxph_OS <- coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Diploid", "Aneuploid", "Polyploid")))) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  filter(!startsWith(term, "race"), !startsWith(term, "proj")) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "OS vs Low-Hypodiploid") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3))))  + ggtitle('(b)')

all_classes_coxph_OS_stage <- coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race + proj + stage, data = forsurv2) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "OS vs Low-Hypodiploid (with stage)") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) + ggtitle('(d)')

all_classes_coxph_OS_vdip <- coxph(Surv(OS.time, OS) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = factor(Class, levels = c('Diploid', 'Near-Haploid', 'Low-Hypodiploid', 'Polyploid', 'Aneuploid')))) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "OS vs Diploid (with stage)") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) 


all_classes_coxph_PFS <- coxph(Surv(PFI.time, PFI) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Diploid", "Aneuploid", "Polyploid")))) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  filter(!startsWith(term, "race"), !startsWith(term, "proj")) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "PFS vs Low-Hypodiploid") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) + ggtitle('(c)')

all_classes_coxph_PFS_stage <- coxph(Surv(PFI.time, PFI) ~ Class + age_at_index + gender + race + proj + stage, data = forsurv2) %>%
  tidy(exponentiate = TRUE, conf.int = TRUE) %>%
  filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "PFS vs Low-Hypodiploid (with stage)") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) + ggtitle('(e)')

forsurv_filtered_stage <- forsurv_filtered %>%
  mutate(stage = ifelse(startsWith(ajcc_pathologic_tumor_stage, 'Stage IV'), 'IV', ifelse(startsWith(ajcc_pathologic_tumor_stage, 'Stage III'), 'III', ifelse(startsWith(ajcc_pathologic_tumor_stage, 'Stage II'), 'II', ifelse(startsWith(ajcc_pathologic_tumor_stage, 'Stage I'), 'I', 'Other'))))) %>%
  filter(stage %in% c('I', 'II', 'III', 'IV'))

stereotyped_coxph_os_stage <- coxph(Surv(OS.time, OS) ~ typ + age_at_index + gender + race + stage, data = forsurv_filtered_stage) %>% tidy(exponentiate = TRUE, conf.int = TRUE) %>% filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "OS, Stereotyped vs CIN LH (with stage)") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) + ggtitle('(f)')

stereotyped_coxph_pfs_stage <- coxph(Surv(PFI.time, PFI) ~ typ + age_at_index + gender + race + stage, data = forsurv_filtered_stage) %>% tidy(exponentiate = TRUE, conf.int = TRUE) %>% filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "PFS, Stereotyped vs CIN LH (with stage)") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) + ggtitle('(g)')

hypo_wgd_PFS <- coxph(Surv(PFI.time, PFI) ~ wgd + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class == 'Low-Hypodiploid')) %>% tidy(exponentiate = TRUE, conf.int = TRUE) %>% filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "PFS, WGD in LH") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) 

hypo_wgd_OS <- coxph(Surv(PFI.time, PFI) ~ wgd + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class == 'Low-Hypodiploid')) %>% tidy(exponentiate = TRUE, conf.int = TRUE) %>% 
  filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  filter(term == 'wgdWGD') %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "OS, WGD in LH") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) 

hypo_wgd_PFS_stage <- coxph(Surv(PFI.time, PFI) ~ wgd + age_at_index + gender + race + proj + stage, data = forsurv2 %>% filter(Class == 'Low-Hypodiploid')) %>% tidy(exponentiate = TRUE, conf.int = TRUE) %>% filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "PFS, WGD in LH (with stage)") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) + ggtitle('(i)')


hypo_wgd_OS_stage <- coxph(Surv(OS.time, OS) ~ wgd + age_at_index + gender + race + proj + stage, data = forsurv2 %>% filter(Class == 'Low-Hypodiploid')) %>% tidy(exponentiate = TRUE, conf.int = TRUE) %>% filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "OS, WGD in LH (with stage)") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) + ggtitle('(h)')

hypo_wgd_OS_PCAWG <- coxph(Surv(donor_survival_time, code) ~ wgd_status + donor_age_at_diagnosis + donor_sex + project, data = pcawg_forsurv %>% filter(Class == 'Low-Hypodiploid') %>% separate(project_code, into = c('project', 'country'), sep = '-', remove = F)) %>% tidy(exponentiate = TRUE, conf.int = T) %>% filter(!startsWith(term, "race"), !startsWith(term, "proj"), !startsWith(term, 'stage')) %>%
  filter(term == 'wgd_statuswgd') %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  theme_large() +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  xlab("Term") +
  ylab("Hazard Ratio") +
  labs(subtitle = "OS, WGD in LH - PCAWG") +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) + ggtitle('l')

wgd_effect_hypos_pcawg <- survfit(Surv(donor_survival_time, code) ~ wgd_status, data = pcawg_forsurv %>% filter(Class == 'Low-Hypodiploid'))
surv_wgd_hypo_pcawg <- ggsurvplot(wgd_effect_hypos_pcawg, data = pcawg_forsurv %>% filter(Class == 'Low-Hypodiploid'), ggtheme = theme_large() + theme(legend.text = element_text(size = 16)), pval = T, palette = c("#E41A1C", "#377EB8"), conf.int = TRUE, legend.title = '', legend.labs = c('No WGD', 'WGD')) + ggtitle('(k)')
# 
# suppfig_survival <- (surv_ploidy_classes_pfs$plot | all_classes_coxph_OS | all_classes_coxph_PFS | all_classes_coxph_OS_vdip) / (all_classes_coxph_OS_stage | all_classes_coxph_PFS_stage | stereotyped_coxph_os_stage | stereotyped_coxph_pfs_stage) / (hypo_wgd_OS_stage | hypo_wgd_PFS_stage | hypo_wgd_OS_PCAWG) / (pcawg_hypo_wgd_plot$plot | pcawg_lh_cancer_types)
# ggsave(suppfig_survival, file = 'paper/supp_fig4_CIN_survival.png', width = 20, height = 30)

linreg_cnh <- tcga_cnh %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('Patient')) %>% 
  mutate(Class = factor(Class, levels = c('Low-Hypodiploid', 'Near-Haploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  filter(!is.na(cpe)) %>% do(tidy(lm(CNH ~ Class + proj + cpe, data = .), conf.int = T))  %>% filter(!startsWith(term, 'proj')) %>% ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high, colour = p.value < 0.05)) +
  geom_errorbar() +
  geom_point() +
  coord_flip() +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) +
  theme_large() + xlab('Term') + theme(legend.position = 'bottom') + ggtitle('(b) CNH vs LH') + scale_colour_manual(values = c(`FALSE` = 'grey50', `TRUE` = 'black'))

linreg_segs <- segs_summary %>%
  left_join(tcga_classes[, c("Patient", "GDC_Aliquot")], by = c("GDC_Aliquot")) %>%
  left_join(tcga_purity[, c("Patient", "cpe")], by = c("Patient")) %>%
  filter(!is.na(cpe)) %>%
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Aneuploid", "Diploid", "Polyploid"))) %>%
  do(tidy(lm(log10(num_segments) ~ Class + proj + cpe, data = .), conf.int = T)) %>%
  filter(term != '(Intercept)', !startsWith(term, "proj")) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high, colour = p.value < 0.05)) +
  geom_errorbar() +
  geom_point() +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) +
  coord_flip() +
  theme_large() +
  xlab("Term") +
  theme(legend.position = "bottom") +
  ggtitle("(b) Segment Count (log10) vs Low-Hypodiploid") + scale_colour_manual(values = c('grey50', 'black'))

#segs: LR adjusted for aneuploidy level (n_disomic)
segs_summary %>%
  left_join(tcga_classes[, c("Patient", "GDC_Aliquot", "n_disomic_nosex")], by = c("GDC_Aliquot")) %>%
  left_join(tcga_purity[, c("Patient", "cpe")], by = c("Patient")) %>%
  filter(!is.na(cpe)) %>%
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Aneuploid", "Diploid", "Polyploid"))) %>%
  mutate(LH = Class == 'Low-Hypodiploid') %>%
  do(tidy(lm(log10(num_segments) ~ LH + proj + cpe + n_disomic_nosex, data = .), conf.int = T)) %>%
  filter(term != '(Intercept)', !startsWith(term, "proj"))

#plot
segs_adj_aneuploidy_plot <- segs_summary %>%
  left_join(tcga_classes[, c("Patient", "GDC_Aliquot", "n_disomic_nosex")], by = c("GDC_Aliquot")) %>%
  left_join(tcga_purity[, c("Patient", "cpe")], by = c("Patient")) %>%
  filter(!is.na(cpe)) %>%
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Aneuploid", "Diploid", "Polyploid"))) %>%
  mutate(LH = Class == 'Low-Hypodiploid') %>%
  do(tidy(lm(log10(num_segments) ~ LH + proj + cpe + n_disomic_nosex, data = .), conf.int = T)) %>%
  filter(term != '(Intercept)', !startsWith(term, "proj")) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) +
  coord_flip() +
  theme_large() +
  xlab("Term") +
  theme(legend.position = "bottom") +
  ggtitle("(f) Segment Count (log10) vs LH") + geom_hline(yintercept = 0, colour = 'red')

#CNH: LR adjusted for aneuploidy level (n_disomic)
tcga_cnh %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('Patient')) %>% 
  mutate(Class = factor(Class, levels = c('Low-Hypodiploid', 'Near-Haploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  mutate(LH = Class == 'Low-Hypodiploid') %>%
  filter(!is.na(cpe)) %>% do(tidy(lm(CNH ~ LH + proj + cpe + n_disomic_nosex, data = .), conf.int = T)) 
  
  

cnh_adj_aneuploidy_plot <- tcga_cnh %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('Patient')) %>% 
  mutate(Class = factor(Class, levels = c('Low-Hypodiploid', 'Near-Haploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  mutate(LH = Class == 'Low-Hypodiploid') %>%
  filter(!is.na(cpe)) %>% do(tidy(lm(CNH ~ LH + proj + cpe + n_disomic_nosex, data = .), conf.int = T)) %>%
  filter(term != '(Intercept)', !startsWith(term, "proj")) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) +
  coord_flip() +
  theme_large() +
  xlab("Term") +
  theme(legend.position = "bottom") +
  ggtitle("(e) CNH vs LH") + geom_hline(yintercept = 0, colour = 'red')


#stereotyped: CNH adjusted for purity
tcga_cnh %>% filter(proj %in% enough_lh) %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('Patient')) %>% filter(!is.na(cpe)) %>% mutate(highlight = proj %in% c('ACC', 'KICH')) %>% do(tidy(lm(CNH ~ highlight + cpe, data = .), conf.int = T))  
#stereotyped: num_segs adjusted for purity
segs_summary %>%  left_join(tcga_classes[, c("Patient", "GDC_Aliquot")], by = c("GDC_Aliquot")) %>% filter(proj %in% enough_lh) %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('Patient')) %>% filter(!is.na(cpe)) %>% mutate(highlight = proj %in% c('ACC', 'KICH')) %>% do(tidy(lm(log10(num_segments) ~ highlight + cpe, data = .), conf.int = T))  


#new Fig. S5: survival

survival_fig <- (surv_ploidy_classes_pfs$plot | all_classes_coxph_OS | all_classes_coxph_PFS) / (all_classes_coxph_OS_stage | all_classes_coxph_PFS_stage | stereotyped_coxph_os_stage) / (stereotyped_coxph_pfs_stage | hypo_wgd_OS_stage | hypo_wgd_PFS_stage) / (pcawg_lh_cancer_types | (surv_pcawg_classes$plot  + theme_large()) | (pcawg_hypo_wgd_plot$plot + theme_large()))
ggsave(plot = survival_fig, file = 'paper/supp_CIN_survival_fig_revised.png', width = 25, height = 30)
ggsave(plot = survival_fig, file = 'paper/supp_fig6_CIN_survival_revised.pdf', width = 25, height = 30)

new_supp_fig_CIN <- (mit_wgd_rate_hyper_rescues | linreg_segs | maxprop_violins) / (segs_violins_nondoubledhypos | maxprop_violins_nodoubledhypos) / (chromoanagenesis | chromothripsis_cc_number_plot_highconf)
ggsave(plot = new_supp_fig_CIN, file = 'paper/supp_fig4_CIN_scales_revised.pdf', width = 25, height = 30)

new_supp_fig_CNH <- (cnh_violins_nondoubledhypos | linreg_cnh) / (stereotyped_segs | (mh_tcga_76)) / (cnh_adj_aneuploidy_plot | segs_adj_aneuploidy_plot)
ggsave(plot = new_supp_fig_CNH, file = 'paper/supp_fig5_CIN_CNH_revised.pdf', width = 15, height = 20)
