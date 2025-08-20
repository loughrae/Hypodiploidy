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
  labs(subtitle = "(a) Acute lymphoblastic leukemia (ALL)", x = "Autosome Count", y = "Cases") +
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
cts_brca <- tcga_classes %>% filter(proj == 'BRCA') %>% autosomes_graph('Breast (BRCA)')
cts_kich <- tcga_classes %>% filter(proj == 'KICH') %>% autosomes_graph('Kidney Chromophobe (KICH)')
cts_acc <- tcga_classes %>% filter(proj == 'ACC') %>% autosomes_graph('Adrenocortical (ACC)')
cts_sarc <- tcga_classes %>% filter(proj == 'SARC') %>% autosomes_graph('Sarcoma (SARC)')

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

class_props_facets <- tcga_classes %>%
  group_by(proj) %>%
  mutate(prop_hypo = mean(group != "Other")) %>%
  mutate(set = ifelse(prop_hypo > 0.5, 'Hypodiploid', ifelse(prop_hypo > 0.05, '≥5% Hypodiploid', '<5% Hypodiploid'))) %>%
  mutate(set = factor(set, levels = c('Hypodiploid', '≥5% Hypodiploid', '<5% Hypodiploid'))) %>%
  mutate(Class = fct_rev(factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Diploid", "Aneuploid", "Polyploid")))) %>%
  ggplot(aes(x = reorder(proj, -prop_hypo), fill = Class)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = class_palette) +
  theme_large_classic() +
  theme(axis.text.y = element_text(size = 12)) +
  theme(legend.position = "bottom") +
  facet_wrap(~set, scales = 'free_x') +
  labs(x = "", y = "Proportion of Cases", subtitle = '(b)') +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    strip.background = element_blank(),  
    legend.position = "bottom"
  )

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
#the LH peak
#ALL peaks
mitcn_meta %>% count(n_chr_nosex >= 23 & n_chr_nosex <= 27)
mitcn_meta %>% count(n_chr_nosex >= 32 & n_chr_nosex <= 37) #129
mitcn_meta %>% count(n_chr_nosex >= 28 & n_chr_nosex <= 31)
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


#hypo rates in TCGA blood cancers
tcga_classes %>% group_by(proj) %>% summarize(n_hypos = sum(group != 'Other'), hypo_rate = n_hypos/n()) %>% arrange(hypo_rate) %>% filter(proj %in% c('DLBC', 'LAML'))

