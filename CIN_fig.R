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
  ggplot(aes(x = factor(any_group, levels = c("Other", "LH", "NH")), fill = wgd)) +
  geom_bar(position = "fill") +
  labs(fill = "", x = "", y = "Proportion of Cases with WGD") +
  theme_large_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(`No WGD` = "black", `WGD` = "lightsteelblue")) +
  ggtitle("(a) ALL") 

mit_wgd_rate_hyper_rescues <- mitcn_meta %>%
  mutate(group_rescued = ifelse(id %in% hyper_rescues & (n_chr_nosex / 2) < 28, 'Near-Haploid', ifelse(id %in% hyper_rescues, 'Low-Hypodiploid', group))) %>%
  mutate(high_wgd_clone = ifelse(ploidy > 2.7, "WGD", "Not WGD")) %>%
  group_by(case_id) %>%
  mutate(any_wgd = sum(high_wgd_clone == "WGD" | case_id %in% wgd_cases | id %in% hyper_rescues)) %>%
  mutate(any_group = ifelse(sum(group_rescued == "Near-Haploid") > 0, "Near-Haploid", ifelse(sum(group_rescued == "Low-Hypodiploid") > 0, "Low-Hypodiploid", "Other"))) %>%
  distinct(case_id, .keep_all = T) %>%
  mutate(wgd = ifelse(any_wgd > 0, "WGD", "No WGD")) %>%
  mutate(any_group = ifelse(any_group == "Near-Haploid", "NH", ifelse(any_group == "Low-Hypodiploid", "LH", any_group))) %>%
  ggplot(aes(x = factor(any_group, levels = c("Other", "LH", "NH")), fill = wgd)) +
  geom_bar(position = "fill") +
  labs(fill = "", x = "", y = "Proportion of Cases with WGD") +
  theme_large_classic() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(`No WGD` = "black", `WGD` = "lightsteelblue")) +
  ggtitle("(a)") 

tcga_wgd_rate <- tcga_classes %>%
  mutate(group = ifelse(group == "Near-Haploid", "NH", ifelse(group == "Low-Hypodiploid", "LH", group))) %>%
  ggplot(aes(x = factor(group, levels = c("Other", "LH", "NH")), fill = wgd)) +
  geom_bar(position = "fill") +
  ggtitle("TCGA") +
  labs(fill = "", x = "", y = "Proportion of Cases with WGD") +
  scale_fill_manual(values = c(`No WGD` = "black", `WGD` = "lightsteelblue")) +
  theme_large_classic() +
  theme(legend.position = "bottom")

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
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  ggtitle('(e)') +
  labs(x = "", y = "Copy Number Heterogeneity") +
  scale_fill_manual(values = short_palette)

cnh_violins_nondoubledhypos <- tcga_cnh %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  mutate(Class = factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ggplot(aes(x = Class, y = CNH)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  labs(x = "", y = "Copy Number Heterogeneity") +
  ggtitle('(f)') +
  scale_fill_manual(values = class_palette)


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
  stat_compare_means(hjust = -1) +
  scale_fill_manual(values = c(`FALSE` = "#D9D9D9", `TRUE` = "#8DA0CB")) +
  theme_large() +
  theme(legend.position = "none") +
  xlab("Stereotyped Hypodiploid") +
  ggtitle('(g)')

## number of copy number segments 
segs_summary <- fasc %>%
  mutate(len = (as.numeric(End) + 1) - as.numeric(Start)) %>%
  group_by(proj, GDC_Aliquot, Class, wgd, group) %>%
  summarize(num_segments = n()) %>%
  mutate(Class = factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ungroup()

segs_violins <- segs_summary %>%
  ggplot(aes(x = Class, y = num_segments)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  scale_y_continuous(trans = "log10") +
  labs(x = "", y = "Number of copy number segments") + 
  scale_fill_manual(values = class_palette) +
  ggtitle('(c)')

segs_violins_nondoubledhypos <- segs_summary %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  ggplot(aes(x = Class, y = num_segments)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  scale_y_continuous(trans = "log10") +
  labs(x = "", y = "Number of copy number segments") + 
  ggtitle('(c)') +
  scale_fill_manual(values = class_palette)

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
  stat_compare_means() +
  scale_fill_manual(values = c(`FALSE` = "#D9D9D9", `TRUE` = "#8DA0CB")) +
  theme_large() +
  theme(legend.position = "none") +
  ggtitle('(g)') +
  xlab("Stereotyped Hypodiploid") + ylab('Number of CN segments')

## max chromosome proportion
chr_prop_summary <- fasc %>%
  mutate(len = (as.numeric(End) + 1) - as.numeric(Start)) %>%
  group_by(proj, GDC_Aliquot, Chromosome, Class, wgd, group) %>%
  summarize(prop_max = max(len) / sum(len)) %>%
  mutate(Class = factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ungroup() 
  
maxprop_violins_func <- function(df, sub) {
  df %>%
    ggplot(aes(x = Class, y = prop_max)) +
    geom_violin(aes(fill = Class)) +
    geom_boxplot(alpha = 0.2) +
    stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
    theme_large(base_size = 14) +
    theme(legend.position = "none") +
    labs(x = '', y = "Proportion of chromosome in longest CN segment") +
    ggtitle(sub) +
    scale_fill_manual(values = class_palette)
  
}
maxprop_violins <- chr_prop_summary %>% maxprop_violins_func(sub = '(b)')
  
maxprop_violins_nodoubledhypos <- chr_prop_summary %>%
  filter(group == "Other" | wgd == "No WGD") %>%
 maxprop_violins_func(sub = '(d)')


# Chromothripsis bar 
chromothripsis <- tcga_classes %>%
  left_join(thrip, by = c("Patient" = "tcga_id")) %>%
  filter(!is.na(is_chromoanagenesis)) %>%
  mutate(chromothripsis = is_chromoanagenesis == 1) %>%
  mutate(Class = factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  ggplot(aes(x = Class, fill = chromothripsis)) +
  geom_bar(position = "fill") +
  labs(fill = "Chromothripsis", x = "", y = "Proportion of Cases") +
  theme_large_classic() +
  scale_fill_manual(values = c(`FALSE` = "black", `TRUE` = "lightsteelblue")) +
  theme(legend.position = "bottom") +
  ggtitle('(e)')

## Survival Analysis 
forsurv <- tcga_classes %>%
  left_join(clin, by = c('Patient' = 'submitter_id')) %>%
  mutate(last = case_when(vital_status == 'Dead' ~ days_to_death, vital_status == 'Alive' ~    days_to_last_follow_up)) %>%
  filter(!is.na(last)) %>%
  mutate(code = case_when(vital_status == 'Dead' ~ 2, vital_status == 'Alive' ~ 1)) %>%
  mutate(Class = factor(Class, levels = c('Diploid', 'Near-Haploid', 'Low-Hypodiploid', 'Aneuploid', 'Polyploid'))) 


hypo_effect <- survfit(Surv(last, code) ~ Class, data = forsurv %>% mutate(Class = as.character(Class)) %>% mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>% mutate(Class = factor(Class, levels = c('Diploid', 'NH', 'LH', 'Aneuploid', 'Polyploid'))))
labels_wrapped <- sub("Class=", "", names(hypo_effect$strata))
surv_ploidy_classes <- ggsurvplot(hypo_effect, data = forsurv, ggtheme = theme_large() + theme(legend.text = element_text(size = 16)), legend.title = '', legend.labs = labels_wrapped)

wgd_effect_hypos <- survfit(Surv(last, code) ~ wgd, data = forsurv %>% filter(Class == 'Low-Hypodiploid'))
surv_wgd_hypo <- ggsurvplot(wgd_effect_hypos, data = forsurv %>% filter(Class == 'Low-Hypodiploid'), ggtheme = theme_large() + theme(legend.text = element_text(size = 16)), pval = T, palette = c("#E41A1C", "#377EB8"), conf.int = TRUE, legend.title = '', legend.labs = c('No WGD', 'WGD'))

# surv: stable vs CIN

forsurv_filtered <- forsurv %>%
  mutate(Class = as.character(Class)) %>%
  filter(Class %in% c("Low-Hypodiploid"), proj %in% enough_lh) %>%
  mutate(typ = ifelse(proj %in% c("KICH", "ACC"), "Stable", "CIN")) 

km_fit <- survfit(Surv(last, code) ~ typ, data = forsurv_filtered)

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

## MH score sensitivity in TCGA

# < 76 autosomes 
mh_tcga_76 <- tcga_scores %>% left_join(tcga_classes, by = c('GDC_Aliquot')) %>% filter(group != 'Other', wgd == 'WGD', n_chr_nosex <= 76) %>% group_by(proj) %>% summarize(n_dh = n(), sens = mean(diff > 0)) %>% filter(n_dh >= 15) %>% ggplot(aes(x = reorder(proj, sens), y = sens)) + geom_col(fill = 'black', colour = 'black') + theme_large() + labs(x = '', y = 'MH score > 0', subtitle = 'Up to 76 autosomes') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + ggtitle('(h)') 
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
ggsave(plot = fig3_hypo_paper, file = 'paper/hypo_fig4.pdf', width = 25, height = 30)

supp_fig3 <- (mit_wgd_rate_hyper_rescues | maxprop_violins) / (segs_violins_nondoubledhypos | maxprop_violins_nodoubledhypos) / (chromothripsis | cnh_violins_nondoubledhypos) / (stereotyped_segs | (mh_tcga_76))
ggsave(plot = supp_fig3, file = 'paper/supp_fig3_CIN.png', width = 25, height = 30)
ggsave(plot = supp_fig3, file = 'paper/hypo_supp_fig3_CIN.pdf', width = 25, height = 30)

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
coxph(Surv(last, code) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid'))) %>% tidy()
# LH v Diploid - controlling for cancer type - p = 0.0237
coxph(Surv(last, code) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid'))) %>% tidy()
# NH v Diploid - p = 0.0251
coxph(Surv(last, code) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Near-Haploid', 'Diploid'))) %>% tidy()
# LH v Aneuploid - p = 0.00819
coxph(Surv(last, code) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Aneuploid'))) %>% tidy()
# LH v Aneuploid - controlling for cancer type - p = 0.0179
coxph(Surv(last, code) ~ Class + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Aneuploid'))) %>% tidy()
# n.s. with polyploid
coxph(Surv(last, code) ~ Class + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class %in% c('Low-Hypodiploid', 'Polyploid'))) %>% tidy()
# doubled v non-doubled LH - n.s. (0.192)
coxph(Surv(last, code) ~ wgd + age_at_index + gender + race, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class == 'Low-Hypodiploid')) %>% tidy()
# doubled v non-doubled LH by proj - n.s. (0.588)
coxph(Surv(last, code) ~ wgd + age_at_index + gender + race + proj, data = forsurv %>% mutate(Class = as.character(Class)) %>% filter(Class == 'Low-Hypodiploid')) %>% tidy()
# stable v CIN hypodiploids
coxph(Surv(last, code) ~ typ + age_at_index + gender + race, data = forsurv_filtered) %>% tidy()

