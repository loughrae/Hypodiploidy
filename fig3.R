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

enough_lh <- tcga_classes %>% filter(group == 'Low-Hypodiploid') %>% count(proj) %>% filter(n >= 15) %>% pull(proj) 


## summarise WGD and hypo rates?


## fig3a: WGD rate by hypodiploidy status
mit_wgd_rate <- mitcn_meta %>%
  mutate(high_wgd_clone = ifelse(ploidy > 2.7, "WGD", "Not WGD")) %>%
  group_by(case_id) %>%
  mutate(any_wgd = sum(high_wgd_clone == "WGD" | case_id %in% wgd_cases | case_id %in% rescues)) %>%
  mutate(any_group = ifelse(sum(group == "Near-Haploid") > 0, "Near-Haploid", ifelse(sum(group == "Low-Hypodiploid") > 0, "Low-Hypodiploid", "Other"))) %>%
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
  labs(subtitle = "(b)", x = "Low-Hypodiploidy Rate", y = "WGD Rate (excluding LH)") 

tcga_cnh <- tcga_classes %>%
  left_join(cnh, by = c("Patient")) %>%
  filter(!is.na(CNH)) 

## fig3e: CNH by ploidy class (violins)
cnh_violins <- tcga_cnh %>%
  mutate(Class = factor(Class, levels = c("Near-Haploid", "Low-Hypodiploid", "Near-Diploid", "Other", "WGD-high"))) %>%
  ggplot(aes(x = Class, y = CNH)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  labs(subtitle = "(e)", x = "", y = "Copy Number Heterogeneity") +
  scale_fill_manual(values = class_palette)

## fig3f: CNH vs hypo rate
hypo_cnh_CI <- tcga_cnh %>%
  filter(!proj %in% c("KICH", "ACC")) %>%
  group_by(proj) %>%
  summarise(wgd_rate = sum(wgd == "WGD") / n(), hypo_rate = mean(group == "Low-Hypodiploid"), median_cnh = median(CNH[Class %in% c("Near-Diploid", "Other")])) %>%
  ggplot(aes(x = hypo_rate, y = median_cnh)) +
  geom_point() +
  geom_smooth(size = 1, method = "lm", colour = "black", fill = "lightsteelblue") +
  stat_cor() +
  geom_label_repel(aes(label = proj)) +
  xlab("Low-Hypodiploidy Rate") +
  ylab("Median copy number heterogeneity") +
  theme_large(base_size = 18) +
  labs(subtitle = "(d)")

segs_violins <- fasc %>%
  mutate(len = (as.numeric(End) + 1) - as.numeric(Start)) %>%
  group_by(proj, GDC_Aliquot, Class, wgd, group) %>%
  summarize(num_segments = n(), prop_max = max(len) / sum(len)) %>%
  mutate(Class = factor(Class, levels = c('Near-Haploid', 'Low-Hypodiploid', 'Near-Diploid', 'Other', 'WGD-high'))) %>%
  ggplot(aes(x = Class, y = num_segments)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
  theme_large(base_size = 18) +
  theme(legend.position = "none") +
  scale_y_continuous(trans = "log10") +
  labs(subtitle = "(????)", x = "", y = "Number of copy number segments") + 
  scale_fill_manual(values = class_palette)

hypo_segs_CI <- fasc %>%
  mutate(len = (as.numeric(End) + 1) - as.numeric(Start)) %>%
  group_by(proj, GDC_Aliquot, Class, wgd, group) %>%
  summarize(num_segments = n()) %>%
  filter(!proj %in% c("KICH", "ACC")) %>%
  group_by(proj) %>%
  summarize(wgd_rate = sum(wgd == "WGD") / n(), hypo_rate = mean(group == "Low-Hypodiploid"), median_nseg_sample_normalploidy = median(num_segments[Class %in% c("Near-Diploid", "Other")])) %>%
  ggplot(aes(x = hypo_rate, y = median_nseg_sample_normalploidy)) +
  geom_point() +
  geom_smooth(size = 1, method = "lm", colour = "black", fill = 'lightsteelblue') +
  stat_cor() +
  geom_label_repel(aes(label = proj)) +
  xlab("Hypodiploidy Rate") +
  ylab("Median Segment Count") +
  theme_large(base_size = 18) +
  labs(subtitle = "(e)")

stereotyped <- tcga_cnh %>%
  filter(Class == "Low-Hypodiploid", proj %in% enough_lh) %>%
  mutate(highlight = proj %in% c("ACC", "KICH")) %>%
  ggplot(aes(x = highlight, y = CNH, fill = highlight)) +
  geom_violin() +
  geom_boxplot(alpha = 0.2) +
  geom_jitter(colour = "black", alpha = 0.2, height = 0) +
  stat_compare_means() +
  scale_fill_manual(values = c(`FALSE` = "#D9D9D9", `TRUE` = "#8DA0CB")) +
  theme_large() +
  theme(legend.position = "none") +
  xlab("Stereotyped Hypodiploid")


## Survival Analysis 
forsurv <- tcga_classes %>%
  left_join(clin, by = c('Patient' = 'submitter_id')) %>%
  mutate(last = case_when(vital_status == 'Dead' ~ days_to_death, vital_status == 'Alive' ~    days_to_last_follow_up)) %>%
  filter(!is.na(last)) %>%
  mutate(code = case_when(vital_status == 'Dead' ~ 2, vital_status == 'Alive' ~ 1)) %>%
  mutate(Class = factor(Class, levels = c('Near-Diploid', 'Near-Haploid', 'Low-Hypodiploid', 'Other', 'WGD-high'))) 


hypo_effect <- survfit(Surv(last, code) ~ Class, data = forsurv)
surv_ploidy_classes <- ggsurvplot(hypo_effect, data = forsurv)

## coxPH


fig3_hypo_paper <- (((mit_wgd_rate | tcga_wgd_rate) | wgd_v_hypo) + plot_layout(widths = c(1,1,2))) / (segs_violins | hypo_segs_CI) / (cnh_violins | hypo_cnh_CI) / (stereotyped | as.ggplot(surv_ploidy_classes$plot))

ggsave(plot = fig3_hypo_paper, file = 'paper/fig3_hypo_paper.png', width = 25, height = 30)