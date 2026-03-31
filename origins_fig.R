source('setup.R')
tcga_classes <- fread('TCGA_ploidy_classes.tsv')

print('Making Figure 4')

## load mutations 
mafs <- fread('combined_mafs.txt') %>%
  mutate(Patient = paste(TCGA, Centre, Indiv, sep = '-')) %>% 
  group_by(Patient) %>% mutate(dups = n_distinct(Tumor_Sample_UUID)) %>% filter(dups == 1) %>% select(-dups) %>%
  left_join(tcga_classes, by = c('Patient')) %>% 
  filter(!is.na(wgd)) %>%
  ungroup()

## count total mutations per sample

enough_lh <- tcga_classes %>% filter(group == 'Low-Hypodiploid') %>% count(proj) %>% filter(n >= 15) %>% pull(proj) 
enough_dip <- tcga_classes %>% filter(Class == 'Diploid') %>% count(proj) %>% filter(n >= 15) %>% pull(proj) 

tcga_mu <- mafs %>%
  group_by(Patient, ploidy, Class, proj) %>%
  summarize(total_mu = n(), total_non_syn_mu = sum(Variant_Classification != "Silent")) %>%
  ungroup() 

purity <- fread("~/Downloads/Aran2015_purity.csv", skip = 3) %>%
  clean_names() %>%
  mutate(Patient = str_extract(sample_id, "^[^-]+-[^-]+-[^-]+"))

tcga_purity <- tcga_classes %>% left_join(purity, by = c('Patient')) %>% add_count(Patient, name = 'purities_per_patient') %>% filter(purities_per_patient == 1) %>% filter(!is.na(cpe))

mu_perploidy_reg <- tcga_mu %>%
  filter(Class %in% c("Diploid", "Low-Hypodiploid")) %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  ungroup() %>%
  do(tidy(lm(mu_per_ploidy ~ Class + proj, data = .))) # near diploids have higher mutations_per_ploidy after controlling for cancer type but that's actually a MEAN bc of hypermutation

#controlling for purity
tcga_mu %>%
  left_join(tcga_purity[, c("Patient", "cpe")], by = c("Patient")) %>%
  filter(!is.na(cpe)) %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  mutate(Class = factor(Class, levels = c("Low-Hypodiploid", "Near-Haploid", "Diploid", "Aneuploid", "Polyploid"))) %>%
  do(tidy(lm(mu_per_ploidy ~ Class + proj + cpe, data = .)))

total_mus_by_class <- tcga_mu %>%
  left_join(tcga_classes, by = c("Class", "ploidy", "proj", "Patient")) %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  { 
    counts <- count(., Class)
  ggplot(., aes(x = Class, y = total_mu, fill = Class)) +
  geom_violin() +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Diploid", "Low-Hypodiploid"), c("Diploid", "Near-Haploid"))) +
  theme_large() +
  xlab("") +
  labs(subtitle = '(a)') +
  scale_y_continuous(trans = "log10") +
  ylab("Total Mutations (log10)") +
  scale_fill_manual(values = class_palette) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = function(x) {
      n <- counts$n[match(x, counts$Class)]
      paste0(x, "\n(n = ", n, ")")
    }
  )
  }

total_mutations_pval_table <- tcga_mu %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  mutate(log10_total_mu = log10(total_mu), log10_mupp = log10(mu_per_ploidy)) %>%
  rstatix::pairwise_wilcox_test(
    log10_total_mu ~ Class,
    comparisons = list(c("Diploid", "Low-Hypodiploid"), c("Diploid", "Near-Haploid")),
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


total_mutations_pval_table_plot <- ggtexttable(
  total_mutations_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(total_mutations_pval_table_plot, file = 'paper/tables/total_mutations_pval_table.pdf', width = 4, height = 2, units = 'in')


mu_per_ploidy_by_class <- tcga_mu %>%
  left_join(tcga_classes, by = c("Class", "ploidy", "proj", "Patient")) %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  { 
    counts <- count(., Class)
    ggplot(., aes(x = Class, y = mu_per_ploidy, fill = Class)) +
      geom_violin() +
      geom_boxplot() +
      stat_compare_means(comparisons = list(c("Diploid", "Low-Hypodiploid"), c("Diploid", "Near-Haploid"))) +
      theme_large() +
      xlab("") +
      labs(subtitle = '(b)') +
      scale_y_continuous(trans = "log10") +
      ylab("Ploidy-Corrected Mutation Count (log10)") +
      scale_fill_manual(values = class_palette) +
      theme(legend.position = "none") +
      scale_x_discrete(labels = function(x) {
        n <- counts$n[match(x, counts$Class)]
        paste0(x, "\n(n = ", n, ")")
      }
      )
  }

mu_pp_pval_table <- tcga_mu %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  mutate(log10_total_mu = log10(total_mu), log10_mupp = log10(mu_per_ploidy)) %>%
  rstatix::pairwise_wilcox_test(
    log10_mupp ~ Class,
    comparisons = list(c("Diploid", "Low-Hypodiploid"), c("Diploid", "Near-Haploid")),
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


mu_pp_pval_table_plot <- ggtexttable(
  mu_pp_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(mu_pp_pval_table_plot, file = 'paper/tables/mu_pp_pval_table.pdf', width = 4, height = 2, units = 'in')



mu_per_ploidy_by_tissue <- tcga_mu %>%
  left_join(tcga_classes, by = c("Class", "ploidy", "proj", "Patient")) %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  filter(Class %in% c("Low-Hypodiploid", "Diploid")) %>%
  filter(proj %in% intersect(enough_lh, enough_dip)) %>%
  ggplot(aes(x = Class, y = mu_per_ploidy)) +
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(~proj) +
  scale_y_continuous(trans = "log10") +
  theme_large()

## regression to detect mutations depleted/enriched in LH and poly vs diploid (only for COSMIC genes which are in mafs)

cosmic <- read.table('~/Downloads/Census_allWed Oct 25 16 41 44 2023.tsv', sep = '\t', header = T) %>% clean_names() #may need to require presence of gene in mafs

mafs_pivot <- mafs %>%
  filter(Hugo_Symbol %in% cosmic$gene_symbol) %>% #filter to COSMIC genes
  filter(Variant_Classification != "Silent") %>% #non-synonymous
  distinct(Patient, Hugo_Symbol) %>%
  mutate(dummy = 1) %>% 
  pivot_wider(names_from = Hugo_Symbol, values_from = dummy, values_fill = list(dummy = 0)) 

status <- tcga_mu %>%
  left_join(mafs_pivot, by = c("Patient")) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>% #patients with no mutations in COSMIC genes
  pivot_longer(-c(Patient, total_mu, total_non_syn_mu, proj, Class, ploidy), names_to = "gene", values_to = "mutated", values_drop_na = F) %>% #NB: the mutate across above is necessary for this to work
  left_join(tcga_purity[, c('Patient', 'cpe')], by = c('Patient')) %>%
  filter(!is.na(cpe))

mu_lh <- status %>%
  filter(Class %in% c("Diploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c('Diploid', 'Low-Hypodiploid'))) %>%
  group_by(gene) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 5) %>%
  do(tidy(glm(mutated ~ Class + total_non_syn_mu + cpe + proj, family = 'binomial', data = .))) %>% ungroup() # 

mu_poly <- status %>%
  filter(Class %in% c("Diploid", "Polyploid")) %>%
  mutate(Class = factor(Class, levels = c('Diploid', 'Polyploid'))) %>%
  group_by(gene) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 5) %>%
  do(tidy(glm(mutated ~ Class + total_non_syn_mu + cpe + proj, family = 'binomial', data = .))) %>% ungroup() # 

mu_poly %>% filter(term == 'ClassPolyploid') %>% mutate(bh = p.adjust(p.value, method = 'BH')) %>% filter(bh < 0.05) %>% filter(gene %in% c('ACVR2A', 'BRAF', 'CTCF', 'TP53', 'RPL22', 'ARID1A'))


mu_hap_v_lh <- status %>%
  filter(Class %in% c("Near-Haploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c('Low-Hypodiploid', 'Near-Haploid'))) %>%
  group_by(gene) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 3) %>% #lower here note
  do(tidy(glm(mutated ~ Class + total_non_syn_mu + cpe + proj, family = 'binomial', data = .))) %>% 
  ungroup()

mu_hap_v_lh %>% filter(term == 'ClassNear-Haploid') %>% mutate(bh = p.adjust(p.value, method = 'BH')) %>% filter(bh < 0.05) #0 genes
mu_hap_v_lh %>% filter(term == 'ClassNear-Haploid') %>% mutate(bh = p.adjust(p.value, method = 'BH')) %>% filter(p.value < 0.05) #10 genes

mu_lh_v_aneuploid <- status %>%
  filter(Class %in% c("Aneuploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c('Aneuploid', 'Low-Hypodiploid'))) %>%
  group_by(gene) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 5) %>% 
  do(tidy(glm(mutated ~ Class + total_non_syn_mu + cpe + proj, family = 'binomial', data = .))) %>% 
  ungroup()

mu_lh_v_aneuploid %>% filter(term == 'ClassLow-Hypodiploid') %>% mutate(bh = p.adjust(p.value, method = 'BH')) %>% filter(bh < 0.05) #ARID1A, CTCF and RPL22 are nominally diff though (depleted in LH)


## fig4a: volcano plot of mutations in LH
lh_volcano <- mu_lh %>%
  filter(term == "ClassLow-Hypodiploid") %>%
  mutate(bh = p.adjust(p.value, method = "BH")) %>%
  mutate(colour_ind = bh < 0.05) %>%
  filter(estimate > -10) %>% #remove outliers
  mutate(logp = -1*log10(p.value)) %>%
  mutate(gene_label = ifelse(colour_ind == F, "", paste0(gene, '\n', round(logp, 1)))) %>%
  mutate(gene_label = ifelse(logp > 10, paste0(gene_label, '*'), gene_label)) %>%
  mutate(logp = pmin(logp, 10)) %>% #cap TP53
  ggplot(aes(colour = colour_ind, x = estimate, y = logp)) +
  geom_point() +
  geom_label_repel(aes(label = gene_label)) +
  scale_colour_manual(values = c(`FALSE` = "darkgreen", `TRUE` = "darkblue")) +
  theme_large() +
  labs(colour = "FDR < 0.05", subtitle = '(a)') +
  labs(y = expression(-log[10](p.value)), x = "Log Odds Ratio") 

aneu_volcano <- mu_lh_v_aneuploid %>%
  filter(term == "ClassLow-Hypodiploid") %>%
  mutate(bh = p.adjust(p.value, method = "BH")) %>%
  mutate(colour_ind = bh < 0.05) %>%
  filter(estimate > -10) %>% #remove outliers
  mutate(logp = -1*log10(p.value)) %>%
  mutate(gene_label = ifelse(colour_ind == F, "", paste0(gene, '\n', round(logp, 1)))) %>%
  mutate(gene_label = ifelse(logp > 10, paste0(gene_label, '*'), gene_label)) %>%
  mutate(logp = pmin(logp, 10)) %>% #cap TP53
  ggplot(aes(colour = colour_ind, x = estimate, y = logp)) +
  geom_point() +
  geom_label_repel(aes(label = gene_label)) +
  scale_colour_manual(values = c(`FALSE` = "darkgreen", `TRUE` = "darkblue")) +
  theme_large() +
  labs(colour = "FDR < 0.05", subtitle = '(c)') +
  labs(y = expression(-log[10](p.value)), x = "Log Odds Ratio") 


## fig4b: mutations in LH vs Poly
lh_v_poly_mus <- mu_lh %>%
  filter(term == "ClassLow-Hypodiploid") %>%
  left_join(mu_poly, by = c("gene"), suffix = c("_lh", "_poly")) %>%
  filter(term_poly == "ClassPolyploid") %>%
  mutate(bh_lh = p.adjust(p.value_lh, method = "BH"), bh_poly = p.adjust(p.value_poly, method = "BH")) %>%
  filter(estimate_lh > -10) %>%
  mutate(colour_ind = ifelse(bh_lh < 0.05 & bh_poly < 0.05, "Both", ifelse(bh_lh < 0.05, "Hypodiploids", ifelse(bh_poly < 0.05, "Polyploids", ifelse(bh_lh >= 0.05 & bh_poly >= 0.05, "n.s.", "what"))))) %>%
  mutate(colour_ind = factor(colour_ind, levels = c('Both', 'Hypodiploids', 'Polyploids', 'n.s.'))) %>%
  mutate(gene_label = ifelse(colour_ind == "n.s.", "", gene)) %>%
  ggplot(aes(x = estimate_lh, y = estimate_poly)) +
  geom_point(aes(colour = colour_ind)) +
  geom_smooth(method = "lm", colour = "black", fill = "lightsteelblue") +
  stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top",
           p.accuracy = 1e-100) +  
  theme_large() +
  labs(colour = '', subtitle = '(c)') +
  geom_label_repel(aes(label = gene_label, colour = colour_ind)) +
  xlab("Log Odds Ratio in Hypodiploids") +
  ylab("Log Odds Ratio in Polyploids") +
  scale_colour_manual(values = c(`Both` = "darkgreen", `Polyploids` = "darkblue", `n.s.` = "lightgrey", `Hypodiploids` = 'red'))

## fig4c: plot TP53 non-synonymous mutation rate
p53_barchart <- mafs %>%
  group_by(wgd, Class, group, Patient, Sample) %>%
  summarize(p53 = sum(Hugo_Symbol == "TP53" & Variant_Classification != 'Silent')) %>%
  group_by(Class) %>%
  mutate(perc_p53 = sum(p53 > 0) / n()) %>%
  mutate(Class = ifelse(Class == "Low-Hypodiploid", "LH", ifelse(Class == "Near-Haploid", "NH", Class))) %>%
  mutate(Class = factor(Class, levels = c("NH", "LH", "Diploid", "Aneuploid", "Polyploid"))) %>%
  mutate(TP53 = ifelse(p53 > 0, 'Mut', 'WT')) %>%
  mutate(TP53 = factor(TP53, levels = c('WT', 'Mut'))) %>%
  ggplot(aes(x = Class, fill = TP53)) +
  geom_bar(position = "fill") +
  xlab("") +
  labs(subtitle = '(b)') +
  ylab("Proportion of Cases") +
  theme_large_classic() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c(`WT` = "black", `Mut` = "lightsteelblue")) +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


## fig4d: Microsatellite instability
msi <- fread('~/Downloads/mantis_msi.csv') %>% clean_names()

msi_proc <- tcga_classes %>% 
  left_join(msi, by = c("Patient" = "case_id")) %>%
  filter(!is.na(mantis_score)) %>%
  filter(proj %in% c("UCEC", "COAD", "STAD", "ACC")) %>%
  mutate(proj = factor(proj, levels = c("COAD", "UCEC", "STAD", "ACC"))) %>%
  filter(Class %in% c("Diploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c("Diploid", "Low-Hypodiploid"))) %>%
  mutate(Class = ifelse(Class == "Diploid", "Dip", "LH")) %>%
  group_by(proj)

msi_bar <- msi_proc %>%
  ggplot(aes(x = Class, fill = mantis_score > 0.4)) +
  geom_bar(position = "fill") +
  facet_wrap(~proj, nrow = 1) +
  theme_large() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c(`FALSE` = "black", `TRUE` = "lightsteelblue")) +
  ylab("Proportion of Cases") +
  xlab("") +
  labs(fill = "MANTIS > 0.4", subtitle = '(g)') +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

msi_proc %>%
  left_join(tcga_purity[, c("Patient", "cpe")], by = c("Patient")) %>%
  filter(!is.na(cpe)) %>%
  mutate(high = as.numeric(mantis_score > 0.4)) %>%
  group_by(proj) %>%
  do(tidy(glm(high ~ Class + cpe, family = "binomial", data = .)))

msi_proc %>%
  left_join(tcga_purity[, c("Patient", "cpe")], by = c("Patient")) %>%
  filter(!is.na(cpe)) %>%
  mutate(high = as.numeric(mantis_score > 0.4)) %>%
  group_by(proj) %>%
  do(tidy(glm(mantis_score ~ Class + cpe, data = .)))

msi_tests <- msi_proc %>%
  mutate(msi_h = mantis_score > 0.4) %>%
  group_by(proj, Class) %>%
  summarise(n = n(), high = sum(msi_h), low = n() - high, prop = high / n) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = Class, values_from = c(n, high, low, prop)) %>%
  rowwise() %>%
  mutate(
    test = list(fisher.test(matrix(c(high_Dip, low_Dip, high_LH, low_LH), nrow = 2))),
    p_value = test$p.value,
    ratio = prop_LH / prop_Dip
  ) %>%
  select(proj, prop_Dip, prop_LH, ratio, p_value)
  
msi_violin <- msi_proc %>%
  ggplot(aes(x = Class, y = mantis_score, colour = Class)) +
  geom_jitter(height = 0) +
  geom_boxplot(alpha = 0.2) +
  facet_wrap(~proj, nrow = 1) +
  theme_large() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c(`Dip` = "darkgreen", `LH` = "darkblue")) +
  ylab("MANTIS score") +
  xlab("") +
  labs(subtitle = '(d)') +
  theme(legend.position = "none") +
  stat_compare_means() +
  geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black") +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


## Hypoxia 

hyp <- fread('~/Downloads/41588_2018_318_MOESM3_ESM.txt') %>% mutate(Patient = gsub('.', '-', patient_id, fixed = T)) %>% #7791 pts, all distinct
  left_join(tcga_classes, by = c('Patient')) %>% 
  filter(!is.na(Class))

# rna_uuids <- UUIDtoUUID(unique(rag_cpm$V3))
# rna_barcodes <- UUIDtoBarcode(rna_uuids$cases.case_id)
# rna_uuids %>% left_join(rna_barcodes, by = c('cases.case_id' = 'case_id')) %>% add_count(submitter_id) %>% filter(n == 1) %>% write.table('rna_ids_singletons.tsv', sep = '\t', quote = F, row.names = F, col.names = T)

rna_ids <- fread('rna_ids_singletons.tsv')

rag_cpm <- fread("ragnum_cpms_filtered.tsv") %>%
  left_join(rna_ids, by = c("V3" = "file_id")) %>%
  filter(!is.na(submitter_id)) %>% #no multiple RNA samples
  filter(submitter_id %in% tcga_classes$Patient) %>%
  left_join(tcga_classes, by = c('submitter_id' = 'Patient'))

cancers <- fread('~/Downloads/bhandari_cancers.tsv', header = F)
#CPM scores and Bhandari Ragnum scores correlate at r = 0.97 (in the 20 cancers they studied). That'll do. Definitely a bunch of outliers tho which is pretty fuckin weird.
rag_cpm %>% filter(proj %in% cancers$V1) %>% group_by(V2) %>% mutate(gene_med = median(V4)) %>% mutate(add = ifelse(V4 > gene_med, 1, -1)) %>% select(V2, V4, submitter_id, proj, gene_med, add) %>% group_by(submitter_id, proj) %>% summarize(ragnum_score_pan = sum(add)) %>% left_join(hyp, by = c('submitter_id' = 'Patient')) %>% filter(!is.na(Ragnum_hypoxia_score_pan_cancer)) %>% ggplot(aes(x = ragnum_score_pan, y = Ragnum_hypoxia_score_pan_cancer)) + geom_point() + geom_smooth(method = 'lm') + stat_cor()


rag_scores <- rag_cpm %>% #all cancers now
  group_by(V2) %>%
  mutate(gene_med = median(V4)) %>%
  mutate(add = ifelse(V4 > gene_med, 1, -1)) %>%
  select(V2, V4, submitter_id, proj, gene_med, add, Class, wgd, group, ploidy) %>%
  group_by(submitter_id, proj, Class, wgd, group, ploidy) %>%
  summarize(ragnum_score_pan = sum(add)) %>%
  ungroup()


rag_scores %>% ggplot(aes(x = reorder(proj, ragnum_score_pan, FUN = median), y = ragnum_score_pan)) + geom_boxplot()

ploidy_rates <- fread('tissue_ploidy_rates.tsv')

hypo_v_ragnum <- rag_scores %>%
  group_by(proj) %>%
  filter(!proj %in% c("ACC", "KICH")) %>%
  summarize(wgd_rate = sum(wgd == "WGD") / n(), hypo_rate = mean(group == "Low-Hypodiploid"), avg_hyp_nohypos = median(ragnum_score_pan[group != "Low-Hypodiploid"]), avg_hyp_nowgd = median(ragnum_score_pan[wgd == "No WGD"])) %>%
  ggplot(aes(x = avg_hyp_nohypos, y = hypo_rate)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", colour = "black", fill = "lightsteelblue") +
  stat_cor() +
  xlab("Median hypoxia score among non-LH tumours") +
  ylab("Hypodiploidy Rate") +
  labs(subtitle = '(g)') +
  theme_large() +
  geom_label_repel(aes(label = proj))

bts <- 12

pairwise_comparisons <- list(
  c("Diploid", "Polyploid"),
  c("Diploid", "Low-Hypodiploid"),
  c("Polyploid", "Low-Hypodiploid"),
  c("Aneuploid", "Low-Hypodiploid"),
  c("Near-Haploid", "Low-Hypodiploid")
)

hyp_violins <- rag_scores %>% 
  mutate(Class = factor(Class, levels = c('Near-Haploid', 'Low-Hypodiploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  ggplot(aes(x = Class, y = ragnum_score_pan)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  theme_large() +
  xlab("") +
  ylab("Ragnum Hypoxia Score") +
  theme(legend.position = "none", axis.text.x = element_text(size = bts), axis.title.y = element_text(size = bts + 2)) +
  ggtitle("(e)") +
  scale_fill_manual(values = class_palette) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

#compute exact p values for supplementary table
hypoxia_pval_table <- rag_scores %>%
  rstatix::pairwise_wilcox_test(
    ragnum_score_pan ~ Class,
    comparisons = pairwise_comparisons,
    p.adjust.method = "none"
  ) %>%
  select(group1, group2, p) %>%
  rename(Comparison = group1, P_Value = p) %>%
  mutate(Comparison = paste(Comparison, "vs", group2)) %>%
  select(-group2)


hypoxia_pval_table_plot <- ggtexttable(
  hypoxia_pval_table,
  rows = NULL,             
  theme = ttheme("classic") 
)

ggsave(hypoxia_pval_table_plot, file = 'paper/tables/hypoxia_pval_table.pdf', width = 4, height = 2, units = 'in')


hyp_ploidy <- rag_scores %>% ggplot(aes(x = ploidy, y = ragnum_score_pan)) +
  geom_point(alpha = 0.2) +
  geom_smooth(colour = 'black', fill = 'lightsteelblue') +
  theme_large() +
  ylab("Ragnum Hypoxia Score") +
  xlab("Tumour Ploidy") +
  theme(axis.text.x = element_text(size = bts), axis.title.y = element_text(size = bts + 2), axis.title.x = element_text(size = bts + 2)) +
  ggtitle("(f)") 

stable_hyp <- rag_scores %>%
  filter(proj %in% enough_lh) %>%
  mutate(proj2 = ifelse(proj %in% c("ACC", "KICH"), "Stereotyped", "CIN")) %>%
  filter(Class == "Low-Hypodiploid") %>%
  ggplot(aes(x = reorder(proj2, ragnum_score_pan, FUN = median), y = ragnum_score_pan, fill = proj2)) +
  geom_violin() +
  geom_boxplot(alpha = 0.2) +
  stat_compare_means(method = 'wilcox.test', aes(label = after_stat(paste0('p = ', format(p, scientific = TRUE, digits = 4))))) + 
  theme_large() +
  theme(legend.position = "none") +
  xlab("") +
  labs(subtitle = '(h)') +
  ylab("Ragnum hypoxia score") +
  scale_fill_manual(values = c(`Stereotyped` = "#3288bd", `CIN` = "#1b9e77")) +
  geom_text(data = . %>% count(proj2), aes(x = proj2, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


stable_p53 <- mafs %>%
  filter(proj %in% enough_lh) %>%
  group_by(wgd, Class, group, Patient, Sample, proj) %>%
  summarize(p53 = sum(Hugo_Symbol == "TP53" & Variant_Classification != 'Silent')) %>%
  filter(Class == 'Low-Hypodiploid') %>%
  mutate(proj2 = ifelse(proj %in% c('ACC', 'KICH'), 'Stereotyped', 'CIN')) %>%
  mutate(proj2 = factor(proj2, levels = c('Stereotyped', 'CIN'))) %>%
  mutate(TP53 = ifelse(p53 > 0, 'Mut', 'WT')) %>%
  mutate(TP53 = factor(TP53, levels = c('WT', 'Mut'))) %>%
  ungroup() %>% 
  ggplot(aes(x = proj2, fill = TP53)) +
  geom_bar(position = "fill") +
  xlab("") +
  labs(subtitle = '(i)') +
  ylab("Proportion of Cases") +
  theme_large_classic() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c(`WT` = "black", `Mut` = "lightsteelblue")) +
  geom_text(data = . %>% count(proj2), aes(x = proj2, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


a1 <- hyp %>%
  group_by(proj) %>%
  mutate(wgd_rate = sum(wgd == "WGD") / n(), poly_rate = sum(Class == "Polyploid") / n(), hypo_rate = sum(Class == "Low-Hypodiploid") / n()) %>%
  ungroup() %>%
  pivot_longer(cols = Buffa_hypoxia_score_pan_cancer:Seigneuric2_hypoxia_score_pan_cancer, names_to = "Score", values_to = "Hypoxia") %>%
  group_by(Score, proj, wgd_rate, poly_rate, hypo_rate) %>%
  summarize(avg_hypoxia = median(Hypoxia), avg_hyp_nonpoly = median(Hypoxia[Class != "Polyploid"]), avg_hyp_nowgd = median(Hypoxia[wgd == "No WGD"]), avg_hyp_nonhypo = median(Hypoxia[Class != "Low-Hypodiploid"])) %>%
  pivot_longer(avg_hypoxia:avg_hyp_nonhypo, names_to = "condition", values_to = "hypoxia")

a2 <- hyp %>%
  group_by(proj) %>%
  mutate(wgd_rate = sum(wgd == "WGD") / n(), poly_rate = sum(Class == "Polyploid") / n(), hypo_rate = sum(Class == "Low-Hypodiploid") / n()) %>%
  ungroup() %>%
  pivot_longer(cols = Buffa_hypoxia_score_pan_cancer:Seigneuric2_hypoxia_score_pan_cancer, names_to = "Score", values_to = "Hypoxia") %>%
  group_by(Score, proj, wgd_rate, poly_rate, hypo_rate) %>%
  summarize(avg_hypoxia = median(Hypoxia), avg_hyp_nonpoly = median(Hypoxia[Class != "Polyploid"]), avg_hyp_nowgd = median(Hypoxia[wgd == "No WGD"]), avg_hyp_nonhypo = median(Hypoxia[Class != "Low-Hypodiploid"])) %>%
  pivot_longer(wgd_rate:hypo_rate, names_to = "ploidy_rate", values_to = "rate")

diffs <- a1 %>%
  left_join(a2, by = c("Score", "proj"), relationship = "many-to-many") %>%
  group_by(Score, ploidy_rate, condition) %>%
  summarize(co = cor(hypoxia, rate), p.val = cor.test(hypoxia, rate)$p.value) %>%
  filter(ploidy_rate == "hypo_rate" & condition %in% c("avg_hyp_nonhypo", "avg_hypoxia") | (ploidy_rate == "poly_rate" & condition %in% c("avg_hyp_nonpoly", "avg_hypoxia"))) %>%
  mutate(ploidy_rate = ifelse(ploidy_rate == "hypo_rate", "Hypodiploidy", "Polyploidy")) %>%
  mutate(condition = ifelse(condition == "avg_hypoxia", "All cases", ifelse(condition == "avg_hyp_nonpoly", "Non-Polyploid Cases", "Non-Hypodiploid Cases"))) %>%
  mutate(condition_pretty = ifelse(condition == 'All cases', 'All cases', ifelse(condition %in% c('Non-Polyploid Cases', 'Non-Hypodiploid Cases'), 'Adjusted', 'what'))) %>%
  ggplot(aes(x = co, y = reorder(sub("(_).*", "", Score), co, FUN = median), colour = condition_pretty)) +
  geom_point(size = 5, aes(shape = p.val < 0.05)) +
  facet_wrap(~ploidy_rate, scales = "free") +
  theme_large() +
  theme(legend.position = "bottom") +
  xlab("Correlation between Hypoxia and Extreme Ploidy") +
  ylab("Hypoxia Score") +
  xlim(c(NA, 1)) +
  labs(subtitle = '(d)') +
  scale_colour_manual(values = c(Adjusted = 'darkgreen', `All cases` = 'darkblue')) + labs(colour = '')

## Mutational Signatures ##


activities <- fread('~/PhD/TCGA/Assignment_Solution/Activities/Assignment_Solution_Activities.txt')

vities <- activities %>% 
  mutate(type = str_extract(Samples, "[^_]+")) %>% 
  mutate(Samples = toupper(sub(".*(tcga*)", "\\1", Samples))) %>% 
  separate(Samples, into = c('tcga', 'bcr', 'pt', NA, NA, NA, NA)) %>% 
  mutate(Patient = paste(tcga, bcr, pt, sep = '-')) %>%
  left_join(tcga_classes, by = c('Patient')) %>%
  filter(!is.na(Class)) %>%
  left_join(clin, by = c('Patient' = 'bcr_patient_barcode')) %>% 
  filter(!is.na(age_at_initial_pathologic_diagnosis)) %>%
  mutate(polyploidy = Class == 'Polyploid', hypo = Class == 'Low-Hypodiploid') %>% 
  left_join(tcga_mu, by = c('Patient', 'proj', 'Class')) %>% 
  filter(!is.na(total_mu)) %>%
  pivot_longer(SBS1:SBS99, names_to = 'Signature', values_to = 'Activity') %>%
  left_join(tcga_purity[, c('Patient', 'cpe')], by = c('Patient')) %>% 
  filter(!is.na(cpe))

sig_hypo_pan <- vities %>%
  group_by(Signature) %>%
  do(tidy(glm(Activity ~ hypo + age_at_initial_pathologic_diagnosis + total_mu + proj + cpe, data = .))) %>%
  filter(term == "hypoTRUE") %>%
  ungroup() %>%
  mutate(bh = p.adjust(p.value, method = "BH"))

mutsigs_plot <- sig_hypo_pan %>% ggplot(aes(x = estimate, y = (-1) * log10(p.value), colour = bh < 0.05)) +
  geom_point() +
  geom_text_repel(aes(label = Signature)) +
  theme_large() +
  ylab("-log10(p)") +
  theme(legend.position = "bottom") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c(`FALSE` = "darkgreen", `TRUE` = "darkblue")) +
  labs(colour = "FDR < 0.05", subtitle = '(d)') 

#fig4 <- (lh_volcano | p53_barchart) / (lh_v_poly_mus | msi_violin) / (hyp_violins | hyp_ploidy) / (hypo_v_ragnum | diffs)
fig4 <- (lh_volcano | p53_barchart) / (lh_v_poly_mus | msi_violin) / (hyp_violins | hyp_ploidy | hypo_v_ragnum) / (stable_hyp | stable_p53)

ggsave(plot = fig4, file = 'paper/hypo_fig5_origins.png', width = 25, height = 30)
ggsave(plot = fig4, file = 'paper/fig5_origins_revised.pdf', width = 25, height = 30)

hypo_v_ragnum_incloutliers <- rag_scores %>% group_by(proj) %>% 
  summarize(wgd_rate = sum(wgd == 'WGD')/n(), hypo_rate = mean(group == 'Low-Hypodiploid'), avg_hyp_nohypos = median(ragnum_score_pan[group != 'Low-Hypodiploid']), avg_hyp_nowgd = median(ragnum_score_pan[wgd == 'No WGD'])) %>%
  ggplot(aes(x = avg_hyp_nohypos, y = hypo_rate)) + 
  geom_point(size = 3) + geom_smooth(method = 'lm', colour = 'black', fill = 'lightsteelblue') + stat_cor() + xlab('Median hypoxia score among non-LH tumours') + ylab('Hypodiploidy Rate') + theme_bw() + geom_label_repel(aes(label = proj))

wgd_v_ragnum <- rag_scores %>% group_by(proj) %>% filter(!proj %in% c('ACC', 'KICH')) %>%
  summarize(wgd_rate = sum(wgd == 'WGD')/n(), hypo_rate = mean(group == 'Low-Hypodiploid'), avg_hyp_nohypos = median(ragnum_score_pan[group != 'Low-Hypodiploid']), avg_hyp_nowgd = median(ragnum_score_pan[wgd == 'No WGD'])) %>%
  ggplot(aes(x = avg_hyp_nowgd, y = wgd_rate)) + 
  geom_point(size = 3) + geom_smooth(method = 'lm', colour = 'black', fill = 'lightsteelblue') + stat_cor() + xlab('Median hypoxia score among non-WGD tumours') + ylab('WGD Rate') + theme_bw() + geom_label_repel(aes(label = proj))

hyp_violins_nondoubledhypos <- rag_scores %>% 
  filter(group == 'Other' | wgd == 'No WGD') %>%
  mutate(Class = factor(Class, levels = c('Near-Haploid', 'Low-Hypodiploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  ggplot(aes(x = Class, y = ragnum_score_pan)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  theme_large() +
  xlab("") +
  ylab("Ragnum Hypoxia Score") +
  theme(legend.position = "none", axis.text.x = element_text(size = bts), axis.title.y = element_text(size = bts + 2)) +
  labs(subtitle = '(a)') +
  scale_fill_manual(values = class_palette) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 

  

msi_violin_nondoubledhypos <- msi_proc %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  ggplot(aes(x = Class, y = mantis_score, colour = Class)) +
  geom_jitter(height = 0) +
  geom_boxplot(alpha = 0.2) +
  facet_wrap(~proj, nrow = 1) +
  theme_large() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c(`Dip` = "darkgreen", `LH` = "darkblue")) +
  ylab("MANTIS score") +
  xlab("") +
  labs(subtitle = '(h)') +
  theme(legend.position = "none") +
  stat_compare_means() +
  geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black") +
  geom_text(data = . %>% count(Class), aes(x = Class, y = -Inf, label = paste0("n = ", n)), vjust = -0.5, size = 5, inherit.aes = FALSE) 


tissue_hypoxia <- rag_scores %>% ggplot(aes(x = reorder(proj, ragnum_score_pan, FUN = median), y = ragnum_score_pan)) +
  geom_boxplot() +
  theme_large() +
  xlab("") +
  ylab("Ragnum Hypoxia Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(subtitle = '(f)') +
  geom_text(data = . %>% count(proj), aes(x = proj, y = -Inf, label = paste0("(", n, ")")), vjust = -0.5, size = 5, inherit.aes = FALSE) 

#check for replication rate confounding

pi_cpm <- fread('~/PhD/Hypodiploidy/all_PIs.tsv') #sent over from server/transcriptome/
pis_cpm <- pi_cpm %>% group_by(V3, cancer) %>% summarize(PI_cpm = median(V4))

#hypodiploidy rate is weakly correlated with replication rate (r = 0.36, p = 0.044)
hypodiploidy_v_PI <- rna_ids %>%
  left_join(pis_cpm, by = c("file_id" = "V3")) %>%
  left_join(tcga_classes, by = c("submitter_id" = "Patient")) %>%
  filter(!is.na(Class)) %>%
  group_by(proj) %>%
  summarize(hypo_rate = mean(group == "Low-Hypodiploid"), med_pi_cpm = median(PI_cpm)) %>%
  filter(!proj %in% c("KICH", "ACC")) %>% #extreme outliers in hypodiploidy rate
  ggplot(aes(x = med_pi_cpm, y = hypo_rate)) +
  geom_point() +
  stat_cor() +
  geom_smooth(method = "lm") +
  geom_text_repel(aes(label = proj)) +
  labs(x = "Median Proliferative Index", y = "Low-hypodiploidy Rate") +
  theme_large() + 
  labs(subtitle = '(e)')

#regression against diploid, controlling for PI and cancer type
rag_scores %>% left_join(rna_ids, by = c('submitter_id')) %>% left_join(pis_cpm, by = c('file_id' = 'V3')) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid')) %>% do(tidy(glm(ragnum_score_pan ~ Class + PI_cpm + proj, data = .))) 
#now against aneuploid
rag_scores %>% left_join(rna_ids, by = c('submitter_id')) %>% left_join(pis_cpm, by = c('file_id' = 'V3')) %>% filter(Class %in% c('Low-Hypodiploid', 'Aneuploid')) %>% do(tidy(glm(ragnum_score_pan ~ Class + PI_cpm + proj, data = .))) 

figS4 <- (total_mus_by_class | mu_per_ploidy_by_class) / (msi_bar | msi_violin_nondoubledhypos) / (tissue_hypoxia) / (hyp_violins_nondoubledhypos | diffs | hypodiploidy_v_PI)
ggsave(plot = figS4, file = 'paper/hypo_supp_fig4_origins.png', width = 25, height = 30)
ggsave(plot = figS4, file = 'paper/hypo_supp_fig4_origins.png', width = 25, height = 30)

## In text ##
#hypoxia scores by class
rag_scores %>% group_by(Class) %>% summarize(median(ragnum_score_pan))

sig_hypo_pan %>% filter(bh < 0.05) 

# p53 mutation rate
mafs %>%
  group_by(wgd, Class, group, Patient, Sample) %>%
  summarize(p53 = sum(Hugo_Symbol == "TP53" & Variant_Classification != 'Silent')) %>%
  group_by(Class) %>%
  mutate(perc_p53 = sum(p53 > 0) / n()) %>% distinct(Class, .keep_all = T)

#p53 mutation rates in LHs by stereotyped vs CIN
mafs %>%
  filter(proj %in% enough_lh) %>%
  group_by(wgd, Class, group, Patient, Sample, proj) %>%
  summarize(p53 = sum(Hugo_Symbol == "TP53" & Variant_Classification != 'Silent')) %>%
  filter(Class == 'Low-Hypodiploid') %>%
  group_by(proj) %>% 
  summarize(perc_p53 = mean(p53 > 0))  %>% mutate(grp = ifelse(proj %in% c('ACC', 'KICH'), proj, 'CIN')) %>% group_by(grp) %>% summarize(median(perc_p53))


## Reviewer Comments ##

#germline associations with hypodiploidy
germline <- fread('~/PhD/TSGs/PCA_pathVar_integrated_filtered_adjusted.tsv') %>%
  filter(CharGer_Classification %in% c('Pathogenic', 'Likely Pathogenic'))

gc <- germline %>% count(bcr_patient_barcode, HUGO_Symbol)

gs <- tcga_classes %>%
  select(Patient) %>%
  left_join(gc, by = c("Patient" = "bcr_patient_barcode")) %>%
  mutate(n_mutations = replace_na(n, 0)) %>%
  complete(Patient, HUGO_Symbol, fill = list(n_mutations = 0)) %>%
  left_join(tcga_classes, by = c("Patient")) %>%
  select(Patient:n_mutations, Class, wgd, proj) %>%
  mutate(mutated = as.numeric(n_mutations > 0))

mu_lh_germline <- gs %>% filter(Class %in% c("Diploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c('Diploid', 'Low-Hypodiploid'))) %>%
  group_by(HUGO_Symbol) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 5) %>%
  do(tidy(glm(mutated ~ Class + proj, family = 'binomial', data = .))) %>% ungroup()

mu_lh_germline %>% filter(term == 'ClassLow-Hypodiploid') %>% mutate(bh = p.adjust(p.value, method = 'BH')) %>% filter(p.value < 0.05)


# do mutational signatures differ between stable and unstable hypodiploids? 
#10d does not

sbs94 <- vities %>%
  filter(
    Signature == "SBS94",
    Class == "Low-Hypodiploid",
    proj %in% enough_lh
  ) %>%
  mutate(
    in_ACC_KICH = proj %in% c("ACC", "KICH"),
    active = Activity > 0
  ) %>%
  count(in_ACC_KICH, active) %>% 
  pivot_wider(names_from = active, values_from = n, values_fill = 0) %>% 
  column_to_rownames(var = "in_ACC_KICH") %>% 
  as.matrix() %>%
  fisher.test() %>%
  tidy()

sbs10d <- vities %>%
  filter(
    Signature == "SBS10d",
    Class == "Low-Hypodiploid",
    proj %in% enough_lh
  ) %>%
  mutate(
    in_ACC_KICH = proj %in% c("ACC", "KICH"),
    active = Activity > 0
  ) %>%
  count(in_ACC_KICH, active) %>% 
  pivot_wider(names_from = active, values_from = n, values_fill = 0) %>% 
  column_to_rownames(var = "in_ACC_KICH") %>% 
  as.matrix() %>%
  fisher.test() %>%
  tidy()

stable_v_unstable_mutsigs_plot <- vities %>%
  filter(Signature %in% c("SBS10d", "SBS94")) %>%
  filter(Class == "Low-Hypodiploid", proj %in% enough_lh) %>%
  mutate(highlight = proj %in% c('ACC', 'KICH')) %>%
  group_by(Signature) %>% #otherwise counts are excessive
  ggplot(aes(x = highlight, fill = Activity > 0)) +
  geom_bar(position = "fill") +
  facet_wrap(~Signature) +
  theme_large() +
  labs(subtitle = '(e)') +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("darkgreen", "darkblue")) +
  xlab("Stable Hypodiploid") +
  geom_text(data = . %>% count(highlight), aes(x = highlight, y = 0, label = paste0("n = ", n)), vjust = 1, size = 5, inherit.aes = FALSE) +
  ylab('Proportion of Cases with Activity > 0')


# Copy Number Signatures from Steele et al (2022)

cn_sigs <- fread('Steele_CNsigs_Attributions.csv')

enrichment_df_all <- cn_sigs %>%
  left_join(tcga_classes, by = c("Sample" = "Patient")) %>%
  filter(!is.na(Class)) %>%
  add_count(Class, name = 'n_samples') %>%
  pivot_longer(CN1:CN21, names_to = "Signature", values_to = "Activity") %>%
  group_by(Class, n_samples, Signature) %>%
  summarise(fraction_active = mean(Activity > 0), .groups = "drop") %>%
  group_by(Signature) %>%
  mutate(global_fraction = mean(fraction_active),
         enrichment = fraction_active - global_fraction) %>%
  ungroup() %>%
  mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class)))

cnrich_all <- enrichment_df_all %>% ggplot(aes(x = paste0(Class, '\n(', n_samples, ')'), y = Signature, fill = enrichment)) +
  geom_tile() +
  labs(subtitle = "All cases") +
  geom_text(aes(label = round(enrichment, 2))) + xlab('Class') + theme_large() + labs(subtitle = '(f)')

enrichment_df_nondoubledhypos <- cn_sigs %>%
  left_join(tcga_classes, by = c("Sample" = "Patient")) %>%
  filter(!is.na(Class)) %>%
  filter(group == 'Other' | wgd == 'No WGD') %>%
  add_count(Class, name = 'n_samples') %>%
  pivot_longer(CN1:CN21, names_to = "Signature", values_to = "Activity") %>%
  group_by(Class, n_samples, Signature) %>%
  summarise(fraction_active = mean(Activity > 0), .groups = "drop") %>%
  group_by(Signature) %>%
  mutate(global_fraction = mean(fraction_active), #should maybe do 'NOT in the group' to avoid penalising big groups
         enrichment = fraction_active - global_fraction) %>%
  ungroup() %>%
  mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class)))


cnrich_ndh <- enrichment_df_nondoubledhypos %>% ggplot(aes(x = paste0(Class, '\n(', n_samples, ')'), y = Signature, fill = enrichment)) +
  geom_tile() +
  labs(subtitle = "Excluding doubled hypodiploids") +
  geom_text(aes(label = round(enrichment, 2))) + xlab('Class') + theme_large()

cnrich_onlystablehypos <- cn_sigs %>%
  left_join(tcga_classes, by = c("Sample" = "Patient")) %>%
  filter(!is.na(Class)) %>%
  filter(Class != 'Low-Hypodiploid' | proj %in% c('ACC', 'KICH')) %>%
  add_count(Class, name = 'n_samples') %>%
  pivot_longer(CN1:CN21, names_to = "Signature", values_to = "Activity") %>%
  group_by(Class, n_samples, Signature) %>%
  summarise(fraction_active = mean(Activity > 0), .groups = "drop") %>%
  group_by(Signature) %>%
  mutate(global_fraction = mean(fraction_active),
         enrichment = fraction_active - global_fraction) %>%
  ungroup() %>%   mutate(Class = ifelse(Class == 'Low-Hypodiploid', 'LH', ifelse(Class == 'Near-Haploid', 'NH', Class))) %>%
  ggplot(aes(x = paste0(Class, '\n(', n_samples, ')'), y = Signature, fill = enrichment)) + geom_tile() + geom_text(aes(label = round(enrichment, 2))) + labs(subtitle = 'Excluding unstable hypodiploids') + xlab('Class') + theme_large()

CN_enrichment <- cnrich_all | cnrich_ndh | cnrich_onlystablehypos

## Hypoxia: adjusting for cancer type and purity
#all
rag_scores %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('submitter_id' = 'Patient')) %>% filter(!is.na(cpe)) %>% do(tidy(glm(ragnum_score_pan ~ Class + proj + cpe, data = .)))

#non doubled
rag_scores %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('submitter_id' = 'Patient')) %>% filter(group == 'Other' | wgd == 'No WGD') %>% filter(!is.na(cpe)) %>% do(tidy(glm(ragnum_score_pan ~ Class + proj + cpe, data = .)))

#stereotyped hypodiploids have lower hypoxia, controlling for purity
rag_scores %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('submitter_id' = 'Patient')) %>% filter(!is.na(cpe)) %>% filter(Class == 'Low-Hypodiploid', proj %in% enough_lh) %>% mutate(highlight = proj %in% c('ACC', 'KICH')) %>% do(tidy(glm(ragnum_score_pan ~ highlight + cpe, data = .)))

## hypoxia by class, adjusted for cancer type and purity, regression plot:
linreg_hypoxia <- rag_scores %>% left_join(tcga_purity[, c('Patient', 'cpe')], by = c('submitter_id' = 'Patient')) %>% 
  filter(!is.na(cpe)) %>% 
  mutate(Class = factor(Class, levels = c('Low-Hypodiploid', 'Near-Haploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  do(tidy(glm(ragnum_score_pan ~ Class + proj + cpe, data = .), conf.int = T)) %>% 
  mutate(term = ifelse(term == 'cpe', 'Purity', term)) %>% 
  filter(term != '(Intercept)') %>%
  filter(!startsWith(term, 'proj')) %>% ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high, colour = p.value < 0.05)) +
  geom_errorbar() +
  geom_point() +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) +
  coord_flip() +
  theme_large() + xlab('Term') + theme(legend.position = 'bottom') + scale_colour_manual(values = c('grey50', 'black')) +
  labs(subtitle = '(b) Hypoxia vs Low-Hypodiploid')

#LR controlling for level of aneuploidy
rag_scores %>% left_join(tcga_purity[, c('Patient', 'cpe', 'n_disomic_nosex')], by = c('submitter_id' = 'Patient')) %>% 
  filter(!is.na(cpe)) %>% 
  mutate(Class = factor(Class, levels = c('Low-Hypodiploid', 'Near-Haploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  mutate(LH = Class == 'Low-Hypodiploid') %>%
  do(tidy(glm(ragnum_score_pan ~ LH + proj + cpe + n_disomic_nosex, data = .), conf.int = T))

rag_adj_aneuploidy_plot <- rag_scores %>% left_join(tcga_purity[, c('Patient', 'cpe', 'n_disomic_nosex')], by = c('submitter_id' = 'Patient')) %>% 
  filter(!is.na(cpe)) %>% 
  mutate(Class = factor(Class, levels = c('Low-Hypodiploid', 'Near-Haploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  mutate(LH = Class == 'Low-Hypodiploid') %>%
  do(tidy(glm(ragnum_score_pan ~ LH + proj + cpe + n_disomic_nosex, data = .), conf.int = T)) %>%
filter(term != '(Intercept)', !startsWith(term, "proj")) %>%
  ggplot(aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar() +
  geom_point() +
  geom_text_repel(aes(label = paste0("p = ", scales::scientific(p.value, digits = 3)))) +
  coord_flip() +
  theme_large() +
  xlab("Term") +
  theme(legend.position = "bottom") +
  ggtitle("(c) Hypoxia vs LH") + geom_hline(yintercept = 0, colour = 'red')

## new supplementary figures

figS4_mutations <- (total_mus_by_class | mu_per_ploidy_by_class) / (aneu_volcano | mutsigs_plot | stable_v_unstable_mutsigs_plot) / CN_enrichment / (msi_bar | msi_violin_nondoubledhypos) 
ggsave(plot = figS4_mutations, file = 'paper/supp_fig7_origins_mutations_revised.pdf', width = 25, height = 30)

figS4_hypoxia <- ((hyp_violins_nondoubledhypos) | (linreg_hypoxia | rag_adj_aneuploidy_plot)) / (diffs | hypodiploidy_v_PI) / tissue_hypoxia
ggsave(plot = figS4_hypoxia, file = 'paper/supp_fig8_origins_hypoxia_revised.pdf', width = 25, height = 30)
