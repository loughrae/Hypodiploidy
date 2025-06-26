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

mu_perploidy_reg <- tcga_mu %>%
  filter(Class %in% c("Diploid", "Low-Hypodiploid")) %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  ungroup() %>%
  do(tidy(lm(mu_per_ploidy ~ Class + proj, data = .))) # near diploids have higher mutations_per_ploidy after controlling for cancer type but that's actually a MEAN bc of hypermutation

total_mus_by_class <- tcga_mu %>%
  left_join(tcga_classes, by = c("Class", "ploidy", "proj", "Patient")) %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  ggplot(aes(x = Class, y = total_mu, fill = Class)) +
  geom_violin() +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Diploid", "Low-Hypodiploid"), c("Diploid", "Near-Haploid"))) +
  theme_large() +
  xlab("") +
  labs(subtitle = '(a)') +
  scale_y_continuous(trans = "log10") +
  ylab("Total Mutations (log10)") +
  scale_fill_manual(values = class_palette) +
  theme(legend.position = "none")

mu_per_ploidy_by_class <- tcga_mu %>%
  left_join(tcga_classes, by = c("Class", "ploidy", "proj", "Patient")) %>%
  mutate(mu_per_ploidy = total_mu / ploidy) %>%
  ggplot(aes(x = Class, y = mu_per_ploidy, fill = Class)) +
  geom_violin() +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Diploid", "Low-Hypodiploid"), c("Diploid", "Near-Haploid"))) +
  theme_large() +
  labs(subtitle = '(b)') +
  xlab("") +
  scale_y_continuous(trans = "log10") +
  ylab("Ploidy-Corrected Mutation Count (log10)") +
  scale_fill_manual(values = class_palette) +
  theme(legend.position = "none")

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
  pivot_longer(-c(Patient, total_mu, total_non_syn_mu, proj, Class, ploidy), names_to = "gene", values_to = "mutated", values_drop_na = F) #NB: the mutate across above is necessary for this to work

mu_lh <- status %>%
  filter(Class %in% c("Diploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c('Diploid', 'Low-Hypodiploid'))) %>%
  group_by(gene) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 5) %>%
  do(tidy(glm(mutated ~ Class + total_non_syn_mu + proj, family = 'binomial', data = .))) %>% ungroup()

mu_poly <- status %>%
  filter(Class %in% c("Diploid", "Polyploid")) %>%
  mutate(Class = factor(Class, levels = c('Diploid', 'Polyploid'))) %>%
  group_by(gene) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 5) %>%
  do(tidy(glm(mutated ~ Class + total_non_syn_mu + proj, family = 'binomial', data = .))) %>% 
  ungroup()

mu_hap_v_lh <- status %>%
  filter(Class %in% c("Near-Haploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c('Low-Hypodiploid', 'Near-Haploid'))) %>%
  group_by(gene) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 3) %>% #lower here note
  do(tidy(glm(mutated ~ Class + total_non_syn_mu + proj, family = 'binomial', data = .))) %>% 
  ungroup()

mu_hap_v_lh %>% filter(term == 'ClassNear-Haploid') %>% mutate(bh = p.adjust(p.value, method = 'BH')) %>% filter(bh < 0.05) #0 genes
mu_hap_v_lh %>% filter(term == 'ClassNear-Haploid') %>% mutate(bh = p.adjust(p.value, method = 'BH')) %>% filter(p.value < 0.05) #10 genes

mu_lh_v_aneuploid <- status %>%
  filter(Class %in% c("Aneuploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c('Aneuploid', 'Low-Hypodiploid'))) %>%
  group_by(gene) %>%
  mutate(n_mutated = sum(mutated == 1)) %>%
  filter(n_mutated >= 5) %>% 
  do(tidy(glm(mutated ~ Class + total_non_syn_mu + proj, family = 'binomial', data = .))) %>% 
  ungroup()

mu_lh_v_aneuploid %>% filter(term == 'ClassLow-Hypodiploid') %>% mutate(bh = p.adjust(p.value, method = 'BH')) %>% filter(bh < 0.05) 


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
  stat_cor() +
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
  ggplot(aes(x = Class, fill = p53 > 0)) +
  geom_bar(position = "fill") +
  xlab("") +
  labs(fill = "Mutated TP53", subtitle = '(b)') +
  ylab("Proportion of Cases") +
  theme_large_classic() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c(`FALSE` = "black", `TRUE` = "lightsteelblue"))

## fig4d: Microsatellite instability
msi <- fread('~/Downloads/mantis_msi.csv') %>% clean_names()

msi_proc <- tcga_classes %>% 
  left_join(msi, by = c("Patient" = "case_id")) %>%
  filter(!is.na(mantis_score)) %>%
  filter(proj %in% c("UCEC", "COAD", "STAD", "ACC")) %>%
  mutate(proj = factor(proj, levels = c("COAD", "UCEC", "STAD", "ACC"))) %>%
  filter(Class %in% c("Diploid", "Low-Hypodiploid")) %>%
  mutate(Class = factor(Class, levels = c("Diploid", "Low-Hypodiploid"))) %>%
  mutate(Class = ifelse(Class == "Diploid", "Dip", "LH")) 

msi_bar <- msi_proc %>%
  ggplot(aes(x = Class, fill = mantis_score > 0.4)) +
  geom_bar(position = "fill") +
  facet_wrap(~proj, nrow = 1) +
  theme_large() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c(`FALSE` = "black", `TRUE` = "lightsteelblue")) +
  ylab("Proportion of Cases") +
  xlab("") +
  labs(fill = "MANTIS > 0.4", subtitle = '(c)')

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
  geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black") 

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
  theme_bw() +
  geom_label_repel(aes(label = proj))

bts <- 12

pairwise_comparisons <- list(
  c("Diploid", "Polyploid"),
  c("Diploid", "Low-Hypodiploid"),
  c("Polyploid", "Low-Hypodiploid"),
  c("Aneuploid", "Low-Hypodiploid")
)

hyp_violins <- rag_scores %>% 
  mutate(Class = factor(Class, levels = c('Near-Haploid', 'Low-Hypodiploid', 'Diploid', 'Aneuploid', 'Polyploid'))) %>%
  ggplot(aes(x = Class, y = ragnum_score_pan)) +
  geom_violin(aes(fill = Class)) +
  geom_boxplot(alpha = 0.2) +
  theme_bw() +
  xlab("") +
  ylab("Ragnum Hypoxia Score") +
  theme(legend.position = "none", axis.text.x = element_text(size = bts), axis.title.y = element_text(size = bts + 2)) +
  ggtitle("(e)") +
  scale_fill_manual(values = class_palette) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) 


hyp_ploidy <- rag_scores %>% ggplot(aes(x = ploidy, y = ragnum_score_pan)) +
  geom_point(alpha = 0.2) +
  geom_smooth(colour = 'black', fill = 'lightsteelblue') +
  theme_minimal() +
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
  stat_compare_means() +
  theme_large() +
  theme(legend.position = "none") +
  xlab("") +
  labs(subtitle = '(h)') +
  ylab("Ragnum hypoxia score") +
  scale_fill_manual(values = c(`Stereotyped` = "#3288bd", `CIN` = "#1b9e77"))

stable_p53 <- mafs %>%
  filter(proj %in% enough_lh) %>%
  group_by(wgd, Class, group, Patient, Sample, proj) %>%
  summarize(p53 = sum(Hugo_Symbol == "TP53" & Variant_Classification != 'Silent')) %>%
  filter(Class == 'Low-Hypodiploid') %>%
  mutate(proj2 = ifelse(proj %in% c('ACC', 'KICH'), 'Stereotyped', 'CIN')) %>%
  mutate(proj2 = factor(proj2, levels = c('Stereotyped', 'CIN'))) %>%
  ggplot(aes(x = proj2, fill = p53 > 0)) +
  geom_bar(position = "fill") +
  xlab("") +
  labs(fill = "Mutated TP53", subtitle = '(i)') +
  ylab("Proportion of Cases") +
  theme_large_classic() +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c(`FALSE` = "black", `TRUE` = "lightsteelblue"))

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
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Correlation between Hypoxia and Extreme Ploidy") +
  ylab("Hypoxia Score") +
  xlim(c(NA, 1)) +
  labs(subtitle = '(g)') +
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
  left_join(clin, by = c('Patient' = 'submitter_id')) %>% 
  filter(!is.na(age_at_index)) %>%
  mutate(polyploidy = Class == 'Polyploid', hypo = Class == 'Low-Hypodiploid') %>% 
  left_join(tcga_mu, by = c('Patient', 'proj', 'Class')) %>% 
  filter(!is.na(total_mu)) %>%
  pivot_longer(SBS1:SBS99, names_to = 'Signature', values_to = 'Activity')

sig_hypo_pan <- vities %>%
  group_by(Signature) %>%
  do(tidy(glm(Activity ~ hypo + age_at_index + total_mu + proj, data = .))) %>%
  filter(term == "hypoTRUE") %>%
  ungroup() %>%
  mutate(bh = p.adjust(p.value, method = "BH"))

#fig4 <- (lh_volcano | p53_barchart) / (lh_v_poly_mus | msi_violin) / (hyp_violins | hyp_ploidy) / (hypo_v_ragnum | diffs)
fig4 <- (lh_volcano | p53_barchart) / (lh_v_poly_mus | msi_violin) / (hyp_violins | hyp_ploidy | hypo_v_ragnum) / (stable_hyp | stable_p53)

ggsave(plot = fig4, file = 'paper/fig4_test.png', width = 25, height = 30)

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
  theme_bw() +
  xlab("") +
  ylab("Ragnum Hypoxia Score") +
  theme(legend.position = "none", axis.text.x = element_text(size = bts), axis.title.y = element_text(size = bts + 2)) +
  labs(subtitle = '(f)') +
  scale_fill_manual(values = class_palette) +
  stat_compare_means(comparisons = pairwise_comparisons, method = "wilcox.test", label = "p.format", tip.length = 0.02) 

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
  labs(subtitle = '(d)') +
  theme(legend.position = "none") +
  stat_compare_means() +
  geom_hline(yintercept = 0.4, linetype = "dashed", colour = "black") 

tissue_hypoxia <- rag_scores %>% ggplot(aes(x = reorder(proj, ragnum_score_pan, FUN = median), y = ragnum_score_pan)) +
  geom_boxplot() +
  theme_large() +
  xlab("") +
  ylab("Ragnum Hypoxia Score") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(subtitle = '(e)')

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
  labs(subtitle = '(h)')

#regression against diploid
rag_scores %>% left_join(rna_ids, by = c('submitter_id')) %>% left_join(pis_cpm, by = c('file_id' = 'V3')) %>% filter(Class %in% c('Low-Hypodiploid', 'Diploid')) %>% do(tidy(glm(ragnum_score_pan ~ Class + PI_cpm, data = .))) 
#regression against aneuploid
rag_scores %>% left_join(rna_ids, by = c('submitter_id')) %>% left_join(pis_cpm, by = c('file_id' = 'V3')) %>% do(tidy(glm(ragnum_score_pan ~ Class + PI_cpm, data = .)))

figS4 <- (total_mus_by_class | mu_per_ploidy_by_class) / (msi_bar | msi_violin_nondoubledhypos) / (tissue_hypoxia) / (hyp_violins_nondoubledhypos | diffs | hypodiploidy_v_PI)
ggsave(plot = figS4, file = 'paper/supp_fig4.png', width = 25, height = 30)

## In text ##
#hypoxia scores by class
rag_scores %>% group_by(Class) %>% summarize(median(ragnum_score_pan))

# no enriched mutational signatures after multiple testing correction
sig_hypo_pan %>% filter(bh < 0.05) 

# p53 mutation rate
mafs %>%
  group_by(wgd, Class, group, Patient, Sample) %>%
  summarize(p53 = sum(Hugo_Symbol == "TP53" & Variant_Classification != 'Silent')) %>%
  group_by(Class) %>%
  mutate(perc_p53 = sum(p53 > 0) / n()) %>% distinct(Class, .keep_all = T)

#p53 mutation rates in LHs by group
mafs %>%
  filter(proj %in% enough_lh) %>%
  group_by(wgd, Class, group, Patient, Sample, proj) %>%
  summarize(p53 = sum(Hugo_Symbol == "TP53" & Variant_Classification != 'Silent')) %>%
  filter(Class == 'Low-Hypodiploid') %>%
  group_by(proj) %>% 
  summarize(perc_p53 = mean(p53 > 0))  %>% mutate(grp = ifelse(proj %in% c('ACC', 'KICH'), proj, 'CIN')) %>% group_by(grp) %>% summarize(median(perc_p53))