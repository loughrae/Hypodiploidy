source('setup.R')

tcga_classes <- fread('TCGA_ploidy_classes.tsv')
tcga_losses <- fread('tcga_losses.tsv')
mitcn_meta <- fread('mitcn_meta.tsv')
mitcn_chr_summary <- fread('mitcn_chr_summary.tsv')

chr_length <- arms %>% group_by(chr) %>% summarize(first = min(st), last = max(en))

mit_arms <- read.table('mitcn_arms.bed', header = F) %>%
  mutate(len = V3 - V2) %>% 
  rename(chr = V1, start = V2, end = V3, id = V4, CN = V5, arm = V9) %>%
  mutate(chr_arm = paste0(sub('chr', '', chr), arm)) %>% #has no X or Y
  left_join(mitcn_meta, by = c('id'))

fa <- fread('fasc_arms.bed') %>% 
  mutate(len = V3 - V2) %>% 
  rename(chr = V1, start = V2, end = V3, GDC_Aliquot = V4, CN = V5, mCN = V6, proj = V7, arm = V11) %>%
  mutate(loh_indicator = ifelse(mCN == 0, 1, 0)) %>% 
  mutate(chr_arm = paste0(sub('chr', '', chr), arm)) %>%
  filter(GDC_Aliquot %in% tcga_classes$GDC_Aliquot)

enough_lh <- tcga_classes %>% filter(group == 'Low-Hypodiploid') %>% count(proj) %>% filter(n >= 15) %>% pull(proj) 

## Chromosome Features ## 
tuson <- fread('~/Downloads/tuson_chrom.csv', skip = 45) %>% clean_names()
ds <- fread('~/Downloads/Collins_rCNV_2022.dosage_sensitivity_scores.tsv') %>% clean_names() %>% rename(gene = number_gene)
ds_rates <- ga %>% left_join(ds, by = c('V4' = 'gene')) %>% filter(!is.na(p_haplo)) %>% group_by(V1) %>% summarize(sum_p_haplo = sum(p_haplo), mean_p_haplo = mean(p_haplo)) 

chr_info <- chr_length %>%
  filter(!chr %in% c('chrX', 'chrY')) %>%
  mutate(length = last / 1000000) %>%
  mutate(location = case_when(
    chr %in% c("chr1", "chr8", "chr11", "chr17", "chr19") ~ "Central",
    chr %in% c("chr13", "chr15", "chr18", "chr21", "chrX") ~ "Peripheral",
    TRUE ~ "Other"
  )) %>%
  mutate(match_chr = as.numeric(sub("chr", "", chr))) %>%
  left_join(tuson, by = c("match_chr" = "chromosome")) %>%
  left_join(ds_rates, by = c("chr" = "V1")) 

print('Making Figure 2')

mit_arm_losses <- mit_arms %>%
  group_by(id, case_id, Sex, Subclone, n_chr_nosex, group, any_lows, any_hyper, n_clones, chr, chr_arm, ploidy) %>%  
  summarize(arm_somy = round(sum(CN*len)/sum(len),0), loss_perc = sum(len[CN < 2])/sum(len)) %>%
  mutate(loss = ifelse(loss_perc > 0.9, 'Lost', 'Retained')) %>%
  mutate(Class = group, proj = 'ALL') %>%
  ungroup() 

mit_chr_losses <- mitcn_chr_summary %>%
  mutate(Class = group, proj = 'ALL') %>%
  ungroup() 

mit_chr_losses_tojoin <- mit_chr_losses %>% 
  mutate(Class = ifelse(ploidy > 2.7, 'Polyploid', Class)) %>% 
  select(id, Class, chr, proj, loss, chr_somy) 

mit_arm_losses_tojoin <- mit_arm_losses %>% 
  mutate(Class = ifelse(ploidy > 2.7, 'Polyploid', Class)) %>% 
  select(id, Class, chr, proj, loss, chr_arm) 

tcga_arm_losses <- fa %>% 
  mutate(prod = CN*len, loh_len = loh_indicator*len) %>%
  group_by(GDC_Aliquot, proj) %>%
  mutate(ploidy = sum(prod)/sum(len)) %>%
  group_by(GDC_Aliquot, chr, chr_arm, proj, ploidy) %>%
  summarize(loh = sum(loh_len)/sum(len), chr_somy = round(sum(prod)/sum(len), 0)) %>% #calculate proportion of LOH per chromosome per sample
  mutate(loss = ifelse(loh > 0.9, 'Lost', 'Retained')) %>% 
  mutate(min_somy = ifelse(loss == 'Retained', 2, pmin(1, chr_somy))) %>% #base_somy = 2 because excluding the sex chromosomes
  ungroup()


joint_chr_losses <- tcga_losses %>%
  left_join(tcga_classes, by = c("GDC_Aliquot", "proj")) %>%
  select(GDC_Aliquot, Class, chr, proj, loss, chr_somy) %>%
  rename(id = GDC_Aliquot) %>%
  rbind(mit_chr_losses_tojoin) %>%
  ungroup() %>%
  left_join(chr_info, by = c('chr'))

joint_arm_losses <- tcga_arm_losses %>%
  left_join(tcga_classes, by = c("GDC_Aliquot", "proj")) %>%
  select(GDC_Aliquot, Class, chr, chr_arm, proj, loss) %>%
  rename(id = GDC_Aliquot) %>%
  rbind(mit_arm_losses_tojoin) %>%
  ungroup() %>%
  left_join(chr_info, by = c('chr'))

## Calculate loss rates ##

joint_chr_loss_rates <- joint_chr_losses %>%
  group_by(Class, chr, length, location, chrom_tsg_og_score, sum_p_haplo, mean_p_haplo) %>%
  summarize(loss_rate = mean(loss == 'Lost')) %>%
  ungroup() 

hypo_arm_loss_rates_tissue <- joint_arm_losses %>%
  filter(Class == 'Low-Hypodiploid') %>%
  group_by(chr_arm, length, location, chrom_tsg_og_score, sum_p_haplo, mean_p_haplo, proj) %>%
  summarize(loss_rate = mean(loss == 'Lost')) %>%
  ungroup() 

hypo_chr_loss_rates <- joint_chr_loss_rates %>% filter(Class == 'Low-Hypodiploid')

#Loss Rates graph
joint_loss_rates_inequality <- hypo_arm_loss_rates_tissue %>%
  filter(proj %in% c(enough_lh, 'ALL')) %>%
  filter(!chr_arm %in% acrocentric_arms) %>%
  group_by(proj) %>% 
  mutate(sha = Entropy(loss_rate), gi = Gini(loss_rate)) 

joint_loss_rates_graph <- joint_loss_rates_inequality %>%
  ggplot(aes(x = reorder(proj, gi), y = factor(chr_arm, levels = chr_levels), fill = loss_rate)) + 
  geom_tile() + scale_fill_viridis_c() + 
  xlab('') + ylab('') + theme_large(base_size = 14) + theme(legend.position = 'right') + labs(fill = 'Lost') + ggtitle('(a)') + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## Clustering heatmap

lossprofiles <- joint_arm_losses %>%
  filter(Class == 'Low-Hypodiploid') %>%
  mutate(loss = as.numeric(loss == "Lost")) %>%
  rename(group = proj) %>%
  select(id, chr_arm, loss, group) %>%
  pivot_wider(names_from = chr_arm, values_from = loss) %>%
  filter(group %in% c(enough_lh, "ALL"))

idg <- lossprofiles %>%
  filter(!group %in% c("NH-ALL", "NH-TCGA")) %>%
  select(id, group)

colours <- distinct_palette(n = length(unique(idg$group)), pal = "brewerPlus")[1:length(unique(idg$group))]
colourdf <- data.frame(Tissue = unique(idg$group), colour = colours)

idg2 <- idg %>%
  left_join(colourdf, by = c("group" = "Tissue")) %>%
  column_to_rownames(var = "id")

toheat <- lossprofiles %>%
  filter(group %in% idg$group) %>%
  select(-group) %>%
  column_to_rownames(var = "id") %>%
  as.matrix() %>%
  proxy::dist(by_rows = T, method = "Jaccard") %>%
  as.matrix()

htmap <- ComplexHeatmap::pheatmap(toheat, col = rev(viridis::viridis(15, option = "D")), show_rownames = F, show_colnames = F, annotation_row = idg2 %>% select(group), annotation_names_row = F, heatmap_legend_param = list(title = "Distance"), annotation_colors = list(group = deframe(distinct(idg2)))) %>% as.ggplot() + ggtitle('(b)')

## jaccard tiles
jacc_tile <- toheat %>%
  as.data.frame() %>%
  rownames_to_column(var = "Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "Distance", values_drop_na = T) %>%
  left_join(idg, by = c("Sample1" = "id")) %>%
  left_join(idg, by = c("Sample2" = "id"), suffix = c("1", "2")) %>%
  filter(Sample1 != Sample2) %>% # remove self-comparisons
  mutate(Similarity = 1 - Distance) %>%
  mutate(comparison = paste(group1, group2, sep = "_"), same = group1 == group2) %>%
  group_by(group1, group2, same) %>%
  summarize(medsim = median(Similarity), std = sd(Similarity)) %>%
  ggplot(aes(x = group1, y = group2, fill = medsim)) +
  geom_tile() +
  geom_text(aes(label = round(medsim, 2))) +
  scale_fill_continuous(type = "viridis") +
  xlab("") +
  ylab("") +
  labs(fill = "Median\nSimilarity") +
  theme_large() +
  theme(legend.position = "right") +
  ggtitle('(c)')

## Nuclear Location
location_palette <- c(
  "Central" = "#1b9e77",  
  "Peripheral"       = "#3288bd",  
  "Other"           = "#a6bddb"  
)

set.seed(123)

nuclear_location <- hypo_chr_loss_rates %>%
  ggplot(aes(x = location, y = loss_rate, fill = location)) + 
  geom_boxplot(alpha = 0.3) +
  geom_jitter(height = 0) + 
  geom_text_repel(aes(label = chr)) +
  theme_large() +
  theme(legend.position = "none") +
  stat_compare_means(hjust = -0.7, size = 4) +
  xlab('Nuclear Location') + 
  ylab('Chromosome Loss Rate') +
  ggtitle('(e)')

### length ###

length_cor <- hypo_chr_loss_rates %>%
  ggplot(aes(x = round(length, 2), y = loss_rate)) +
  geom_point(size = 5, colour = "darkblue") +
  geom_smooth(method = "lm", fill = "lightsteelblue", colour = "black") +
  stat_cor(vjust = 1.5, size = 8) +
  geom_text_repel(size = 6, aes(label = sub("chr", "", chr))) +
  theme_large() +
  xlab("Chromosome Length (Mb)") +
  ylab("Chromosome Loss Rate") +
  ggtitle("(c)")


## tuson ##
tuson_loss <- hypo_chr_loss_rates %>%
  ggplot(aes(x = chrom_tsg_og_score, y = loss_rate)) +
  geom_point(size = 5, colour = 'darkblue') +
  geom_smooth(method = 'lm', fill = 'lightsteelblue', colour = 'black') +
  stat_cor(size = 8) +
  theme_large() +
  xlab("CHROM score (TSG - OG)") +
  ylab("Chromosome Loss Rate") +
  geom_text_repel(size = 6, aes(label = sub('chr', '', chr))) +
  ggtitle('(f)')


## Dosage Sensitivity ##
ds_sum_p_haplo <- hypo_chr_loss_rates %>% 
  ggplot(aes(x = sum_p_haplo, y = loss_rate)) +
  geom_point(size = 5, colour = 'darkblue') +
  geom_smooth(method = 'lm', fill = 'lightsteelblue', colour = 'black') +
  stat_cor(size = 8, label.x = Inf, label.y = Inf, hjust = 1.1, vjust = 1.5) +  theme_large() +
  xlab("Haplosensitivity score (sum)") +
  ylab("Chromosome Loss Rate") +
  geom_text_repel(size = 6, aes(label = sub('chr', '', chr))) +
  ggtitle('(d)')



fig2 <- (joint_loss_rates_graph | htmap) / (length_cor | ds_sum_p_haplo ) / ( nuclear_location | tuson_loss)
#(joint_loss_rates_graph | htmap) /  / ( tuson_loss | nuclear_location)

ggsave(plot = fig2, file = 'paper/fig2_paper.png', width = 25, height = 30)

## Supp Fig 2 ##

nh_all <- mit_arm_losses %>%
  filter(Class %in% c("Near-Haploid", "Low-Hypodiploid")) %>%
  group_by(Class, chr_arm) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  ggplot(aes(x = Class, y = factor(chr_arm, levels = chr_levels), fill = loss_rate)) +
  geom_tile() +
  scale_fill_viridis_c() +
  xlab("") +
  ylab("") +
  theme_large(base_size = 14) +
  theme(legend.position = "right") +
  labs(fill = "Lost") +
  ggtitle("(a)")

location_no17 <- hypo_chr_loss_rates %>%
  filter(chr != 'chr17') %>%
  ggplot(aes(x = location, y = loss_rate, fill = location)) + 
  geom_boxplot(alpha = 0.3) +
  geom_jitter(height = 0) + 
  geom_text_repel(aes(label = chr)) +
  theme_large() +
  theme(legend.position = "none") +
  stat_compare_means(hjust = -0.7, size = 4) +
  xlab('Nuclear Location') + 
  ylab('Chromosome Loss Rate') +
  ggtitle('(d)')

ds_loss_mean <- hypo_chr_loss_rates %>% 
  ggplot(aes(x = mean_p_haplo, y = loss_rate)) +
  geom_point(size = 5, colour = 'darkblue') +
  geom_smooth(method = 'lm', fill = 'lightsteelblue', colour = 'black') +
  stat_cor(size = 8, label.x = Inf, label.y = Inf, hjust = 1.1, vjust = 1.5) +  theme_large() +
  xlab("Haplosensitivity score (mean)") +
  ylab("Chromosome Loss Rate") +
  geom_text_repel(size = 6, aes(label = sub('chr', '', chr))) 

lr_overall <- hypo_chr_loss_rates %>%
  do(tidy(glm(loss_rate ~ sum_p_haplo + location + chrom_tsg_og_score + length, data = .)))

lr_overall_meands <- hypo_chr_loss_rates %>%
  filter(Class == 'Low-Hypodiploid') %>%
  do(tidy(glm(loss_rate ~ mean_p_haplo + location + chrom_tsg_og_score + length, data = .)))

lr_tissues <- joint_chr_losses %>%
  filter(Class == 'Low-Hypodiploid') %>%
  filter(proj %in% c(enough_lh, 'ALL')) %>%
  group_by(chr, proj, sum_p_haplo, location, chrom_tsg_og_score, length) %>% 
  summarize(loss_rate = mean(loss == 'Lost')) %>%
  group_by(proj) %>%
  do(tidy(glm(loss_rate ~ sum_p_haplo + location + chrom_tsg_og_score + length, data = .)))

lr_tissues_figure <- lr_tissues %>%
  mutate(term = case_when(term == 'sum_p_haplo' ~ 'p(HI)', term == 'chrom_tsg_og_score' ~ 'TSG - OG', TRUE ~ term)) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = ifelse(p.value >= 0.05, 'n.s.', ifelse(p.value < 0.05 & estimate > 0, '+', '-'))) %>%
  ggplot(aes(x = proj, y = term, colour = p.value < 0.05)) +
  geom_point(size = 4) +
  theme_large_classic() +
  theme(legend.position = "bottom") +
  labs(y = "", x = "Cancer Type") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle('(e)')


degree_aneu <- joint_chr_losses %>% #group 0 will be NaN, that's fine
  group_by(id, Class) %>%
  mutate(n_losses = sum(loss == "Lost")) %>%
  filter(proj != 'ALL') %>% #ALL losses are calculated differently (no LOH)
  filter(Class %in% c('Low-Hypodiploid', 'Aneuploid', 'Near-Haploid')) %>%
  mutate(losses_group = ifelse(Class %in% c("Low-Hypodiploid", "Near-Haploid"), "6+", as.character(n_losses)))

lr_aneuploidy <- degree_aneu %>%
  group_by(chr, losses_group, sum_p_haplo, location, chrom_tsg_og_score, length) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  group_by(losses_group) %>%
  do(tidy(glm(loss_rate ~ sum_p_haplo + location + chrom_tsg_og_score + length, data = .))) #all TSGs

loss_group_figure <- lr_aneuploidy %>%
  mutate(term = case_when(term == 'sum_p_haplo' ~ 'p(HI)', term == 'chrom_tsg_og_score' ~ 'TSG - OG', TRUE ~ term)) %>%
  filter(term != "(Intercept)", losses_group != "0") %>%
  mutate(sig = ifelse(p.value >= 0.05, 'n.s.', ifelse(p.value < 0.05 & estimate > 0, '+', '-'))) %>%
  ggplot(aes(x = losses_group, y = term, colour = p.value < 0.05)) +
  geom_point(size = 4) +
  theme_large_classic() +
  theme(legend.position = "bottom") +
  labs(y = "", x = "Number of chromosomes with LOH") +
  ggtitle('(f)')


arm_cors <- joint_arm_losses %>%
  left_join(degree_aneu %>% distinct(id, n_losses, losses_group), by = c("id")) %>%
  filter(!is.na(losses_group)) %>%
  group_by(losses_group, chr_arm) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  mutate(arm = str_extract(chr_arm, ".{1}$")) %>%
  mutate(chr = sub(".{1}$", "", chr_arm)) %>%
  ungroup() %>%
  select(-chr_arm) %>%
  pivot_wider(names_from = arm, values_from = loss_rate) %>%
  filter(!is.na(p)) %>%
  ggplot(aes(x = p, y = q, colour = as.factor(losses_group), group = as.factor(losses_group))) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor(size = 4) +
  geom_label_repel(aes(label = chr)) +
  theme_large() +
  theme(legend.position = "right") +
  labs(colour = "# LOH chromosomes", y = "% of cases losing q-arm", x = "% of cases losing p-arm") + 
  ggtitle('(b)')

aneuploid_v_hypo <- joint_arm_losses %>%
  filter(proj %in% c(enough_lh, "ALL")) %>%
  group_by(Class, proj, chr_arm) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  filter(Class %in% c("Aneuploid", "Low-Hypodiploid")) %>%
  pivot_wider(names_from = Class, values_from = loss_rate) %>%
  clean_names() %>%
  mutate(grp = low_hypodiploid > 0.5) %>%
  ggplot(aes(y = aneuploid, x = low_hypodiploid, group = grp)) +
  geom_point() +
  geom_text_repel(aes(label = paste(chr_arm, proj))) +
  geom_smooth(method = "lm") +
  stat_cor(size = 8)

suppfig2 <- ( nh_all | arm_cors) / (jacc_tile) / (location_no17 | lr_tissues_figure | loss_group_figure)
ggsave(plot = suppfig2, file = 'paper/supp_fig2.png', width = 27, height = 30)


## Stats in text ##

# non-uniformity
joint_arm_losses %>%
  filter(Class == "Low-Hypodiploid", proj %in% enough_lh) %>%
  group_by(proj, chr_arm) %>%
  summarize(n_losses = sum(loss == "Lost"), prop_lost = n_losses / n()) %>%
  group_by(proj) %>%
  summarize(tidy(chisq.test(n_losses))) %>%
  arrange(desc(p.value)) # all p < 1e-6

#non uniformity for whole chromosomes
joint_chr_losses %>%
  filter(Class == "Low-Hypodiploid", proj %in% c(enough_lh, 'ALL')) %>%
  group_by(proj, chr) %>%
  summarize(n_losses = sum(loss == "Lost"), prop_lost = n_losses / n()) %>%
  group_by(proj) %>%
  summarize(sha = Entropy(prop_lost), gi = Gini(prop_lost), tidy(chisq.test(n_losses))) %>%
  arrange(desc(statistic)) # all p < 1e-4 

#uniformity for NH-ALL
joint_chr_losses %>%
  filter(Class == "Near-Haploid", proj == 'ALL') %>%
  group_by(proj, chr) %>%
  summarize(n_losses = sum(loss == "Lost"), prop_lost = n_losses / n()) %>%
  ungroup() %>%
  summarize(tidy(chisq.test(n_losses))) 

# 77% of retained chromosomes in ALL are 21, 14 and 18
joint_chr_losses %>%
  filter(Class == "Near-Haploid", proj == "ALL") %>%
  group_by(chr) %>%
  summarize(n_losses = sum(loss == "Lost"), loss_rate = mean(loss == "Lost"), n_retained = sum(loss == "Retained")) %>%
  arrange(loss_rate) %>%
  mutate(prop_of_all_retained = n_retained / sum(n_retained)) %>%
  head(3) %>%
  ungroup() %>%
  mutate(sum(prop_of_all_retained))

# chr7 is retained in >90% of LH cases in 10/16 cancer types
joint_chr_losses %>%
  filter(proj %in% c(enough_lh, "ALL")) %>%
  group_by(Class, proj, chr) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  filter(chr == "chr7", Class == "Low-Hypodiploid") %>%
  ungroup() %>%
  count(loss_rate < 0.1)

# most retained chromosomes overall -- note no cancer type filtering
joint_chr_losses %>%
  group_by(Class, chr) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  filter(Class == "Low-Hypodiploid") %>%
  ungroup() %>%
  arrange(loss_rate)

## univariate LR: length alone
degree_aneu %>%
  group_by(chr, losses_group, sum_p_haplo, location, chrom_tsg_og_score, length) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  group_by(losses_group) %>%
  do(tidy(glm(loss_rate ~ length, data = .))) %>% filter(term == 'length')

## univariate LR: location alone
degree_aneu %>%
  group_by(chr, losses_group, sum_p_haplo, location, chrom_tsg_og_score, length) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  group_by(losses_group) %>%
  do(tidy(glm(loss_rate ~ location, data = .))) %>% filter(term == 'locationPeripheral')

## univariate LR: 
degree_aneu %>%
  group_by(chr, losses_group, sum_p_haplo, location, chrom_tsg_og_score, length) %>%
  summarize(loss_rate = mean(loss == "Lost")) %>%
  group_by(losses_group) %>%
  do(tidy(glm(loss_rate ~ sum_p_haplo, data = .))) %>% filter(term == 'sum_p_haplo')

#chr7
joint_chr_losses %>% filter(Class == 'Low-Hypodiploid', proj %in% c(enough_lh, 'ALL')) %>% group_by(chr, proj) %>% summarize(loss_rate = mean(loss == 'Lost')) %>% filter(chr == 'chr7')
joint_chr_losses %>% filter(Class %in% c('Low-Hypodiploid', 'Near-Haploid'), proj != 'ALL') %>% group_by(chr) %>% summarize(loss_rate = mean(loss == 'Lost')) %>% arrange(loss_rate)

#median loss rates by location 
hypo_chr_loss_rates %>% group_by(location) %>% summarize(median(loss_rate))