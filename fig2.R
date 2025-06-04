source('setup.R')

tcga_classes <- fread('TCGA_ploidy_classes.tsv')
tcga_losses <- fread('tcga_losses.tsv')
mitcn_meta <- fread('mitcn_meta.tsv')
mitcn_prep <- fread('mitcn_prep.tsv') 

mit_arms <- read.table('mitcn_arms.bed', header = F) %>%
  mutate(len = V3 - V2) %>% 
  rename(chr = V1, start = V2, end = V3, id = V4, CN = V5, arm = V9) %>%
  mutate(chr_arm = paste0(sub('chr', '', chr), arm)) %>% #has no X or Y
  left_join(mitcn_meta, by = c('id'))

fa <- fread('fasc_arms.bed') %>% 
  mutate(len = V3 - V2) %>% 
  rename(chr = V1, start = V2, end = V3, GDC_Aliquot = V4, CN = V5, mCN = V6, proj = V7, arm = V11) %>%
  mutate(loh_indicator = ifelse(mCN == 0, 1, 0)) %>% 
  mutate(chr_arm = paste0(sub('chr', '', chr), arm)) 

enough_lh <- tcga_classes %>% filter(group == 'Low-Hypodiploid') %>% count(proj) %>% filter(n >= 15) %>% pull(proj) 

print('Making Figure 2')

mit_arm_losses <- mit_arms %>%
  group_by(id, case_id, Sex, Subclone, n_chr_nosex, group, any_lows, any_hyper, n_clones, chr, chr_arm, ploidy) %>%  
  summarize(arm_somy = round(sum(CN*len)/sum(len),0), loss_perc = sum(len[CN < 2])/sum(len)) %>%
  mutate(loss = ifelse(loss_perc > 0.9, 'Lost', 'Retained')) %>%
  mutate(Class = group, proj = 'ALL') %>%
  ungroup() 

mit_chr_losses <- mitcn_prep %>%
  group_by(id, case_id, Sex, Subclone, n_chr_nosex, group, any_lows, any_hyper, n_clones, chr, ploidy) %>%  
  summarize(chr_somy = round(sum(CN*len)/sum(len),0), loss_perc = sum(len[CN < 2])/sum(len)) %>%
  mutate(loss = ifelse(loss_perc > 0.9, 'Lost', 'Retained')) %>%
  mutate(Class = group, proj = 'ALL') %>%
  ungroup() 

mit_chr_losses_tojoin <- mit_chr_losses %>% 
  mutate(Class = ifelse(ploidy > 2.7, 'WGD-high', Class)) %>% #ugh this is a fucking mess
  select(id, Class, chr, proj, loss, chr_somy) 

mit_arm_losses_tojoin <- mit_arm_losses %>% 
  mutate(Class = ifelse(ploidy > 2.7, 'WGD-high', Class)) %>% #ugh this is a fucking mess
  select(id, Class, chr, proj, loss, chr_arm) 

tcga_arm_losses <- fa %>% 
  mutate(len = (end - start) + 1) %>%
  mutate(loh_indicator = ifelse(mCN == 0, 1, 0)) %>%
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
  ungroup()

joint_arm_losses <- tcga_arm_losses %>%
  left_join(tcga_classes, by = c("GDC_Aliquot", "proj")) %>%
  select(GDC_Aliquot, Class, chr, chr_arm, proj, loss) %>%
  rename(id = GDC_Aliquot) %>%
  rbind(mit_arm_losses_tojoin) %>%
  ungroup()

chr_loss_rates_pan_lowhypo <- joint_chr_losses %>%
  filter(Class == 'Low-Hypodiploid') %>%
  group_by(chr) %>%
  summarize(loss_rate = mean(loss == 'Lost'))

arm_loss_rates_pan_lowhypo <- joint_arm_losses %>%
  filter(Class == 'Low-Hypodiploid') %>%
  group_by(chr_arm) %>%
  summarize(loss_rate = mean(loss == 'Lost'))

#Loss Rates graph
joint_loss_rates_inequality <- joint_arm_losses %>%
  filter(Class == 'Low-Hypodiploid') %>%
  filter(proj %in% c(enough_lh, 'ALL')) %>%
  filter(!chr_arm %in% acrocentric_arms) %>%
  group_by(chr_arm, proj) %>% 
  summarize(loss_rate = mean(loss == 'Lost')) %>%
  group_by(proj) %>% 
  mutate(sha = Entropy(loss_rate), gi = Gini(loss_rate)) 

joint_loss_rates_graph <- joint_loss_rates_inequality %>%
  ggplot(aes(x = reorder(proj, gi), y = factor(chr_arm, levels = chr_levels), fill = loss_rate)) + 
  geom_tile() + scale_fill_viridis_c() + 
  xlab('') + ylab('') + theme_large(base_size = 14) + theme(legend.position = 'right') + labs(fill = 'Lost') + ggtitle('(a)') + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# arm cors?

#similarity to dipliod/aneuploid??

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

htmap <- ComplexHeatmap::pheatmap(toheat, col = rev(viridis::viridis(15, option = "D")), show_rownames = F, show_colnames = F, annotation_row = idg2 %>% select(group), annotation_names_row = F, heatmap_legend_param = list(title = "Distance"), annotation_colors = list(group = deframe(distinct(idg2)))) %>% as.ggplot()

## jaccard tiles
jacc_tile <- toheat %>%
  as.data.frame() %>%
  rownames_to_column(var = 'Sample1') %>%
  pivot_longer(-Sample1, names_to = 'Sample2', values_to = 'Distance', values_drop_na = T) %>%
  left_join(idg, by = c('Sample1' = 'id')) %>%
  left_join(idg, by = c('Sample2' = 'id'), suffix = c('1', '2')) %>%
  filter(Sample1 != Sample2) %>% #remove self-comparisons
  mutate(Similarity = 1 - Distance) %>% mutate(comparison = paste(group1, group2, sep = '_'), same = group1 == group2) %>% group_by(group1, group2, same) %>% summarize(medsim = median(Similarity), std = sd(Similarity)) %>% ggplot(aes(x = group1, y = group2, fill = medsim)) + geom_tile() + geom_text(aes(label = round(medsim, 2))) + scale_fill_continuous(type = 'viridis') + xlab('') + ylab('') + labs(fill = 'Median\nSimilarity') + theme_large() + theme(legend.position = 'right')

## Nuclear Location
location_palette <- c(
  "Central" = "#1b9e77",  # Teal green
  "Peripheral"       = "#3288bd",  # Rich blue
  "Other"           = "#a6bddb"   # Light blue-gray (fallback group)
)

set.seed(123)
jitter_width <- 0.4
jittered_data <- chr_loss_rates_pan_lowhypo %>%
  mutate(location = case_when(
    chr %in% c("chr1", "chr8", "chr11", "chr17", "chr19") ~ "Central",
    chr %in% c("chr13", "chr15", "chr18", "chr21", "chrX") ~ "Peripheral",
    TRUE ~ "Other"
  )) %>%
  mutate(jitter_x = as.numeric(factor(location)) + runif(n(), -jitter_width, jitter_width))

nuclear_location <- ggplot(jittered_data, aes(x = location, y = loss_rate)) +
  geom_boxplot(aes(fill = location), alpha = 0.2) +
  geom_point(aes(x = jitter_x, colour = location), size = 5) +
  geom_text_repel(size = 6, 
                  aes(x = jitter_x, label = sub('chr', '', chr)),
                  nudge_y = 0.02,             
                  box.padding = 0.4,         
                  point.padding = 0.5,        
                  max.overlaps = Inf         
  ) +
  theme_large() +
  theme(legend.position = "none") +
  stat_compare_means(hjust = -0.7) +
  xlab('Nuclear Location') + ylab('Chromosome Loss Rate') 
nuclear_location


### length ###
chr_length <- arms %>% group_by(chr) %>% summarize(first = min(st), last = max(en))

length_cor <- chr_loss_rates_pan_lowhypo %>%
  left_join(chr_length, by = c('chr')) %>%
  mutate(length = last/1000000) %>%
  ggplot(aes(x = round(length, 2), y = loss_rate)) +
  geom_point(size = 5, colour = 'darkblue') + geom_smooth(method = 'lm', fill = 'lightsteelblue', colour = 'black') + stat_cor(vjust = 1.5) + 
  geom_text_repel(size = 6, aes(label = sub('chr', '', chr))) +
  theme_large() + xlab('Chromosome Length (Mb)') + ylab('Chromosome Loss Rate') 

tuson <- fread('~/Downloads/tuson_chrom.csv', skip = 45) %>% clean_names()

## tuson ##
tuson_loss <- chr_loss_rates_pan_lowhypo %>%
  mutate(match_chr = as.numeric(sub("chr", "", chr))) %>%
  left_join(tuson, by = c("match_chr" = "chromosome")) %>%
  ggplot(aes(x = chrom_tsg_og_score, y = loss_rate)) +
  geom_point(size = 5, colour = 'darkblue') +
  geom_smooth(method = 'lm', fill = 'lightsteelblue', colour = 'black') +
  stat_cor() +
  theme_large() +
  xlab("CHROM score (TSG - OG)") +
  ylab("Chromosome Loss Rate") +
  geom_text_repel(size = 6, aes(label = sub('chr', '', chr))) 


## DS ##
ds <- fread('~/Downloads/Collins_rCNV_2022.dosage_sensitivity_scores.tsv') %>% clean_names() %>% rename(gene = number_gene)
ds_rates <- ga %>% left_join(ds, by = c('V4' = 'gene')) %>% filter(!is.na(p_haplo)) %>% group_by(V1) %>% summarize(sum_p_haplo = sum(p_haplo), mean_p_haplo = mean(p_haplo)) 

ds_loss <- chr_loss_rates_pan_lowhypo %>% 
  left_join(ds_rates, by = c("chr" = "V1")) %>%
  ggplot(aes(x = sum_p_haplo, y = loss_rate)) +
  geom_point(size = 5, colour = 'darkblue') +
  geom_smooth(method = 'lm', fill = 'lightsteelblue', colour = 'black') +
  stat_cor(label.x = Inf, label.y = Inf, hjust = 1.1, vjust = 1.5) +  theme_large() +
  xlab("Haplosensitivity score") +
  ylab("Chromosome Loss Rate") +
  geom_text_repel(size = 6, aes(label = sub('chr', '', chr))) 

(joint_loss_rates_graph | htmap) / (length_cor | nuclear_location ) / ( tuson_loss | ds_loss)
#(joint_loss_rates_graph | htmap) /  / ( tuson_loss | nuclear_location)

ggsave('paper/fig2_paper.png', width = 25, height = 30)