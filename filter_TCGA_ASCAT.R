source('setup.R')

print('Filtering TCGA')

#### Import GDC metadata ####

meta <- read.table('../cancer_wgd/GDCquery_results_allcancers.txt', header = TRUE)
ffpe <- read.table('../cancer_wgd/ffpe.txt', header = TRUE, sep = '\t') #downloaded using GDCquery from TCGAbiolinks, as described in cancer_wgd/GDCqueries.R

#### Recommended exclusions by ASCAT team ####

ascat_contam_swap <- "TCGA-04-1351, TCGA-04-1371, TCGA-06-0178, TCGA-57-1586, TCGA-V5-AASX-b" 
ascat_germline <- "TCGA-05-4403, TCGA-06-2559, TCGA-12-3651, TCGA-25-1877, TCGA-26-5135, TCGA-2G-AAEW, TCGA-77-6842-a, TCGA-86-8054, TCGA-91-A4BD, TCGA-97-8179, TCGA-A5-A2K4, TCGA-A8-A07L, TCGA-AB-2812, TCGA-AB-2843, TCGA-AB-2847, TCGA-AB-2855, TCGA-AB-2856, TCGA-AB-2887, TCGA-AB-2891, TCGA-AB-2909, TCGA-AB-2913, TCGA-AB-2924, TCGA-AB-2944, TCGA-AB-2954, TCGA-AB-2997, TCGA-AB-3008, TCGA-AC-A2QJ, TCGA-AH-6547, TCGA-AP-A0LQ, TCGA-AY-4071, TCGA-BH-A1EW, TCGA-BJ-A28W-a, TCGA-BK-A139-e, TCGA-BM-6198, TCGA-BT-A20U, TCGA-C5-A1MK, TCGA-CG-5717, TCGA-CH-5753, TCGA-CH-5769, TCGA-CV-6441, TCGA-D8-A146, TCGA-D8-A1XC, TCGA-DD-AAEH, TCGA-DJ-A2Q2, TCGA-DM-A1D6, TCGA-E2-A15G, TCGA-EE-A29D, TCGA-EY-A1G8, TCGA-FD-A62N, TCGA-FS-A4F0, TCGA-GS-A9U4, TCGA-P5-A72U, TCGA-UZ-A9PJ" 
ascat_exclus <- "TCGA-05-4403, TCGA-06-2559, TCGA-12-3651, TCGA-25-1877, TCGA-26-5135, TCGA-2G-AAEW, TCGA-77-6842-a, TCGA-86-8054, TCGA-91-A4BD, TCGA-97-8179, TCGA-A5-A2K4, TCGA-A8-A07L, TCGA-AB-2812, TCGA-AB-2843, TCGA-AB-2847, TCGA-AB-2855, TCGA-AB-2856, TCGA-AB-2887, TCGA-AB-2891, TCGA-AB-2909, TCGA-AB-2913, TCGA-AB-2924, TCGA-AB-2944, TCGA-AB-2954, TCGA-AB-2997, TCGA-AB-3008, TCGA-AC-A2QJ, TCGA-AH-6547, TCGA-AP-A0LQ, TCGA-AY-4071, TCGA-BH-A1EW, TCGA-BJ-A28W-a, TCGA-BK-A139-e, TCGA-BM-6198, TCGA-BT-A20U, TCGA-C5-A1MK, TCGA-CG-5717, TCGA-CH-5753, TCGA-CH-5769, TCGA-CV-6441, TCGA-D8-A146, TCGA-D8-A1XC, TCGA-DD-AAEH, TCGA-DJ-A2Q2, TCGA-DM-A1D6, TCGA-E2-A15G, TCGA-EE-A29D, TCGA-EY-A1G8, TCGA-FD-A62N, TCGA-FS-A4F0, TCGA-GS-A9U4, TCGA-P5-A72U, TCGA-UZ-A9PJ, TCGA-04-1351, TCGA-04-1371, TCGA-06-0178, TCGA-57-1586, TCGA-V5-AASX-b" 
ascat_exclus <- gsub(',', '', trimws(scan(text = ascat_exclus, what = ','))) 

#### Import and format combined ASCAT segment data ####

ascat <- fread('~/cat_all_AS.txt', header = TRUE) %>%
  filter(GDC_Aliquot != 'GDC_Aliquot') #remove interspersed headers from cat

#convert UUIDs to barcodes using TCGAutils to query GDC
aliquots <- unique(ascat$GDC_Aliquot) #11104 rows
codes <- UUIDtoBarcode(aliquots, from_type = 'aliquot_ids') #11104 rows

codes <- codes %>%
  separate(col = 'portions.analytes.aliquots.submitter_id', into = c('TCGA', 'TSS', 'Individual', 'Sample', 'Portion', 'Plate', 'Center'), sep = '-', remove = FALSE) %>%
  mutate(Patient = paste(TCGA, TSS, Individual, sep = '-')) %>%
  mutate(Specimen = paste(Patient, Sample, sep = '-')) %>%
  separate(col = Sample, into = c('Sample.Type', 'Vial'), sep = 2) #11104 rows

codes %>% fwrite('TCGA_codes.tsv', quote = F, sep = '\t', col.names = T, row.names = F)

#### Filter ASCAT samples ####

filtered_codes <- codes %>%
  filter(Sample.Type %in% c('01', '03', '09')) %>% #remove metastatic and recurrent; keep solid primary tumour and primary blood cancers
  filter(!Patient %in% ascat_exclus) %>% 
  filter(Specimen %in% ffpe[ffpe$is_ffpe == FALSE,]$submitter_id) %>%  #remove FFPE samples
  arrange(Vial, desc(Plate), Portion) %>%  #Plate sorted as suggested by Broad Institute
  left_join(meta, by = c('Patient')) %>%
  distinct(Patient, .keep_all = TRUE)  

filtered_ascat <- ascat %>% #removing sex chromosomes here
  filter(GDC_Aliquot %in% filtered_codes$portions.analytes.aliquots.aliquot_id) %>% #keep only preferred samples
  mutate_at(c('Start', 'End', 'Copy_Number', 'Major_Copy_Number', 'Minor_Copy_Number'), as.numeric)  %>% #convert these columns to numeric
  left_join(filtered_codes, by = c('GDC_Aliquot' = 'portions.analytes.aliquots.aliquot_id')) %>% #get barcodes 
  separate(project, into = c('tcga', 'proj'), sep = '-', remove = F) %>%
  mutate_at(c('Start', 'End', 'Minor_Copy_Number', 'Major_Copy_Number', 'Copy_Number'), as.numeric) %>%
  select(Chromosome, Start, End, GDC_Aliquot, Copy_Number, Minor_Copy_Number, proj) %>%
  filter(!Chromosome %in% c('chrX', 'chrY'))

filtered_ascat %>% fwrite('filtered_ascat.tsv', sep = '\t', quote = F, col.names = T, row.names = F)

filtered_ascat %>% 
  mutate(Start = Start - 1) %>%
  arrange(Chromosome, Start, End) %>%
  fwrite('filtered_ascat.bed', sep = '\t', quote = F, col.names = F, row.names = F)

tcga_losses <- filtered_ascat %>% 
  mutate(len = (End - Start) + 1) %>%
  rename(CN = Copy_Number, mCN = Minor_Copy_Number, chr = Chromosome, start = Start, end = End) %>%
  mutate(loh_indicator = ifelse(mCN == 0, 1, 0)) %>%
  mutate(prod = CN*len, loh_len = loh_indicator*len) %>%
  group_by(GDC_Aliquot, proj) %>%
  mutate(ploidy = sum(prod)/sum(len)) %>% #calculate mean copy number per sample 
  group_by(GDC_Aliquot, chr, proj, ploidy) %>%
  summarize(loh = sum(loh_len)/sum(len), chr_somy = round(sum(prod)/sum(len), 0)) %>% #calculate proportion of LOH per chromosome per sample
  mutate(loss = ifelse(loh > 0.9, 'Lost', 'Retained')) %>% #this should also cover nullisomy
  mutate(min_somy = ifelse(loss == 'Retained', 2, pmin(1, chr_somy))) %>% #base_somy = 2 because excluding the sex chromosomes
  ungroup()


tcga_classes <- tcga_losses %>% 
  group_by(GDC_Aliquot, proj, ploidy) %>% #regroup by sample
  summarize(n_loh_nosex = sum(loss == 'Lost'), n_chr_nosex = sum(chr_somy), min_chr_nosex = sum(min_somy), n_disomic_nosex = sum(chr_somy == 2), n_nullisomic_nosex = sum(chr_somy == 0)) %>%
  mutate(group = ifelse(min_chr_nosex < 28, 'Near-Haploid', ifelse(min_chr_nosex < hypo_threshold, 'Low-Hypodiploid', 'Other'))) %>% #LH should have LOH of >= 6 chromosomes assuming no nullisomies
  left_join(ascat_medicc[, c('sample', 'WGD_status')], by = c('GDC_Aliquot' = 'sample')) %>% 
  mutate(wgd = WGD_status) %>%
  mutate(Class = ifelse(group != 'Other', group, 
                        ifelse(wgd == 'WGD', 'Polyploid', 
                               ifelse(n_disomic_nosex == 22, 'Diploid', 'Aneuploid')))) %>% #create 5 ploidy classes
  ungroup() %>%
  filter(n_chr_nosex >= 22, min_chr_nosex >= 22, n_nullisomic_nosex == 0) %>% #there were 2 samples with min_chr_nosex < 22 and 12 samples with a nullisomic autosome
  left_join(codes[, c('portions.analytes.aliquots.aliquot_id', 'Patient')], by = c('GDC_Aliquot' = 'portions.analytes.aliquots.aliquot_id')) #all accounted for
  

tcga_classes %>% write.table('TCGA_ploidy_classes.tsv', sep = '\t', col.names = T, row.names = F, quote = F)

tcga_losses %>%
  filter(GDC_Aliquot %in% tcga_classes$GDC_Aliquot) %>%
  write.table("tcga_losses.tsv", sep = "\t", quote = F, col.names = T, row.names = F)


tissue_rates <- tcga_classes %>%
  group_by(proj) %>%
  summarise(wgd_rate = sum(wgd == "WGD") / n(), 
            hypo_rate = mean(group == "Low-Hypodiploid"), 
            wgd_rate_nohypo = sum(wgd == "WGD" & group == "Other") / sum(group == "Other")) 

tissue_rates %>% fwrite('tissue_ploidy_rates.tsv', quote = F, sep = '\t', )
  
