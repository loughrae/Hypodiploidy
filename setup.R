setwd('~/PhD/Hypodiploidy')

## load packages
library(data.table)
library(tidyverse)

#getting-data packages
library(janitor)
library(TCGAbiolinks)
library(TCGAutils)
library(openxlsx)

#visualisation/analysis packages
library(ggpubr)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(ggplotify)
library(RColorBrewer)
library(survival)
library(survminer)
library(broom)
library(DescTools)
library(microViz)


select <- dplyr::select
rename <- dplyr::rename
options(scipen = 20)

## settings
hypo_threshold <- 39 #I use < so min_chr can be at most 38, i.e. has lost at least 6 chromosomes

## chromosome arm annotations
chr_levels <- paste0(rep(c(1:22, 'X', 'Y'), each = 2), c('p', 'q')) #for nice viz
acrocentric_arms <- paste0(c(13, 14, 15, 21, 22), 'p') #these arms don't have TCGA data
GRCh38.bands <- read.table('~/Downloads/GRCh38.bands', header = T)
arms <- GRCh38.bands %>% #reformat to get chromosome arm positions and lengths
  mutate(chr = paste0('chr', chr)) %>% 
  select(chr, start, end, name) %>% 
  mutate(arm = substr(name, start = 1, stop = 1)) %>% 
  group_by(chr, arm) %>% 
  summarize(st = min(start), en = max(end)) %>% 
  arrange(chr, st, en) %>% 
  select(chr, st, en, arm) 

## gene annotations 
ga <- read.table('~/Downloads/genes_arms.bed', sep = '\t') %>% add_count(V4) %>% filter(n == 1)

## MEDICC WGD calls based on ASCAT from GDCquery
ascat_medicc <- fread('ASCAT_MEDICC_calls.tsv', header = T)

## clinical data
clin <- read.table("~/Downloads/clin_query.txt", sep = "\t", header = T) %>% 
  select(submitter_id, age_at_index, gender, race, disease, vital_status, days_to_death, days_to_last_follow_up)

## themes 
theme_large_base <- function(base_theme, base_size = 18, base_family = "") {
  base_theme(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(size = base_size),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      strip.text = element_text(size = base_size)
    )
}

theme_large <- function(base_size = 18, base_family = "") {
  theme_large_base(theme_bw, base_size, base_family)
}

theme_large_classic <- function(base_size = 18, base_family = "") {
  theme_large_base(theme_classic, base_size, base_family)
}

class_palette <- c(
  "Near-Haploid" = "#1b9e77", 
  "Low-Hypodiploid"    = "#66c2a5", 
  "Polyploid"       = "#3288bd", 
  "Aneuploid"           = "#a6bddb" ,  
  "Diploid" = 'lightblue'
)

pairwise_comparisons <- list(
  c("Diploid", "Polyploid"),
  c("Diploid", "Low-Hypodiploid"),
  c("Polyploid", "Low-Hypodiploid")
)
