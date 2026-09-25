# Hypodiploidy
Scripts for paper on karyotypic evolution of hypodiploid tumours

**1. Processing Karyotypes **

All samples in the Mitelman Database on 25/09/2023 with morphology 'acute lymphoblastic leukaemia/lymphoblastic lymphoma' were downloaded via Google BigQuery. 

- 'run_cytoconverter_mitALL.R': long-form karyotypes are extracted and run through CytoConverter, which produces a table of copy number events (gains and losses). This takes ~30 minutes. - - 'prep_cytoconverter_output.R': filters out errors and samples with unknown chromosome numbers, and accounts for diploid samples with no copy number changes
- 'process_mitelman_cns.R': annotates clonal structures and computes statistics for each sample e.g losses
- 'make_medicc.sh' and 'runmedicc.sh': multi-clone samples with at least one LH/NH clone are run through MEDICC to make WGD calls

**2. Processing TCGA Data **

ASCAT copy number segment data was downloaded for all available TCGA cases from the GDC using TCGAbiolinks in R. 

- 'filter_TCGA_ASCAT.R': filters and deduplicates samples and computes loss statistics

**3. Figure 1: Prevalence **

- 'setup.R': is called by the figure scripts, sets up thresholds and figure defaults
- 'prevalence_fig.R': Figure 1: graphs and in-text statistics describing rates of hypodiploidy across cancers

**4. Figure 2: The Masked Hypodiploidy Score **

- 'MH_score_heuristic.R': sets up 'true positive' (MEDICC-called WGDs from multi-clone cases) and 'true negative' (genotyped hyperdiploids from Woodward et al 2023) cases and computes accuracy measures. Counts putative hidden hypodiploids in the Mitelman database.
- 'MH_score_fig.R': graphs and in-text statistics describing MH score distributions in hypodiploids vs hyperdiploids, and ranges of accuracy; processes TARGET-ALL dataset as another source of negatives (non-hypodiploid ALL cases)

**5. Figure 3: Patterns of Chromosome Loss **

- bedtools: maps copy number changes to chromosome arms
- 'patterns_fig.R': computes chromosome and chromosome arm loss rates and compares them to chromosome features (TSG/oncogene density, dosage sensitivity, length, nuclear location) - graphs and in-text statistics

**6. Figure 4: Generalised Chromosomal Instability **

- 'CIN_fig.R': compares hypodiploidy status to WGD rate, intra-tumour copy number heterogeneity (van Dijk et al), level of smaller-scale copy number changes, and survival; processes PCAWG dataset as an independent source of survival data - graphs and in-text statistics

**7. Figure 5: Origins of Chromosomal Instability **

- 'origins_fig.R': compares hypodiploidy status to mutational profiles (TP53 enrichment, MSI status, mutational signatures) and hypoxia (between ploidy classes and across cancer types) -- graphs and in-text statistics
  



