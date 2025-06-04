set -euo pipefail
set -x

Rscript run_cytoconverter_mitALL.R

## process cytoconverter results to get Mitelman-ALL copy numbers
Rscript prep_cytoconverter_output.R
rm all_mit/*.bed
bash handle_mit_diploids.sh mit_diploids.tsv all_mit #create a pseudo CN change for the diploid cases
bash convert_cyto.sh all_cytoconverted.bed all_mit all_mitelman_CNs.bed 
Rscript process_mitelman_cns.R #this outputs mitcn_all_for_arms.bed, mitcn_prep.tsv, mitcn_meta.tsv and mitALL_preMEDICC.bed 

## Run MEDICC on multi-clone ALL samples with at least one LH or NH clone to make WGD calls
bash makemedicc.sh mitALL_preMEDICC.bed mitALL_for_MEDICC mitALL combined_mitALL_forMEDICC.txt
source ~/anaconda3/etc/profile.d/conda.sh
conda activate medicc_env
bash runmedicc.sh combined_mitALL_forMEDICC.txt MITMED > runmedicc.log 2>&1
conda deactivate

## process TCGA copy number data (ASCAT seg files)
Rscript filter_TCGA_ASCAT.R 


## MH score heuristic
Rscript MH_score_heuristic.R

## do analyses and make figures
Rscript fig1.R 

bedtools intersect -a filtered_ascat.bed -b hg38_arms.bed  -wb > fasc_arms.bed  
bedtools intersect -a mitcn_all_for_arms.bed -b hg38_arms.bed  -wb > mitcn_arms.bed
Rscript fig2.R

Rscript fig3.R
Rscript fig4.R
