# Heiland-et-al-2022
R code and accompanying files to reproduce miRNA target investigation and pathway enrichment analysis as performed for Heiland et al 2022

File descriptions:

"mirna335-r-workflow-20220622.R" 
R code to reproduce miRNA target investigation and pathway enrichment analysis. Recreates Supplementary Table 1, Supplementary Table 2, and Figure 2A

"f_validation_db_200527.txt" 

Experimentally validated miRNA-target interactions collated from miRTarbase V7 (Chou et al., 2018) and TarBase V8 (Karagkouni et al., 2018) as specified in the paper. Required around lines 150-165 in R code. 
**File is 94MB, and cannot be uploaded to Github. Please email niamhmconnolly@rcsi.com to request.**

"ALL_EPILEPSY_GENES.csv" 

Genes/proteins implicated in epilepsy - in-house database collated from various sources as specified in the paper. Required around line 375-385 in R code.

"mrna-brain-expressed.txt" 

Genes/proteins expressed in the brain - downloaded from Protein Atlas, as specified in the paper. Required around line 385-400 in R code.
