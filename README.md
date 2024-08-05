# Nanopore m6A analysis  
The codes were used for comprehensive studying of m6A modifications in muscle-invasive bladder cancer.  
## step1: Differential m6A modifications  
[DENA](https://github.com/weir12/DENA) was used for differential m6A modifications.  
The code was in STEP1_DENA.sh  
## step2: Jaspar_enrichment Analysis  
The figure 2 was used [Jasper](https://jaspar2020.genereg.net/enrichment/):  
> JASPAR_enrich.sh oneSetBg /data/hg38 /DENA/DENAListposttest_001.sort.bed /Draw/Anno/QiWeiWu_TTTNNNMeRatio_rm0_3col.sort.bed JASPAR_enrichment http://jaspar.genereg.net/api/v1/matrix/ 5  
## step3: Prognostic model Analysis  
The Prognostic model with m6A altered genes was contributed by STEP2_model.R  
