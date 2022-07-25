# Description of code for manuscript "A Mendelian randomization-based exploration of red blood cell distribution width and mean corpuscular volume with risk of hemorrhagic strokes"
Authors: Jundong Liu,<sup>1</sup> Elizabeth L Chou,<sup>2</sup> Kui Kai Lau,<sup>34</sup> Peter YM Woo,<sup>5</sup> Tsz Kin Wan,<sup>6</sup> Ruixuan Huang,<sup>6</sup> Kei Hang Katie Chan<sup>167*</sup>

####

1.	Department of Biomedical Sciences, City University of Hong Kong, Hong Kong SAR, China
2.	Massachusetts General Hospital, Boston, MA, USA 
3.	Division of Neurology, Department of Medicine, The University of Hong Kong, Hong Kong SAR, China 
4.	The State Key Laboratory of Brain and Cognitive Sciences, The University of Hong Kong, Hong Kong SAR, China 
5.	Department of Neurosurgery, Kwong Wah Hospital, Hong Kong SAR, China
6.	Department of Electrical Engineering, City University of Hong Kong, Hong Kong SAR, China
7.	Department of Epidemiology, Centre for Global Cardiometabolic Health, Brown University, RI, USA 

​	*Corresponding author

####

## Requirements
This codes are test on Rstudio(2021.09.1+372 "Ghost Orchid" Release (8b9ced188245155642d024aa3630363df611088a, 2021-11-08)), requiring following packages.

Rstudio: TwoSampleMR, readxl, tidyverse, RadialMR, boot, LDlinkR,  dplyr, MVMR, MRPRESSO, boot.pval, MendelianRandomization,  GWAS.utils, phenoscanner, cause, RMVMR, RColorBrewer,  cowplot, gridExtra, MASS, glmnet, quantreg, robustbase,  Hmisc, forestplot, qvalue, fdrtool, rlist, pipeR, vcfR, stringr, data.table

# Data source

1. Medical Research Council Integrative Epidemiology Unit OpenGWAS, https://gwas.mrcieu.ac.uk/
2. UK Biobank database, https://www.ukbiobank.ac.uk/
3. UK Biobank - Neale lab, http://www.nealelab.is/uk-biobank
4. FINNGEN, https://www.finngen.fi/fi

# Steps of analysis
1. Genetic correlation: preporcessing_data_for_phenoSpD.R
1. UMR: preporcessing_data_for_umr.R, umr_estimates.R, umr_tool.R
2. MVMR: mvmr_estimates.R, mvmr_tool.R
3. Coloc: read_vcf.R, coloc.R
3. Mediation analysis: mediation analysis.R

