# Transcriptomics Analysis on the role of NF2 in early cardiac mesodermal lineage specification

This is repository contains R scripts for the analysis of NF2 WT and KO transcriptomics across cardiomyocyte differentiation time-course (Day 0, 1, 5 and 14).

Tradeseq were used for Trajectory inference analysis of bulk RNA-sequencing data, and the TRIAGE tool were used for the prioritization of transcription regulatory factors from differentially expressed genes. 

Software References:

[1] Van den Berge K, Roux de Bezieux H, Street K, Saelens W, Cannoodt R, Saeys Y, et al. Trajectory-based differential expression analysis for single-cell sequencing data. Nat Commun. 2020;11(1):1201.

[2] Shim WJ, Sinniah E, Xu J, Vitrinel B, Alexanian M, Andreoletti G, et al. Conserved Epigenetic Regulatory Logic Infers Genes Governing Cell Identity. Cell Syst. 2020;11(6):625-39 e13.

[3] Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010;26(1):139-40.






## Local R implementation

Clone the project

```bash
  git clone https://github.com/Foo-Lab/NF2_TranscriptomicsAnalysis.git
```

Go to the project directory or access using R-studio (Set working directory) - with R-3.5.1+

RScripts can be found in the source folder.


### R Scripts 
Below are the following scripts for the respective analysis in the manuscript:
```
CRISPR_Scatterplot.r - Plotting and labeling of scatterplots from CRISPR Screen data
CRISPR_Tileplot.r - Plotting of tileplots from CRISPR Screen data
DifferentialAnalysis_edgeR.r - Differential gene expression analysis 
TRIAGE_Visualization.R - Visualization of TRIAGE data
TrajectoryInference.r - Visualization of bulk pseudotime data
```

### Session 
```
sessionInfo()

```

```
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] LC_COLLATE=English_Singapore.1252  LC_CTYPE=English_Singapore.1252    LC_MONETARY=English_Singapore.1252
[4] LC_NUMERIC=C                       LC_TIME=English_Singapore.1252    

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tradeSeq_1.12.0        ggpointdensity_0.1.0   viridis_0.6.4          viridisLite_0.4.2      dplyr_1.1.2           
 [6] stringr_1.5.0          reshape_0.8.9          scales_1.2.1           goseq_1.50.0           geneLenDataBase_1.34.0
[11] BiasedUrn_2.0.11       RColorBrewer_1.1-3     VennDiagram_1.7.3      futile.logger_1.4.3    gplots_3.1.3          
[16] ggplot2_3.4.3          edgeR_3.40.2           limma_3.54.2          

loaded via a namespace (and not attached):
 [1] nlme_3.1-160                bitops_1.0-7                matrixStats_0.63.0          bit64_4.0.5                
 [5] filelock_1.0.2              progress_1.2.2              httr_1.4.7                  GenomeInfoDb_1.34.9        
 [9] tools_4.2.2                 utf8_1.2.3                  R6_2.5.1                    KernSmooth_2.23-20         
[13] mgcv_1.8-41                 DBI_1.1.3                   BiocGenerics_0.44.0         colorspace_2.1-0           
[17] withr_2.5.2                 gridExtra_2.3               tidyselect_1.2.0            prettyunits_1.2.0          
[21] bit_4.0.5                   curl_5.0.0                  compiler_4.2.2              cli_3.6.1                  
[25] Biobase_2.58.0              formatR_1.14                xml2_1.3.4                  DelayedArray_0.23.2        
[29] rtracklayer_1.58.0          caTools_1.18.2              pbapply_1.7-2               slingshot_2.6.0            
[33] rappdirs_0.3.3              digest_0.6.31               Rsamtools_2.14.0            XVector_0.38.0             
[37] pkgconfig_2.0.3             MatrixGenerics_1.10.0       dbplyr_2.4.0                fastmap_1.1.1              
[41] rlang_1.1.1                 rstudioapi_0.15.0           RSQLite_2.3.1               BiocIO_1.8.0               
[45] generics_0.1.3              BiocParallel_1.32.6         gtools_3.9.4                RCurl_1.98-1.12            
[49] magrittr_2.0.3              GO.db_3.16.0                GenomeInfoDbData_1.2.9      Matrix_1.6-0               
[53] Rcpp_1.0.10                 munsell_0.5.0               S4Vectors_0.36.2            fansi_1.0.4                
[57] TrajectoryUtils_1.6.0       lifecycle_1.0.4             stringi_1.7.12              yaml_2.3.7                 
[61] SummarizedExperiment_1.28.0 zlibbioc_1.44.0             plyr_1.8.8                  BiocFileCache_2.6.1        
[65] blob_1.2.4                  parallel_4.2.2              crayon_1.5.2                lattice_0.20-45            
[69] splines_4.2.2               Biostrings_2.66.0           GenomicFeatures_1.50.4      hms_1.1.3                  
[73] KEGGREST_1.38.0             locfit_1.5-9.7              pillar_1.9.0                igraph_1.5.1               
[77] GenomicRanges_1.50.2        rjson_0.2.21                codetools_0.2-18            biomaRt_2.54.1             
[81] stats4_4.2.2                futile.options_1.0.1        XML_3.99-0.14               glue_1.6.2                 
[85] lambda.r_1.2.4              png_0.1-8                   vctrs_0.6.4                 gtable_0.3.4               
[89] cachem_1.0.8                princurve_2.1.6             restfulr_0.0.15             SingleCellExperiment_1.20.1
[93] tibble_3.2.1                GenomicAlignments_1.34.1    AnnotationDbi_1.60.2        memoise_2.0.1              
[97] IRanges_2.32.0

```


## Contact
Please submit a github issue or contact me @ mick.lee@u.nus.edu if you experience any issues/bugs with R scripts.




