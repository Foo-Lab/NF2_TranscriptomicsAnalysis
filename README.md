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

## Contact
Please submit a github issue or contact me @ mick.lee@u.nus.edu if you experience any issues/bugs with R scripts.




