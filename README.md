# Enhancer Drug Targets
Identification of ovarian cancer subtype specific drug targets using epigenomic and transcriptomic data

![Pharmacogenomics Pipeline](https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true)

## 1 Create Enhancer based MAGIC matrix

### Scripts: 
1. CreateEnhancerMAGICMatrix.r
### Inputs:
1. GeneHancer gff file
2. TF ChIP-seq bed files
3. TF ChIP-seq meta data file
### Outputs:
1. MAGIC mtx matrix file
### Description

ChIP-seq TF peak files and meta data can be downloaded from ENCODE experiment matrix (https://www.encodeproject.org/matrix/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=TF+ChIP-seq&target.investigated_as=transcription+factor&assembly=GRCh38&files.file_type=bed+narrowPeak)

Once files are selected, use the Download button to download a text file containing the url of each peak file and a corresponding meta data file. Use curl to download bed.gz files as instructed.

To generate a 'MAGIC matrix' file based on GeneHancer known enhancer locations, use the r script 'CreateEnhancerMAGICMatrix.r'. This requires a copy of the geneHancer database in gff format. The 2017 version can be downloaded directly from GeneCards https://www.genecards.org/, and the latest version can be obtained by request.  

## 2 DifferentialAnalysis of H3K27ac Enhancer Sites

### Scripts: 
1. DifferentialAnalysisOfEnhancerSites.rmd
### Inputs:
1. GeneHancer gff file
2. BAM and corresponding MACS2 H3K27ac ChIP peak file from cancer cell lines
3. BAM and corresponding MACS2 H3K27ac ChIP peak file from cancer tissue
4. BAM and corresponding MACS2 H3K27ac ChIP peak file from healthy cell lines
### Outputs:
1. Upregulated enhancer site text file (containing two columns with list of background enhancers and upregulated enhancers)
### Description

## 3 Identify Enriched TF binding in differentially expressed H3K27ac enhancer sites

### Scripts: 
1. MAGIC_1_1.py
### Inputs:
1. MAGIC enhancer MTX file
2. Upregulated enhancer site text file
### Outputs:
1. Summary of upregulated TFs xls
2. Details of upregulated TFs xls
3. For further details see
### Description

## 4 Filter MAGIC results by enriched Motifs

### Scripts: 
1. MotifSearch.rmd
### Inputs:
1. 
### Outputs:
1. 
### Description

## 5 Identify complexes enriched TFs are members of and filter by expression data

### Scripts: 
1. ComplexComponentSearchAndFilter.rmd
2. pharamacoGenomicsFunctions.R
### Inputs:
1. 
### Outputs:
1. 
### Description

## 6 Identify protein chemical interactions and identify compounds

### Scripts: 
1. GetProteinInteractionsAndCompounds.rmd
### Inputs:
1. 
### Outputs:
1. 
### Description
