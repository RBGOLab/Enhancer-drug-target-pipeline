# Enhancer Drug Targets
Scripts used for identification of ovarian cancer subtype specific drug targets using epigenomic and transcriptomic data from **PAPER REFERENCE TO GO HERE**. Example input and output files are provided, and comprehensive isntructions for each script and the pipeline steps which the correspond to are outlined below.

![Pharmacogenomics Pipeline](https://github.com/RBGO-Lab/EnhancerDrugTargetPipeline/blob/main/Figure1.png)
TF target pipeline 

## 1) Create Enhancer based MAGIC matrix (Pipe line Step 5)

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

## 2) Differential Analysis of H3K27ac Enhancer Sites (Pipe line Step 3 + 4)

### Scripts: 
1. DifferentialAnalysisOfEnhancerSites.rmd
### Inputs:
1. GeneHancer gff file
2. BAM and corresponding MACS2 H3K27ac ChIP peak file from cancer cell lines
3. BAM and corresponding MACS2 H3K27ac ChIP peak file from cancer tissue
4. BAM and corresponding MACS2 H3K27ac ChIP peak file from healthy cell lines
5. BAM and corresponding MACS2 H3K27ac ChIP peak file from healthy tissue (if available)
### Outputs:
1. Upregulated enhancer site text file (containing two columns with list of background enhancers and upregulated enhancers)
### Description

Identifies enhancer sites annoatated in GeneHancer database which are signifcantly upregulated in cancer compared with corresponding healthy samples from 
H3K27ac histone ChIP-seq datasets. DiffBind package used for differential analysis (Stark R, Brown G (2011). DiffBind: differential binding analysis of ChIP-Seq peak data, https://bioconductor.org/packages/release/bioc/html/DiffBind.html) . The output from this script is a list of upregulated enhancer locations and a background list of all enhancers within the cancer datasets.   

## 3) Identify Enriched TF binding in differentially expressed H3K27ac enhancer sites (Pipe line Step 5)

### Scripts: 
1. MAGIC_1_1.py
### Inputs:
1. MAGIC enhancer MTX file
2. Upregulated enhancer site text file
### Outputs:
1. Summary of upregulated TFs xls
2. Details of upregulated TFs xls
3. For further details see https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007800 and https://github.com/asroopra/MAGIC
### Description

MAGIC algorithm (Roopra A (2020) MAGIC: A tool for predicting transcription factors and cofactors driving gene sets using ENCODE data. PLOS Computational Biology 16(4): e1007800., https://github.com/asroopra/MAGIC) used to identify TFs with binding significantly enriched in upregulated enhancer sites. 

## 4) Filter MAGIC results by enriched Motifs (Pipe line Step 6)

### Scripts: 
1. MotifSearch.rmd
### Inputs:
1. GeneHancer gff file
2. MAGIC Summary of upregulated TFs xls
3. Upregulated enhancer site text file
### Outputs:
1. csv file with TFs and P value indicating whether its motif is enriched within the upregulated enhancer sites
### Dependencies
1. Local installation of meme suite https://meme-suite.org/meme/doc/install.html?man_type=web
### Description
Each enriched TF identified by MAGIC is put through an empircal test to check if its motif is significantly enriched within the upregulated enhancer sites. This acts as a further filtering step for TF binding.   

## 5) Identify complexes enriched TFs are components, and filter by expression data (Pipe line Step 7 + 8)

### Scripts: 
1. ComplexComponentSearchAndFilter.rmd
2. pharamacoGenomicsFunctions.R
### Inputs:
1. CORUM protein complex database in json format (can be downloaded from http://mips.helmholtz-muenchen.de/corum/#download)
2. List of TFs with enriched binding in csv format
3. Directory containing gene count tables from STAR aligner (Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635. Epub 2012 Oct 25. PMID: 23104886; PMCID: PMC3530905., https://github.com/alexdobin/STAR)
### Outputs:
1. Complex centred CSV file containing complexes which contain enriched TFs and all component proteins
2. Protein centred CSV of all proteins identified in complexes 

### Description
Identifies all complexes in CORUM DB ( Tsitsiridis G, Steinkamp R, Giurgiu M, Brauner B, Fobo G, Frishman G, Montrone C, Ruepp A. CORUM: the comprehensive resource of mammalian protein complexes–2022. Nucleic Acids Res. 2022 Nov 16. doi: 10.1093/nar/gkac1015. ,http://mips.helmholtz-muenchen.de/corum/) which contain one or more of the TFs with enriched binding. Gene count tables from RNA-seq experiments, corresponding to the cancer of interest, are used to filter out complexes which contain 

## 6) Identify protein chemical interactions and identify compounds (Pipe line Step 8)

### Scripts: 
1. GetProteinInteractionsAndCompounds.rmd
### Inputs:
1. Protein centred CSV of all proteins identified in complexes from ComplexComponentSearchAndFilter.rmd
2. Stitch database TSV containing all protein/chemical interactions (can be downloaded from STITCHwebsite http://stitch.embl.de/cgi/download.pl?UserId=5ns76rjpXdqG&sessionId=ENQw3vOuwuFK&species_text=Homo+sapiens)
3. Chembl database (installed locally)
### Outputs:
1. CSV file containing compounds which target proteins in identified complexes including trial status
### Dependencies
1. Local installation of chembl MYSQL database
### Description
Using STITCH database for protein/chemical interaction scores (Kuhn M, von Mering C, Campillos M, Jensen LJ, Bork P. STITCH: interaction networks of chemicals and proteins. Nucleic Acids Res. 2008 Jan;36(Database issue):D684-8. doi: 10.1093/nar/gkm795. Epub 2007 Dec 15. PMID: 18084021; PMCID: PMC2238848., http://stitch.embl.de/cgi/input.pl?UserId=5ns76rjpXdqG&sessionId=ENQw3vOuwuFK) and Chembl database for clincial trial data, find potential compounds which target proteins within identified complexes 

