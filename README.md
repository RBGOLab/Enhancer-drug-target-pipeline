# Enhancer Drug Targets
Identification of ovarian cancer subtype specific drug targets using epigenomic and transcriptomic data

## 1 Create Enhancer based MAGIC matrix

### Scripts: 
1. CreateEnhancerMAGICMatrix.r
### Inputs:
1. GeneHancer gff file
2. TF ChIP-seq bed files
3. TF ChIP-seq meta data file
### Outputs:

ChIP-seq TF peak files and meta data can be downloaded from ENCODE experiment matrix (https://www.encodeproject.org/matrix/?type=Experiment&control_type!=*&status=released&perturbed=false&assay_title=TF+ChIP-seq&target.investigated_as=transcription+factor&assembly=GRCh38&files.file_type=bed+narrowPeak)

Once files are selected, use the Download button to download a text file containing the url of each peak file and a corresponding meta data file. Use curl to download bed.gz files as instructed.

To generate a 'MAGIC matrix' file based on GeneHancer known enhancer locations, use the r script 'CreateEnhancerMAGICMatrix.r'. This requires a copy of the geneHancer database in gff format. The 2017 version can be downloaded directly from GeneCards https://www.genecards.org/, and the latest version can be obtained by request.  

