library(GenomicRanges)
library(ENCODExplorer)
library(rtracklayer)
library(ChIPpeakAnno)


setwd("/data/DaveJames/OC_Epigenetics_Paper/MAGIC_Enhancers/")
# Get genehancer data
gh <- read.table("/data/DaveJames/OC_Epigenetics_Paper/genehancer.gff", stringsAsFactors = F, header = T, sep = "\t")
ghGr <- toGRanges(gh, format = "GFF")
# Get TF files =====================================================================================================================
# encode_df <- get_encode_df()
# Get HG38 results
query_results = queryEncodeGeneric(assembly="GRCh38",
                                   file_format="^bed$", output_type="optimal IDR thresholded peaks",
                                   fixed=FALSE, assay = "TF ChIP-seq")
# 
# # query_results = queryEncodeGeneric(assembly="GRCh38",
# #                                    file_format="^bed$", output_type="optimal IDR thresholded peaks",
# #                                    fixed=FALSE, assay = "TF ChIP-seq", investigated_as = c("transcription factor", "cofactor", 
# #                                                                                            "RNA polymerase complex"))
# 
# # Get Hg19 results
# query_resultsHg19 = queryEncodeGeneric(assembly="hg19",
#                                    file_format="^bed$", output_type="optimal IDR thresholded peaks",
#                                    fixed=FALSE, assay = "TF ChIP-seq", investigated_as = "transcription factor")
# 
files <- data.frame(File = query_results$file_accession, TF = query_results$target, Biosample = query_results$biosample_name)
# #
write.table(files, "downloaded_TF_Tracks_hg19.csv")
write.table(query_results, "query_result_TF_ENCODE.csv")
# #
# # Creatw file for downloading bed files with curl
# dwnLoadFiles <- sapply(files$File, function(x){paste0("https://www.encodeproject.org/files/",x , "/@@download/",x,".bed.gz")})
# write.table(as.data.frame(dwnLoadFiles), "/data/DaveJames/OC_Epigenetics_Paper/MAGIC_Enhancers/ENCODE_EXPLORE/download_Data_hg38.txt", quote =  F, row.names =  F, col.names = F)
# 
# # Load data  =====================================================================================================================
# # Find bed files
pksDr <- "/data/DaveJames/OC_Epigenetics_Paper/MAGIC_Enhancers/PeakFiles/"
pkFn <- list.files(pksDr ,pattern = ".bed\\.gz")
# 
# # Load metadata
metaData <- read.table("/data/DaveJames/OC_Epigenetics_Paper/MAGIC_Enhancers/ENCODE_EXPLORE/query_result_TF_ENCODE.csv", stringsAsFactors = F)


# Make MAGIC matrix  ===============================================================================================================

enNo <- unlist(lapply(1:length(ghGr), function(x){paste0("'GENEHANCER_", x, "'")}))
df <- data.frame(GENE = enNo)
cnt <- 0
for (pk in pkFn){
  cnt <- cnt + 1
  print(paste("On peak set", cnt, pk ))
  # Read in peaks file
  tPks <- read.table((paste0(pksDr, pk)))
  colnames(tPks) <- c("chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
  suppressWarnings(tPks <- toGRanges(tPks))

  #tDf <- data.frame(score = rep(0, length(ghGr)))

  # get overlaps
  suppressWarnings(ovlps <- findOverlaps(ghGr, tPks ))

  # Create the identifier
  idFn <- gsub(".bed.gz", "", pk)
  datEntry <- (metaData[metaData$file_accession == idFn,])
  fac <- datEntry$target
  smp <- datEntry$biosample_name
  id <-  paste0(smp, "_", idFn, ":", fac)

  # Get the scores on the doors
  # st <- Sys.time()
  x <-   unlist(lapply(1:length(ghGr), function(dx){ if(dx %in% ovlps@from) {max(tPks$signalValue[ovlps@to[ovlps@from == dx]])} else {0.0}}))
  # make sum of x = 100000
  x <- x*(100000/sum(x))
  x <- round(x, 2)
  # if (pk == pkFn[1]){
  #
  #   df <- data.frame(tNm = x)
  #   colnames(df) <- id
  # } else {
    df$tNm <- x
    colnames(df) <- c(colnames(df)[1:(ncol(df) -1)], id)
  # }


  # ed <- Sys.time()
  # tm <- ed-st

}

write.table(df, "/data/DaveJames/OC_Epigenetics_Paper/MAGIC/MAGIC/Matrices/enhancerMatrix161221.mtx", quote = F, row.names = F, sep = "\t")

saveRDS(df, "/data/DaveJames/OC_Epigenetics_Paper/MAGIC/MAGIC/MAGIC_EnhancerMatrix_AllTFs.rdata")


