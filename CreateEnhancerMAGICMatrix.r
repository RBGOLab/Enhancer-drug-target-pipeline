library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno)

wd <- "" # location to write results
pksDr <- "" # directory containing TF ChIP-seq bed.gz files downloaded from ENCODE
metaDataDr <- "" # location of TF peak metadata

# Get genehancer data
gh <- read.table("genehancer.gff", stringsAsFactors = F, header = T, sep = "\t") # set path to Genehancer gff
ghGr <- toGRanges(gh, format = "GFF")


## Load data  =====================================================================================================================
## Find bed files
pkFn <- list.files(pksDr ,pattern = ".bed\\.gz")
 
## Load metadata
metaData <- read.table(metaDataDr, stringsAsFactors = F, sep = "\t", header = T)


# Make MAGIC matrix  ===============================================================================================================
#browser()
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
  datEntry <- (metaData[metaData$File.accession == idFn,])
  fac <- datEntry$Experiment.target
  smp <- datEntry$Biosample.term.name
  id <-  paste0(smp, "_", idFn, ":", fac)

  # Get the scores on the doors
  x <-   unlist(lapply(1:length(ghGr), function(dx){ if(dx %in% ovlps@from) {max(tPks$signalValue[ovlps@to[ovlps@from == dx]])} else {0.0}}))
  # make sum of x = 100000
  x <- x*(100000/sum(x))
  x <- round(x, 2)
    df$tNm <- x
    colnames(df) <- c(colnames(df)[1:(ncol(df) -1)], id)

}

write.table(df, file.path(wd, "MAGIC_EnhancerMatrix.mtx"), quote = F, row.names = F, sep = "\t")


