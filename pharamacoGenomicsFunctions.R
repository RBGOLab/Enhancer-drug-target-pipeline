## Function files for analysis of upregulated TF's in promoter regiosn and enahncer regions

# Function for comparing RNA_seq data with MAGIC factors ================================================================================= ###

library(readxl)

facRNA_Ovlp_flag_GFP <- function(rnaFile, facFile, pMAGIC = 0.1, pRNA = 0.05){
  
  ## Load differential analysis from clear cell RNA-seq tissue data
  RNASeq <- read.csv(rnaFile)
  
  ## Load enriched factors from clear cell
  UpRegFacMAGIC <- read_xls(facFile)
  UpRegFacMAGIC$Factor <- gsub(pattern = "3xFLAG-", replacement = "",x =  UpRegFacMAGIC$Factor) # get rid of the 3xFlag- tag
  UpRegFacMAGIC$Factor <- gsub(pattern = "eGFP-", replacement = "",x =  UpRegFacMAGIC$Factor) # get rid of the eGFP- tag
  UpRegFacMAGIC <- UpRegFacMAGIC[UpRegFacMAGIC$`Corrected p` < pMAGIC ,]
  
  
  
  ## Differential analysis p < 0.05
  #browser()
  RNASeq$padj[is.na(RNASeq$padj)] <- 1
  RNASeq_FDRLs0p05 <- RNASeq[RNASeq$padj < pRNA, ]
  
  ovlpFacs <- intersect(RNASeq_FDRLs0p05$Gene.Symbol, UpRegFacMAGIC$Factor)
  
  upRegFacRNASeq <- data.frame(factor = ovlpFacs, 
                               RNASeqfc = RNASeq_FDRLs0p05$log2FoldChange[RNASeq_FDRLs0p05$Gene.Symbol %in% ovlpFacs])
  #browser()
  upRegFacRNASeq$tfEnrichScore <- unlist(lapply(ovlpFacs, function(x){mean(UpRegFacMAGIC$Score[UpRegFacMAGIC$Factor == x])} ))
  
  upRegFacRNASeq <- upRegFacRNASeq[order(upRegFacRNASeq$RNASeqfc, decreasing = T),]
  
}


## Function for finding complexes in CORUM DB containing TF's ========================================================================== ###

library(org.Hs.eg.db)


library(rjson)
cmplxDf <- fromJSON(file = "/data/DaveJames/OC_Epigenetics_Paper/CORUM_Complexes/coreComplexes.json")

## function to get complexes
GetCmplexes <- function(cmplxDf, lst){
  # lst is a list of proteins
  # cmplxDf is a datframe of cmplxes from corum
  
  df <- NULL
 # browser()
  for (lx in 1:length(lst)){
    
    p <- mapIds(org.Hs.eg.db, lst[lx], "ENTREZID", "SYMBOL")
    l <- logical(length(cmplxDf))
    
    
    for (dx in 1:length(cmplxDf)){
      
      if (p %in% unlist(strsplit(cmplxDf[[dx]]$`subunits(Entrez IDs)`, ";")) & cmplxDf[[dx]]$Organism == "Human"){
        l[dx] = TRUE
        
      }
      
      
    }
    #browser()
    tDf <- data.frame(TF = rep(lst[lx], sum(l)), EntrezID = rep(p, sum(l)),
                      Complex_Name = unlist(lapply(which(l),function(x) unlist(cmplxDf[[x]]$ComplexName))), 
                      PubMedID = unlist(lapply(which(l),function(x) unlist(cmplxDf[[x]]$`PubMed ID`))),
                      Proteins = unlist(lapply(which(l),function(x) unlist(cmplxDf[[x]]$`subunits(Gene name)`))), stringsAsFactors = F
    )
    
    if (nrow(tDf) == 0){
      
      tDf <- data.frame(TF = lst[lx], EntrezID = p,
                        Complex_Name = "No complex", 
                        PubMedID = "No complex",
                        Proteins = lst[lx]
                        , stringsAsFactors = F
      ) 
      
    }

    if (is.null(df)){
      df <- tDf
      
    } else {
      
      df <- rbind(df, tDf)
      
    }
    
  }
  
  df
  
}

## Make a protein list from complexes

getProteinsListFromComplexes <- function(cmplxDF){
  
  cmplxDfProteins <- (paste(as.character(cmplxDF$Proteins), collapse = ';'))
  
  cmplxDfProteinsAll <- (unlist(strsplit(cmplxDfProteins, ";")))
  indx <- !duplicated(cmplxDfProteinsAll)
  
  cmplxDfProteins <- unique(unlist(strsplit(cmplxDfProteins, ";")))
  
  complexesCompOf <- character(length(cmplxDfProteins))
  
  for (dx in 1:length(cmplxDfProteins)){


    complexesCompOf[dx] <- paste(cmplxDF$Complex_Name[grepl(cmplxDfProteins[dx], cmplxDF$Proteins)], collapse = " | ")

  }
  # 
  
  #browser()
  
  cmplxDfBsMn <- (paste(as.character(cmplxDF$BaseMean), collapse = '|'))
  cmplxDfBsMn <- unique(unlist(strsplit(cmplxDfBsMn, "\\|")))  
 
  cmplxDflfc <- (paste(as.character(cmplxDF$log2fc), collapse = '|'))
  cmplxDflfc <- unlist(strsplit(cmplxDflfc, "\\|"))[indx]
  #cmplxDflfc <- unique(unlist(strsplit(cmplxDflfc, "\\|")))  
  
  cmplxDfnrmCnts <- (paste(as.character(cmplxDF$nrmExprss), collapse = '|'))
  cmplxDfnrmCnts <- unique(unlist(strsplit(cmplxDfnrmCnts, "\\|"))) 
   
  cmplxDfProteinsDF <- data.frame(Protein = cmplxDfProteins, log2fc = cmplxDflfc, BaseMean = cmplxDfBsMn, 
                                  normCounts = cmplxDfnrmCnts, componentOf = complexesCompOf) 
  
  
}

## RNA-seq DEG data to filter genes with low expression ================================================================================================= ##
tfsRNA_Seq <- function(MAGIC_res, RNA_GnTb, MAGIC_Pth, RNA_GnCntPth, pMAGIC = 0.01, pRNA = 0.05, gnCntThres = 10){

  UpRegFacMAGIC <- read.csv(file.path(MAGIC_Pth, MAGIC_res))$x

  ## Expression analysis ================================================================================================= ##
  
  gnTbPths <- list.files(RNA_GnCntPth, paste0(RNA_GnTb, ".*ReadsPerGene.out.tab"), full.names = T)
  
  cnt <- 1
  # Get the gene count results
  for (fx in gnTbPths){
    
    if (cnt == 1){
      
      gnCnts <- read.table(fx, colClasses = c("character", "integer", "NULL", "NULL"),  row.names = 1, skip = 4)
      
      
    } else {
      
      
      gnCnts <- cbind(gnCnts, read.table(fx, colClasses = c("character", "integer", "NULL", "NULL"),  row.names = 1, skip = 4))
      
    }
    
    cnt <- cnt + 1 
    
  }
  
  rownames(gnCnts) <- gsub("\\..*", "", rownames(gnCnts))
  
  # Get the genes in ENSEMBL form
  # MAGIC_Gns2Ensembl <- ensembldb::select(EnsDb.Hsapiens.v86, keys = UpRegFacMAGIC$Factor, keytype = "SYMBOL", columns = c("GENEID"))
  MAGIC_Gns2Ensembl <- ensembldb::select(EnsDb.Hsapiens.v86, keys = UpRegFacMAGIC, keytype = "SYMBOL", columns = c("GENEID"))
  MAGIC_Gns2Ensembl <- ensembldb::select(EnsDb.Hsapiens.v86, keys = UpRegFacMAGIC, keytype = "SYMBOL", columns = c("GENEID"))
  
  gnCntsMAGIC <- gnCnts[rownames(gnCnts) %in% MAGIC_Gns2Ensembl$GENEID, ]
  
  #browser()
  gnCntsMAGIC$mean <- rowMeans(gnCntsMAGIC)
  gnCntsMAGIC$SYMBOL <- MAGIC_Gns2Ensembl$SYMBOL[na.omit(match(rownames(gnCntsMAGIC) , MAGIC_Gns2Ensembl$GENEID))]

  
  gnCntsMAGIC$GENE_ID_Check <- MAGIC_Gns2Ensembl$GENEID[na.omit(match(rownames(gnCntsMAGIC) , MAGIC_Gns2Ensembl$GENEID))]
  gnCntsMAGIC1 <- gnCntsMAGIC
  gnCntsMAGIC <- gnCntsMAGIC[gnCntsMAGIC$mean > gnCntThres, ] # remove any genes with mean gene count less than 10
  gnCntsMAGIC <- gnCntsMAGIC[!(rowSums(gnCntsMAGIC[, 1:5] < gnCntThres) > 2) , ] # remove any genes with 3 or more sample counts less than 10
  
  
  gnCntsMAGIC_TFRemoved <- gnCntsMAGIC1$SYMBOL[!(gnCntsMAGIC1$SYMBOL %in% gnCntsMAGIC$SYMBOL)]

  gnCntsMAGIC
  
}

## ================================================================================================= ##

checkProteinExpression <- function(lst, RNA_GnCntPth, RNA_GnTb, outputPth = NULL){
  lst1 <- lst
  #browser()
  # Load the gene count tables
  gnTbPths <- list.files(RNA_GnCntPth, paste0(RNA_GnTb, ".*ReadsPerGene.out.tab"), full.names = T)
  
  cnt <- 1
  # Get the gene count results
  for (fx in gnTbPths){
    
    if (cnt == 1){
      
      gnCnts <- read.table(fx, colClasses = c("character", "integer", "NULL", "NULL"),  row.names = 1, skip = 4)
      
      
    } else {
      
      
      gnCnts <- cbind(gnCnts, read.table(fx, colClasses = c("character", "integer", "NULL", "NULL"),  row.names = 1, skip = 4))
      
    }
    
    cnt <- cnt + 1 
    
  }
  #browser()
  rownames(gnCnts) <- gsub("\\..*", "", rownames(gnCnts))
  
  #### Get the ENSEMBL IDs
  symbol2Ensembl <- ensembldb::select(EnsDb.Hsapiens.v86, keys = lst$Proteins, keytype = "SYMBOL", columns = c("GENEID", "SYMBOL"))
  
  gnCntsFilt <- gnCnts[rownames(gnCnts) %in% symbol2Ensembl$GENEID, ]  
  gnCntsFilt$SYMBOL <- symbol2Ensembl$SYMBOL[na.omit(match(rownames(gnCntsFilt) , symbol2Ensembl$GENEID))]
  gnCntsFilt$GENE_ID_Check <- symbol2Ensembl$GENEID[na.omit(match(rownames(gnCntsFilt) , symbol2Ensembl$GENEID))]
  
  gnCntsFilt <- gnCntsFilt[!(rowSums(gnCntsFilt[, 1:(ncol(gnCntsFilt)-2)] < 10) >2), ] # remove any genes with 3 or more counts less than 
  
  
  uniqueCmpl <- unique(lst$PubMedID)
  
  for (dx in 1:length(uniqueCmpl)){
    
    t <- lst[lst$PubMedID == uniqueCmpl[dx] , ]
    
    if (!all(t$Proteins %in% gnCntsFilt$SYMBOL)){ # if any of the proteins have low expression remove them
      
      lst <- lst[!(lst$PubMedID == uniqueCmpl[dx]), ]
      
    }
    
    
  }
  
  #browser()
  lstRm <- unique(lst1$Complex_Name[!(lst1$Complex_Name %in% lst$Complex_Name)  ])
  
  if(!is.null(outputPth)){
  write.table(lstRm, paste0(outputPth, "/", RNA_GnTb, "_ComplexesRemoved.tsv" ), quote = F, sep = "\t")
  }

  lstOut <- cbind(lst, gnCntsFilt[match(lst$Proteins, gnCntsFilt$SYMBOL),])

  # If there all complexes for a TF are removed, then replace it with a no complex result
  if (any(!(unique(lst1$TF) %in% unique(lst$TF) ))){
    
    tfOut <- unique(lst1$TF)[!(unique(lst1$TF) %in% unique(lst$TF))] 
    #tfOutidx <- which(!(unique(lst1$TF) %in% unique(lst$TF)) )
    
    for (ttdx in length(tfOut)){
      
      ensId <- symbol2Ensembl[symbol2Ensembl$SYMBOL == tfOut[ttdx],]$GENEID
      tDf <- data.frame(TF = tfOut[ttdx], EntrezID = 0, Complex_Name = "No Expressed Complexes", PubMedID = "No Expressed Complexes", 
                        Proteins = tfOut[ttdx], V2 = 0, V2 = 0, V2 = 0, V2 = 0, V2 = 0, SYMBOL = tfOut[ttdx], GENE_ID_Check = ensId, stringsAsFactors = F)
      rownames(tDf) <- ensId
      colnames(tDf) <- c("TF", "EntrezID", "Complex_Name", "PubMedID", "Proteins", "V2", "V2", "V2", "V2", "V2", "SYMBOL", "GENE_ID_Check")
      lstOut <- rbind(lstOut, tDf)
    }
  }
  lstOut
}


