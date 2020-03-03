#' @title Converts list gRNAs to primers for cloning into pcfd5 
#'
#' @description This function converts a list gRNAs to primers for cloning into pcfd5 plasmid following the protocol of Port & Bullock (2016). 
#' The input is a csv with the gene name in column 1 and the gRNAs in columns 2 to 5. The output is a spreadsheet with the appropriate primers designed to fill a 96 well plate. 
#' Forward and reverse primers for the same PCR reaction asre assigned to the same well. If you are ordering individual primers in tubes, ignore output columns 1 & 2. 
#' This function requires the Biostrings and tidyr packages.
#'
#' @param filepath
#'
#' @return NULL
#'
#' @examples gRNA2primers('~/Desktop/sample.csv')
#'
#' @export gRNA2primers

library(Biostrings)
library(tidyr)

gRNA2primers <- function(filepath){
  file <- read.csv(filepath, head=T, sep=',')
  
  # CREATE OUTPUT STRUCTURE
  primers_out <- data.frame(matrix(nrow = nrow(gRNA_out), ncol = 7))
  colnames(primers_out) <- c("gene", "PCR1_F", "PCR1_R", "PCR2_F", "PCR2_R", "PCR3_F", "PCR3_R")
  
  for (i in 1:nrow(gRNA_out)){
    ## DEFINE GUIDES
    primers_out[i,1] = as.character(file[i,1])
    g1 = DNAStringSet(as.character(file[i,2]))
    g2 = DNAStringSet(as.character(file[i,3]))
    g3 = DNAStringSet(as.character(file[i,4]))
    g4 = DNAStringSet(as.character(file[i,5]))
    ## DEFINE PRIMERS
    primers_out[i,2] <- as.character(paste0("GCGGCCCGGGTTCGATTCCCGGCCGATGCA", as.character(g1), "GTTTTAGAGCTAGAAATAGCAAG")) ##PCR1fwd
    primers_out[i,3] <- as.character(paste0(as.character(reverseComplement(g2)), "TGCACCAGCCGGGAATCGAACCC")) ##PCR1rev:
    primers_out[i,4] <- as.character(paste0(as.character(g2), "GTTTTAGAGCTAGAAATAGCAAG")) ##PCR2fwd:
    primers_out[i,5] <- as.character(paste0(as.character(reverseComplement(g3)), "TGCACCAGCCGGGAATCGAACCC"))	##PCR2rev:
    primers_out[i,6] <- as.character(paste0(as.character(g3), "GTTTTAGAGCTAGAAATAGCAAG")) ##PCR3fwd:
    primers_out[i,7] <- as.character(paste0("ATTTTAACTTGCTATTTCTAGCTCTAAAAC", as.character(reverseComplement(g4)), "TGCACCAGCCGGGAATCGAACCC")) ##PCR3rev
  }
  
  # CONVERT TO LONG FORMAT
  long <- gather(primers_out, primer, seq, PCR1_F:PCR3_R, factor_key=TRUE)
  long <- long[order(long$gene, long$primer),]
  
  
  # CREATE OUTPUT DATAFRAME
  out <- data.frame(matrix(nrow = nrow(long), ncol = 7))
  colnames(out) <- c("Plate_name", "Well_Position", "Sequence_Name", "Sequence", "Scale", "Purification", "Normalization_Style")
  
  # CALCULATE WELL POSITION
  c = toupper(letters[1:8])
  col = rep(toupper(letters[1:8]), each=24)
  r = rep(seq(1, 12, 3), each=6)
  add = rep(seq(0, 2), each=2)
  radd = rep((r+add), times=8)
  
  # COMBINE DATA INTO THE PRIMER ORDER FORM
  for (ii in 1:nrow(long)){
    out[ii, 1] <- paste0("plate_", ceiling(ii/32) )
    out[ii, 2] <- paste0(col[ii],as.character(radd[ii]))
    out[ii, 3] <- paste0(long$gene[ii], "_", long$primer[ii])
    out[ii, 4] <- long$seq[ii]
    out[ii, 5] <- "25 nmole DNA Plate Oligo"
    out[ii, 6] <- "Standard Desalting"
    out[ii, 7] <- "Concentration and Quality"
  }
    
  # SAVE FILES
  filename <- as.character(paste0(substr(filepath, 1, nchar(filepath)-4), "_primer_plate.csv"))
  write.csv(out, filename, row.names = FALSE)
   
  }
  