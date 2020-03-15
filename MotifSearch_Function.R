

#--------------------------------------------# Toy Sequences #--------------------------------------------#


long_seq <- "ATGTCGATCGCATCGATCGACCGATCGACGCATCGATCAGCTAGCATCGAATGATGCATCGAATATCGTACGATCGATCGAAGCTAGCAGAAT"
short_seq <- "GAA"

#----------------------------------# Advanced Motif Search Function #-------------------------------------#

adv_motif_search <- function(long_seq, short_seq, similarity = NULL) {

  #---Initializing Required Containers---#
  bool_list <- as.list(c())
  simil_list <- list()
  final_mat <- matrix()
  simil_mat <- matrix()
  
#----------Checks for input of optional argument----------#
    
  if(is.null(similarity)) {
    similarity = 1.0
  } else {
    similarity <- similarity
  }

  #----Checks if 1st sequence is longer than 2nd and checks for correct length conditions----#    
  if (nchar(long_seq) > nchar(short_seq)) {
    
    # Checks for: 
    # length of long_seq, must be less than 1000
    # length of long_seq, must be more than 10
    # length of short_seq, must be less than 10
    if (nchar(long_seq) < 1000 & 
        nchar(long_seq) > 10 & 
        nchar(short_seq) < 10 & 
        nchar(short_seq) != 0) 
    {
      #---Determines positions on long seq to slice, by the length of short seq---#
      slice_position <- seq(from = 1, to = nchar(long_seq), by = nchar(short_seq)) 
      
      #---Slicing function to slice the long seq by length of short seq---#
      string_sub <- function(k) substr(long_seq, k, k + (nchar(short_seq)-1)) 
      
      #---Applying slicing function to slice according to a given slice positions---#
      seq_sliced <- sapply(slice_position, string_sub) 
      
      #---Detect sequences in long seq that contain the short seq---#
      motif_posit <- which(seq_sliced == short_seq) 
      
      #---Number of motif matches---#
      match_count <- length(motif_posit) 
      
      #---Splits string slices into individual characters---#
      seq_sliced_2 <- strsplit(seq_sliced, character(0)) 
      short_split <- strsplit(short_seq, character(0))
      
      #---Initializing Binary Matrix---#
      binary_mat <- matrix(data = NA, nrow = length(short_split[[1]]))[,-1]

      #---Generates a binary matrix of match scores---#
      for (i in 1:length(seq_sliced_2)) {
        for (j in 1:length(short_split[[1]])) { # testing each slice against the short seq
          if (seq_sliced_2[[i]][j] == short_split[[1]][j]) { # Gives a score of 1 when nucleotide matches
            bool_list <- append(bool_list, "1")
            
          } else if (seq_sliced_2[[i]][j] != short_split[[1]][j]) { # Gives a score of 0 when nucleotide does not match
            bool_list <- append(bool_list, "0")
          }
        }
        #---Converts the accumulated list into a matrix---#
        binary_mat <- matrix(bool_list, nrow = length(short_split[[1]])) 
        
        #---Convert characters in matrix into numeric values---#
        binary_mat <- apply(binary_mat, c(1,2), as.numeric) 
      }
      
      #---Finding total sum of match scores within each test slice---#
      test_match_count <- as.numeric(colSums(binary_mat)) 
      
      #---Calculates the similarity indexes for all slices and outputs into a list---#
      for (p in 1:length(test_match_count)) {
        simil_index <- test_match_count[p] / length(short_split[[1]]) 
        simil_list <- append(simil_list, simil_index)
      } 
      
      #----------Combining all lists into a presentable data frame-------------#
      
      seq_mat <- t(data.frame(seq_sliced_2, stringsAsFactors = FALSE))
      simil_mat <- t(data.frame(simil_list, stringsAsFactors = FALSE))
      rownames(simil_mat) <- rep(1:nrow(simil_mat))
      colnames(simil_mat) <- "Similarity Index"
      final_mat <- cbind(seq_mat, simil_mat) # matrix with all sequences and similarities
      rownames(final_mat) <- rep(1:nrow(simil_mat))
      
      #---Determines if the similarity scores fall within the similarity threshold---#
      simil_mat <- simil_mat[which(simil_mat >= similarity & simil_mat <= 1.0),]
    }
  } else {
    #---Prints an error if the 1st Sequence is longer than the 2nd Sequence---#
    print("Error: 1st Argument must be longer than 2nd Argument")
  }
  #---Building a Final List that combines all the computed information---#
  final_list <- list(noquote(final_mat) ,noquote(simil_mat), paste("Number of matches = ", sum(simil_mat)))
  
  return(print(final_list)) # matrix of the motif match positions is returned to the user
}
motif_result_toy <- adv_motif_search(long_seq, short_seq, similarity = 1.0)


#---------------------# Applying function on p53 tumor protein in homo sapiens # -------------------------#

#---Homo sapiens tumor protein p53 (TP53), RefSeqGene (LRG_321) on chromosome 17---#

# The following sequence is a segment of DNA taken from the RefSeqGene (LRG_321) on chromosome 17
# of the Tumor Protein p53 (TP53) found in Homo sapiens.

p53_seq <- 
"GATGGGATTGGGGTTTTCCCCTCCCATGTGCTCAAGACTGGCGCTAAAAGTTTTGAGCTTCTCAAAAGTC
TAGAGCCACCGTCCAGGGAGCAGGTAGCTGCTGGGCTCCGGGGACACTTTGCGTTCGGGCTGGGAGCGTG
CTTTCCACGACGGTGACACGCTTCCCTGGATTGGGTAAGCTCCTGACTGAACTTGATGAGTCCTCTCTGA
GTCACGGGCTCTCGGCTCCGTGTATTTTCAGCTCGGGAAAATCGCTGGGGCTGGGGGTGGGGCAGTGGGG
ACTTAGCGAGTTTGGGGGTGAGTGGGATGGAAGCTTGGCTAGAGGGATCATCATAGGAGTTGCATTGTTG
GGAGACCTGGGTGTAGATGATGGGGATGTTAGGACCATCCGAACTCAAAGTTGAACGCCTAGGCAGAGGA
GTGGAGCTTTGGGGAACCTTGAGCCGGCCTAAAGCGTACTTCTTTGCACATCCACCCGGTGCTGGGCGTA
GGGAATCCCTGAAATAAAAGATGCACAAAGCATTGAGGTCTGAGACTTTTGGATCTCGAAACATTGAGAA
CTCATAGCTGTATATTTTAGAGCCCATGGCATCCTAGTGAAAACTGGGGCTCCATTCCGAAATGATCATT
TGGGGGTGATCCGGGGAGCCCAAGCTGCTAAGGTCCCACAACTTCCGGACCTTTGTCCTTCCTGGAGCGA
TCTTTCCAGGCAGCCCCCGGCTCCGCTAGATGGAGAAAATCCAATTGAAGGCTGTCAGTCGTGGAAGTGA"

motif <- "GAA"

motif_result_p53 <- adv_motif_search(p53_seq, motif, 1.0)

#------------ Motif Match Positions above Similarity Threshold ----------------#
# [[2]]
# 84 170 189 203 249 253 
# 1   1   1   1   1   1 

#------------ Number of Matches above Similarity Threshold --------------------#
# [[3]]
# [1] "Number of matches =  6"
