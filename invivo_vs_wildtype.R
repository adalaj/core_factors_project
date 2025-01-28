#!/usr/bin/env Rscript

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/Nathan")

library(tidyverse)
library(data.table)
library(Biostrings)
library(pheatmap)

# I separated each tables and saved them in separate csv files based on their table title
# Task 1: To compare similarity between invivo and wildtype sequences 


reshaped_invivo <- fread("reshaped_S6_invivo.csv", sep = ",", header = TRUE)
#90

reshaped_wildtype<- data.frame( wildtype_identifier = c(">wildtype"),
                                wildtype_sequences = c("GTGTGAGGAAAAGTAGTT"),
                                length = c(18),
                                filename = "wildtype")

#fwrite(reshaped_wildtype, "wildtype_sequence.csv")

sigma<- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE)

attempt2<- data.frame((matrix(nrow = 0, ncol=10)))
colnames(attempt2)<- c("query", "subject", "similarity_perc", "NW_score_NCBI", "aligned_length", "match",
                       "identity_perc_NCBI","mismatch",
                       "aligned_query","aligned_subject")


for (i in 1:nrow(reshaped_invivo)){
  # Define query sequence and identifier
  query_id<- reshaped_invivo$filename[i]
  query_seq <- reshaped_invivo$invivo_sequences[i]
  
  #Loop through the rows for subject sequences
  
  for (j in 1:nrow(reshaped_wildtype)){
    #if (i !=j){ #would make sense if i dont want compare identical query and subject
    # Define subject sequence and identifier
    subject_id <- reshaped_wildtype$filename[j]
    subject_seq <- reshaped_wildtype$wildtype_sequences[j]
    
    alignment<- pwalign::pairwiseAlignment(query_seq, subject_seq,substitutionMatrix = sigma,type="global", gapOpening = 1000, gapExtension = 1000)
    
    query<-  query_id
    subject<-  subject_id
    similarity_perc<- round(pwalign::pid(alignment),2) #similarity percent
    NW_score_NCBI<- pwalign::score(alignment) #NW score
    aligned_length<- width(pwalign::alignedPattern(alignment)) #The 2 objects returned by alignedPattern(x) and alignedSubject(x) are guaranteed to have the same shape (i.e. same length() and width())
    match<- pwalign::nmatch(alignment)
    identity_perc_NCBI<- round((match/aligned_length)*100, 2)#as per ncbi
    mismatch<- pwalign::nmismatch(alignment)
    aligned_query<-as.character(pwalign::alignedPattern(alignment)) #aligned query
    aligned_subject<- as.character(pwalign::alignedSubject(alignment)) #aligned subject
    
    
    
    final<- data.frame(query, subject, similarity_perc, NW_score_NCBI, aligned_length, match,
                       identity_perc_NCBI,mismatch,
                       aligned_query,aligned_subject)
    
    attempt2 <- rbind(attempt2, final)
  }
  
}

fwrite(attempt2, "invivo_vs_wildtype_sequences_global_blast_withoutgaps.csv")
#90


globalblast_4c<- attempt2 %>% select(query, subject, identity_perc_NCBI) #4c - 4 columns
globalblast_4c$identifier <- paste(globalblast_4c$query, globalblast_4c$subject, sep = "&")

similarity_matrix <- data.frame(matrix(NA, nrow=0, ncol=1))

colnames(similarity_matrix)<- t(reshaped_wildtype[,4]) #4 is filename

test<- data.frame(matrix(NA, nrow = 0, ncol = 1)) #additional column is for identifier column. This column will  contain query name
colnames(test)<- c( "identifier")

similarity_matrix <- cbind(test, similarity_matrix) #colnames now increased to 91.


for (i in 1:nrow(reshaped_invivo)){
  listofvalues<- list()
  listofvalues <- append(listofvalues, reshaped_invivo$filename[i])
  for (j in 1:length(reshaped_wildtype$filename)){
    k <- paste(reshaped_invivo$filename[i], reshaped_wildtype$filename[j], sep="&") 
    match<- globalblast_4c[globalblast_4c$identifier== k,] #match is a dataframe that only contain one row which provided condition
    listofvalues <- append(listofvalues, match[1,3])
    #if (nrow(match)<1){ #code debug
    #print(k)
    #print(match)}#append function works by adding things row wise 
    #listofvalues [[i]][i]<- match[1,3] # this didnt work but good to know how to add data points to first entry of list
  }
  
  my_row <- do.call(rbind, listofvalues) # do.call is used to execute a function with a list of arguments.
  similarity_matrix[nrow(similarity_matrix)+1,] <- my_row #rbind will not work because my myrow gives as column 
  
} 

fwrite(similarity_matrix,
       "invivo_vs_wildtype_sequences_similarity_matrix_without_gaps.csv")

# if you plan to read the similarity table again then read as data as data frame 


similarity_matrix<- fread("invivo_vs_wildtype_sequences_similarity_matrix_without_gaps.csv", header= TRUE, sep = ",")
similarity_matrix <- as.data.frame(similarity_matrix)
similarity_matrix_new<- similarity_matrix
row.names(similarity_matrix_new)<- similarity_matrix_new$identifier

similarity_matrix_new<- similarity_matrix_new %>% select(wildtype)

similarity_matrix_new<- as.data.frame(similarity_matrix_new)

similarity_matrix_new<- similarity_matrix_new %>% mutate(across(1, as.numeric))


#whenever you want to plot heatmap perform the above steps, as in read similarity matrix, remove first column, transform column name and make them numeric

#rownames does appear in excel as well in R. Deleted all files.
fwrite(similarity_matrix_new,
       "invivo_vs_wildtype_sequences_without_gaps_similarity_matrix_heatmap_input.csv", sep = ",")



invivo_vs_wildtype<- pheatmap(similarity_matrix_new, 
                               cluster_rows = TRUE, 
                               cluster_cols = FALSE, 
                               display_numbers = FALSE,
                               cellwidth = 50,
                               scale = "none", 
                               main="Comparative Analysis of Sequence Percent Identity between Invivo and wildtype sequence")

#dist_object <- as.dist(distance_matrix)
#Warning message:
#In as.dist.default(distance_matrix) : non-square matrix 
# due to which
#Error in hclust(dist_object, method = "average") : dissimilarities of improper length

ggsave( "invivo_vs_wildtype_sequences_without_gaps_heatmap.tiff", 
        plot = invivo_vs_wildtype, height= 11, width = 15, dpi = 600)

similarity_matrix_new$serial_row_number<- 1:90 #if there is only one column then it assume it as list not data.table or frame
heatmap_order<- similarity_matrix_new[invivo_vs_wildtype$tree_row$order,] # This reorders the rows of similarity_matrix according to the clustering order
heatmap_order$sequences <- row.names(heatmap_order) #designated to column 3
heatmap_order$serial_row_number <- 1:90 #new row number according to heatmap
#replaces column 2
heatmap_order$original_row_order<- invivo_vs_wildtype$tree_row$order
#designated to column 4

heatmap_order2_invivo_vs_wildtype <- heatmap_order %>% select(original_row_order, serial_row_number, sequences, wildtype)
fwrite(heatmap_order2_invivo_vs_wildtype, "invivo_vs_wildtype_sequences_without_gaps_similarity_matrix_heatmap_order_output.csv")
mean(heatmap_order2_invivo_vs_wildtype$wildtype)
#[1] 35.37011

