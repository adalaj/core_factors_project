#!/usr/bin/env Rscript

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/Nathan")

library(tidyverse)
library(data.table)
library(Biostrings)
library(pheatmap)

# I separated each tables and saved them in separate csv files based on their table title
# Task 1: To compare similarity between invitro and invivo sequences 

reshaped_invivo<- fread("reshaped_S6_invivo.csv", sep = ",", header = TRUE)
#90 rows

reshaped_invitro <- fread("reshaped_S3_R7_invitro.csv", sep = ",", header = TRUE)
#516


sigma<- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE)

attempt2<- data.frame((matrix(nrow = 0, ncol=10)))
colnames(attempt2)<- c("query", "subject", "similarity_perc", "NW_score_NCBI", "aligned_length", "match",
                       "identity_perc_NCBI","mismatch",
                       "aligned_query","aligned_subject")


for (i in 1:nrow(reshaped_invitro)){
  # Define query sequence and identifier
  query_id<- reshaped_invitro$filename[i]
  query_seq <- reshaped_invitro$invitro_sequences[i]
  
  #Loop through the rows for subject sequences
  
  for (j in 1:nrow(reshaped_invivo)){
    #if (i !=j){ #would make sense if i dont want compare identical query and subject
    # Define subject sequence and identifier
    subject_id <- reshaped_invivo$filename[j]
    subject_seq <- reshaped_invivo$invivo_sequences[j]
    
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

fwrite(attempt2, "invitro_vs_invivo_sequences_global_blast_withoutgaps.csv")
#46440


globalblast_4c<- attempt2 %>% select(query, subject, identity_perc_NCBI) #4c - 4 columns
globalblast_4c$identifier <- paste(globalblast_4c$query, globalblast_4c$subject, sep = "&")

similarity_matrix <- data.frame(matrix(NA, nrow=0, ncol=90))

colnames(similarity_matrix)<- t(reshaped_invivo[,4]) #4 is filename

test<- data.frame(matrix(NA, nrow = 0, ncol = 1)) #additional column is for identifier column. This column will  contain query name
colnames(test)<- c( "identifier")

similarity_matrix <- cbind(test, similarity_matrix) #colnames now increased to 91.


for (i in 1:nrow(reshaped_invitro)){
  listofvalues<- list()
  listofvalues <- append(listofvalues, reshaped_invitro$filename[i])
  for (j in 1:length(reshaped_invivo$filename)){
    k <- paste(reshaped_invitro$filename[i], reshaped_invivo$filename[j], sep="&") 
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
       "invitro_vs_invivo_sequences_similarity_matrix_without_gaps.csv")

# if you plan to read the similarity table again then read as data as data frame 


similarity_matrix<- fread("invitro_vs_invivo_sequences_similarity_matrix_without_gaps.csv", header= TRUE, sep = ",")
similarity_matrix <- as.data.frame(similarity_matrix)
similarity_matrix_new<- similarity_matrix
row.names(similarity_matrix_new)<- similarity_matrix_new$identifier
similarity_matrix_new<- similarity_matrix_new[,-c(1)]

similarity_matrix_new<- similarity_matrix_new %>% mutate(across(1:90, as.numeric))
#whenever you want to plot heatmap perform the above steps, as in read similarity matrix, remove first column, transform column name and make them numeric

#rownames does appear in excel as well in R. Deleted all files.
fwrite(similarity_matrix_new,
       "invitro_vs_invivo_sequences_without_gaps_similarity_matrix_heatmap_input.csv", sep = ",")




invitro_vs_invivo<- pheatmap(similarity_matrix_new, 
                             cluster_rows = TRUE, 
                             cluster_cols = TRUE, 
                             scale = "none", 
                             main="Comparative Analysis of Sequence Percent Identity between Invitro (R7) and Invivo sequences")

#dist_object <- as.dist(distance_matrix)
#Warning message:
  #In as.dist.default(distance_matrix) : non-square matrix 
# due to which
#Error in hclust(dist_object, method = "average") : dissimilarities of improper length

ggsave( "invitro_vs_invivo_sequences_without_gaps_heatmap.tiff", 
        plot = invitro_vs_invivo, height= 11, width = 15, dpi = 600)


heatmap_order<- similarity_matrix_new[invitro_vs_invivo$tree_row$order,invitro_vs_invivo$tree_col$order] # This reorders the rows of similarity_matrix according to the clustering order
heatmap_order$sequences <- row.names(heatmap_order) #designated to column 91
heatmap_order$serial_row_number <- 1:516 #new row number according to heatmap
#designated to column 92
heatmap_order$original_row_order<- invitro_vs_invivo$tree_row$order
#designated to column 93

heatmap_order2_invitro_vs_invivo <- heatmap_order %>% select(original_row_order, serial_row_number, sequences, 1:90)
fwrite(heatmap_order2_invitro_vs_invivo, "invitro_vs_invivo_sequences_without_gaps_similarity_matrix_heatmap_order_output.csv")



##attempted to make similarity score histogram 
heatmap_order2<- fread("invitro_vs_invivo_sequences_without_gaps_similarity_matrix_heatmap_order_output.csv", header = TRUE, sep = ",")
heatmap_order2_matrix<- heatmap_order2[,-c(1,2,3)]


heatmap_order2_matrix_scores <- as.vector(heatmap_order2_matrix) #list of 90 instead of single vector
heatmap_order2_matrix_scores_flat <- unlist(heatmap_order2_matrix)


table(heatmap_order2_matrix_scores_flat)
heatmap_order2_matrix_scores_flat
#0  5.56 11.11 16.67 22.22 27.78 33.33 38.89 44.44    50 55.56 61.11 66.67 72.22 
#348   928  3060  6080  9132  9133  7656  5307  3110  1161   414    85    22     4


scores_data<- as.data.frame(heatmap_order2_matrix_scores_flat)
colnames(scores_data)<- "identity_perc_score"
fwrite(scores_data, "invitro_vs_invivo_sequences_without_gaps_similarity_matrix_histogram_input.csv")

heatmap_order2_matrix_scores_histogram<-scores_data %>% group_by(scores_data$identity_perc_score) %>% dplyr::count()
colnames(heatmap_order2_matrix_scores_histogram)<- c("identity_perc_score", "frequency")
fwrite(heatmap_order2_matrix_scores_histogram, "invitro_vs_invivo_sequences_without_gaps_similarity_matrix_histogram_frequency_table.csv")



#if you want to see more like, midpoints and breaks then 
#check<- hist(heatmap_order2_matrix_scores), you can navigate other information

#A histogram is plotted, showing how frequently different similarity scores occur.


avg_score <- round(mean(heatmap_order2_matrix_scores_flat),2)
#[1] 27.66821
#abline: Draws a vertical line at the mean of the similarity scores to highlight the central tendency
#lwd is thickness
#lty is dashed for better visualization


invitro_vs_invivo_histo<- ggplot(scores_data, aes(x = heatmap_order2_matrix_scores_flat)) +
  geom_histogram(
    binwidth = 5, 
    fill = "skyblue",
    color = "darkblue", 
    boundary = 0 
  ) +
  geom_vline(aes(xintercept = mean(heatmap_order2_matrix_scores_flat)), color = "red", linewidth = 2, linetype = "dashed")+
  labs (
    title = "Invitro sequences (R7) vs Invivo Percent Identity Score", 
    subtitle = paste("Average percent score: ", avg_score, sep= ""),
    x = "Percent Identity Score",
    y = "Frequency")+
  
  scale_x_continuous(
    limits = c(0,100),
    breaks = seq(0, 100, 10)
  )+
  scale_y_continuous(
    limits = c(0,10000),
    breaks = seq(0, 10000, 2000)
  )+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        legend.position = "top")

ggsave("invitro_vs_invivo_sequences_percent_identity_score_histogram.tiff", 
       plot= invitro_vs_invivo_histo, height = 11, width = 12, dpi=600)




