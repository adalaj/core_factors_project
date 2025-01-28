#!/usr/bin/env Rscript
# Author: Jyoti Devendra Adala under supervision of Dr. Vladimir A Kuznetsov 
# For updates and contributions, visit : https://github.com/adalaj
setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/Nathan")

library(tidyverse)
library(data.table)
library(Biostrings)
library(pheatmap)

# I separated each tables and saved them in separate csv files based on their table title
# Task 1: To compare similarity between invitro and theoretical sequences 

reshaped_theoretical<- fread("reshaped_S2_Rt_invitro.csv", sep = ",", header = TRUE)
#500 rows

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
  
  for (j in 1:nrow(reshaped_theoretical)){
    #if (i !=j){ #would make sense if i dont want compare identical query and subject
    # Define subject sequence and identifier
    subject_id <- reshaped_theoretical$filename[j]
    subject_seq <- reshaped_theoretical$theoretical_sequences[j]
    
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

fwrite(attempt2, "invitro_vs_theoretical_sequences_global_blast_withoutgaps.csv")
#258000


globalblast_4c<- attempt2 %>% select(query, subject, identity_perc_NCBI) #4c - 4 columns
globalblast_4c$identifier <- paste(globalblast_4c$query, globalblast_4c$subject, sep = "&")

similarity_matrix <- data.frame(matrix(NA, nrow=0, ncol=500))

colnames(similarity_matrix)<- t(reshaped_theoretical[,4]) #4 is filename

test<- data.frame(matrix(NA, nrow = 0, ncol = 1)) #additional column is for identifier column. This column will  contain query name
colnames(test)<- c( "identifier")

similarity_matrix <- cbind(test, similarity_matrix) #colnames now increased to 91.


for (i in 1:nrow(reshaped_invitro)){
  listofvalues<- list()
  listofvalues <- append(listofvalues, reshaped_invitro$filename[i])
  for (j in 1:length(reshaped_theoretical$filename)){
    k <- paste(reshaped_invitro$filename[i], reshaped_theoretical$filename[j], sep="&") 
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
       "invitro_vs_theoretical_sequences_similarity_matrix_without_gaps.csv")

# if you plan to read the similarity table again then read as data as data frame 


similarity_matrix<- fread("invitro_vs_theoretical_sequences_similarity_matrix_without_gaps.csv", header= TRUE, sep = ",")
similarity_matrix <- as.data.frame(similarity_matrix)
similarity_matrix_new<- similarity_matrix
row.names(similarity_matrix_new)<- similarity_matrix_new$identifier
similarity_matrix_new<- similarity_matrix_new[,-c(1)]

similarity_matrix_new<- similarity_matrix_new %>% mutate(across(1:500, as.numeric))
#whenever you want to plot heatmap perform the above steps, as in read similarity matrix, remove first column, transform column name and make them numeric

#rownames does appear in excel as well in R. Deleted all files.
fwrite(similarity_matrix_new,
       "invitro_vs_theoretical_sequences_without_gaps_similarity_matrix_heatmap_input.csv", sep = ",")




invitro_vs_theoretical<- pheatmap(similarity_matrix_new, 
                             cluster_rows = TRUE, 
                             cluster_cols = TRUE, 
                             scale = "none", 
                             main="Comparative Analysis of Sequence Percent Identity between Invitro (R7) and theoretical sequences")

#dist_object <- as.dist(distance_matrix)
#Warning message:
  #In as.dist.default(distance_matrix) : non-square matrix 
# due to which
#Error in hclust(dist_object, method = "average") : dissimilarities of improper length

ggsave( "invitro_vs_theoretical_sequences_without_gaps_heatmap.tiff", 
        plot = invitro_vs_theoretical, height= 11, width = 15, dpi = 600)


heatmap_order<- similarity_matrix_new[invitro_vs_theoretical$tree_row$order,invitro_vs_theoretical$tree_col$order] # This reorders the rows of similarity_matrix according to the clustering order
heatmap_order$sequences <- row.names(heatmap_order) #designated to column 501

heatmap_order$serial_row_number <- 1:516 #new row number according to heatmap
#designated to column 502
heatmap_order$original_row_order<- invitro_vs_theoretical$tree_row$order
#designated to column 503

heatmap_order2_invitro_vs_theoretical <- heatmap_order %>% select(original_row_order, serial_row_number, sequences, 1:500)
fwrite(heatmap_order2_invitro_vs_theoretical, "invitro_vs_theoretical_sequences_without_gaps_similarity_matrix_heatmap_order_output.csv")


##attempted to make similarity score histogram 
heatmap_order2<- fread("invitro_vs_theoretical_sequences_without_gaps_similarity_matrix_heatmap_order_output.csv", header = TRUE, sep = ",")
heatmap_order2_matrix<- heatmap_order2[,-c(1,2,3)]


heatmap_order2_matrix_scores <- as.vector(heatmap_order2_matrix) #list of 500 instead of single vector
heatmap_order2_matrix_scores_flat <- unlist(heatmap_order2_matrix)


table(heatmap_order2_matrix_scores_flat)
heatmap_order2_matrix_scores_flat
#   0  5.56 11.11 16.67 22.22 27.78 33.33 38.89 44.44    50 55.56 61.11 66.67 72.22 77.78 
#1582  9231 25654 43870 54290 50712 36851 21145  9557  3684  1094   221   106     2     1


scores_data<- as.data.frame(heatmap_order2_matrix_scores_flat)
colnames(scores_data)<- "identity_perc_score"
fwrite(scores_data, "invitro_vs_theoretical_sequences_without_gaps_similarity_matrix_histogram_input.csv")

heatmap_order2_matrix_scores_histogram<-scores_data %>% group_by(scores_data$identity_perc_score) %>% dplyr::count()
colnames(heatmap_order2_matrix_scores_histogram)<- c("identity_perc_score", "frequency")
fwrite(heatmap_order2_matrix_scores_histogram, "invitro_vs_theoretical_sequences_without_gaps_similarity_matrix_histogram_frequency_table.csv")



#if you want to see more like, midpoints and breaks then 
#check<- hist(heatmap_order2_matrix_scores), you can navigate other information

#A histogram is plotted, showing how frequently different similarity scores occur.


avg_score <- round(mean(heatmap_order2_matrix_scores_flat),2)
#[1] 24.89853
#abline: Draws a vertical line at the mean of the similarity scores to highlight the central tendency
#lwd is thickness
#lty is dashed for better visualization


invitro_vs_theoretical_histo<- ggplot(scores_data, aes(x = heatmap_order2_matrix_scores_flat)) +
  geom_histogram(
    binwidth = 5, 
    fill = "skyblue",
    color = "darkblue", 
    boundary = 0 
  ) +
  geom_vline(aes(xintercept = mean(heatmap_order2_matrix_scores_flat)), color = "red", linewidth = 2, linetype = "dashed")+
  labs (
    title = "Invitro sequences (R7) vs Theoretical Percent Identity Score", 
    subtitle = paste("Average percent score: ", avg_score, sep= ""),
    x = "Percent Identity Score",
    y = "Frequency")+
  
  scale_x_continuous(
    limits = c(0,100),
    breaks = seq(0, 100, 10)
  )+
  scale_y_continuous(
    limits = c(0,55000),
    breaks = seq(0, 55000, 10000)
  )+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        legend.position = "top")

ggsave("invitro_vs_theoretical_sequences_percent_identity_score_histogram.tiff", 
       plot= invitro_vs_theoretical_histo, height = 11, width = 12, dpi=600)





