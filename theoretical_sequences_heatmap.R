#!/usr/bin/env Rscript

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/Nathan")

library(tidyverse)
library(data.table)
library(Biostrings)
library(pheatmap)




# To see if ncbi blastn is there or not : Sys.which("blastn")
# This is the word document provided by bruce: CF Selex Supplementary Tables.docx
# I separated each tables and saved them in separate csv files based on their table title
# Task 1: To compare similarity within theoretical sequences 

theoretical<- fread("S2_Rt_invitro.csv", sep = ",", header = TRUE)


#make a new data table that will have identifier and its sequence next to each other. 
#Transform theoretical table 


reshaped_theoretical <- data.frame(
  Identifier = theoretical[seq(1, nrow(theoretical), 2), 1], # Extract odd rows as identifiers
  Sequence = theoretical[seq(2, nrow(theoretical), 2), 1]   # Extract even rows as sequence
)

colnames(reshaped_theoretical)<- c("theoretical_identifier", "theoretical_sequences")

reshaped_theoretical$length<- nchar(reshaped_theoretical$theoretical_sequences)

reshaped_theoretical$filename<- gsub(">", "", reshaped_theoretical$theoretical_identifier)

reshaped_theoretical$filename<- gsub("[\\|/\\(\\)]", "-", reshaped_theoretical$filename)

reshaped_theoretical$filename<- sub("-$", "", reshaped_theoretical$filename)

#here sequence are in lowercase, so making it uppercase

reshaped_theoretical$theoretical_sequences <- toupper(reshaped_theoretical$theoretical_sequences)

anyDuplicated(reshaped_theoretical$filename)
#0

fwrite(reshaped_theoretical, "reshaped_S2_Rt_invitro.csv")

sigma<- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE)
attempt2<- data.frame((matrix(nrow = 0, ncol=10)))
colnames(attempt2)<- c("query", "subject", "similarity_perc", "NW_score_NCBI", "aligned_length", "match",
                       "identity_perc_NCBI","mismatch",
                       "aligned_query","aligned_subject")

for (i in 1:nrow(reshaped_theoretical)){
  # Define query sequence and identifier
  query_id<- reshaped_theoretical$filename[i]
  query_seq <- reshaped_theoretical$theoretical_sequences[i]
  
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




fwrite(attempt2, "theoretical_sequences_global_blast_withoutgaps.csv")
#250000

##to make similarity heatmap: ##task 1: rotate the global_blast in matrix form, where we have 516 col and 516 rows showing 266256 combinations

globalblast_4c<- attempt2 %>% select(query, subject, identity_perc_NCBI) #4c - 4 columns
globalblast_4c$identifier <- paste(globalblast_4c$query, globalblast_4c$subject, sep = "&")

identifier_list <- list(reshaped_theoretical$filename) #if list is removed it will save as values with 516 strings.
#it is a list of list

similarity_matrix <- data.frame(matrix(NA, nrow=0, ncol=500))

colnames(similarity_matrix)<- t(reshaped_theoretical[,4]) #4 is filename

test<- data.frame(matrix(NA, nrow = 0, ncol = 1)) #additional column is for identifier column. This column will  contain query name
colnames(test)<- "theoretical_sequences"

similarity_matrix <- cbind(test, similarity_matrix) #colnames now increased to 91.


for (i in 1: length (identifier_list[[1]])){
  listofvalues<- list()
  listofvalues <- append(listofvalues, identifier_list[[1]][i])
  for (j in 1: length(identifier_list[[1]])){
    k <- paste(identifier_list[[1]][i], identifier_list[[1]][j], sep="&") 
    match<- globalblast_4c[globalblast_4c$identifier== k,] #match is a dataframe that only contain one row which provided condition
    listofvalues <- append(listofvalues, match[1,3])
    #if (nrow(match)<1){ #code debug
    #print(k)
    #print(match)}#append function works by adding things row wise 
    #listofvalues [[i]][i]<- match[1,3] # this didnt work but good to know how to add data points to first entry of list
  }
  
  my_row <- do.call(rbind, listofvalues) # do.call is used to execute a function with a list of arguments.
  similarity_matrix[nrow(similarity_matrix)+1,] <- my_row #rbind will not work becuase my myrow gives as column 
  
} 


fwrite(similarity_matrix,
       "theoretical_sequences_similarity_matrix_without_gaps.csv")

# if you plan to read the similarity table again then read as data as data frame 


similarity_matrix<- fread("theoretical_sequences_similarity_matrix_without_gaps.csv", header= TRUE, sep = ",")
similarity_matrix <- as.data.frame(similarity_matrix)
similarity_matrix_new<- similarity_matrix[,-c(1)]
row.names(similarity_matrix_new)<- t(colnames(similarity_matrix_new)) #if as.data.frame is done then only this code will work.

similarity_matrix_new<- similarity_matrix_new %>% mutate(across(1:500, as.numeric))
#whenever you want to plot heatmap perform the above steps, as in read similarity matrix, remove first column, transform column name and make them numeric

#rownames does appear in excel as well in R. Deleted all files.
fwrite(similarity_matrix_new,
       "theoretical_sequences_without_gaps_similarity_matrix_heatmap_input.csv", sep = ",")




#pheatmap(similarity_matrix_new, cluster_rows = TRUE, cluster_cols = FALSE)
# here the clustering of column is happening randomly. not sure of order. Also didnt see diagonal line





# Convert similarity matrix to distance matrix
# A similarity of 100% (perfect match) gives a distance of 0
# A similarity of 0% (no match) gives a distance of 1
distance_matrix <- 1 - (similarity_matrix_new / 100)
fwrite(distance_matrix, "theoretical_sequences_without_gaps_similarity_matrix_heatmap_distance_input.csv")


# Convert distance matrix to dist object
dist_object <- as.dist(distance_matrix)

#Perform hierarchical clustering
row_clustering <- hclust(dist_object, method = "average")  # Choose method: average, complete, single, etc


#For two clusters R and S, the single linkage returns the minimum distance between two points i and j such that i belongs to R and j belongs to S.
#For two clusters R and S, the complete linkage returns the maximum distance between two points i and j such that i belongs to R and j belongs to S.
#For two clusters R and S, first for the distance between any data-point i in R and any data-point j in S and then the arithmetic mean of these distances are calculated.
#Average Linkage returns this value of the arithmetic mean.


theoretical_sequence_heatmap_without_gaps<- pheatmap(
  similarity_matrix_new,                  # Original similarity matrix
  cluster_rows = row_clustering,     # Use custom clustering for rows
  cluster_cols = row_clustering,     # (Optional) Same clustering for columns
  scale = "none",                    # Don't scale data
  main = "Comparative Analysis of Sequence Percent Identity in theoretical_sequences"
)

theoretical_sequence_heatmap_without_gaps_using_distance<- pheatmap(
  distance_matrix,   #it will be same but colors will be changed because similar sequence have short distance, least as 0               
  cluster_rows = row_clustering,     
  cluster_cols = row_clustering,     
  scale = "none",                    
  main = "Comparative Analysis of Sequence Similarity distance in Theoretical_sequences")



ggsave( "theoretical_sequences_without_gaps_heatmap.tiff", 
        plot = theoretical_sequence_heatmap_without_gaps, height= 11, width = 15, dpi = 300)

ggsave( "theoretical_sequences_without_gaps_heatmap_distance.png", 
        plot = theoretical_sequence_heatmap_without_gaps_using_distance, height= 30, width = 30, dpi = 300)

sorted_data<- similarity_matrix[theoretical_sequence_heatmap_without_gaps$tree_row$order,] # This reorders the rows of similarity_matrix according to the clustering order

#The first column will provide the sequence identifier name based on first column of sorted data 
heatmap_order<- sorted_data %>% select("theoretical_sequences") #captures first row

sorted_data4 <-similarity_matrix[theoretical_sequence_heatmap_without_gaps$tree_row$order, #This reorders the rows of similarity_matrix according to the clustering order
                                 sorted_data$theoretical_sequences]#reorder the columns based on theoretical_sequences

heatmap_order <- cbind(heatmap_order, sorted_data4)
heatmap_order$original_row_order <- row.names(heatmap_order) #get original row number based on similarity_matrix
#designated to column 502

heatmap_order$serial_row_number <- 1:500 #new row number according to heatmap
#designated to column 503

heatmap_order2_theoretical <- heatmap_order %>% select(original_row_order, serial_row_number, theoretical_sequences, 1:500)
fwrite(heatmap_order2_theoretical, "theoretical_sequences_without_gaps_similarity_matrix_heatmap_order_output.csv")




##attempted to make similarity score histogram 
heatmap_order2<- fread("theoretical_sequences_without_gaps_similarity_matrix_heatmap_order_output.csv", header = TRUE, sep = ",")
heatmap_order2<- as.data.frame(heatmap_order2)
heatmap_order2_matrix<- heatmap_order2[,-c(1,2,3)]
heatmap_order2_matrix_scores <- heatmap_order2_matrix[lower.tri(heatmap_order2_matrix, diag = FALSE)] 
#124750

#the number of lower triangular part of a square matrix, excluding teh diagonal is calculated by uisng: n(n-1)/2 where n is matrix dimension which is 516.
#lower.tri: The elements below the diagonal are TRUE, and the rest (diagonal and above) are FALSE.
#Using this logical matrix, we extract only the lower triangular elements of the similarity matrix, excluding the diagonal.
# Since similarity matrices are symmetric, the upper and lower triangular parts are identical. We only need one of them to avoid redundancy.



scores_data<- as.data.frame(heatmap_order2_matrix_scores)
colnames(scores_data)<- c("identity_perc_score")

fwrite(scores_data, "theoretical_sequences_without_gaps_similarity_matrix_histogram_input.csv")

table(scores_data)
#identity_perc_score
#0  5.56 11.11 16.67 22.22 27.78 33.33 38.89 44.44    50 55.56 61.11 66.67 72.22 
#681  4273 11984 21130 26752 25068 17703 10215  4596  1682   529   112    18     7 

heatmap_order2_matrix_scores_histogram<-scores_data %>% group_by(scores_data$identity_perc_score) %>% dplyr::count()
colnames(heatmap_order2_matrix_scores_histogram)<- c("theoretical_identity_perc_score", "frequency")
fwrite(heatmap_order2_matrix_scores_histogram, "theoretical_sequences_without_gaps_similarity_matrix_histogram_frequency_table.csv")



#if you want to see more like, midpoints and breaks then 
#check<- hist(heatmap_order2_matrix_scores), you can navigate other information

#A histogram is plotted, showing how frequently different similarity scores occur.


avg_score <- round(mean(heatmap_order2_matrix_scores),2)
#[1] 24.96
#abline: Draws a vertical line at the mean of the similarity scores to highlight the central tendency
#lwd is thickness
#lty is dashed for better visualization

#in excel, if you take entire matrix and do average =AVERAGE(D2:SY517) it gives = 36.45997484


theoretical_histo<- ggplot(scores_data, aes(x = heatmap_order2_matrix_scores)) +
  geom_histogram(
    binwidth = 5, 
    fill = "skyblue",
    color = "darkblue", 
    boundary = 0 
  ) +
  geom_vline(aes(xintercept = mean(heatmap_order2_matrix_scores)), color = "red", linewidth = 2, linetype = "dashed")+
  labs (
    title = "Theoretical sequences Percent Identity Score", 
    subtitle = paste("Average percent score: ", avg_score, sep= ""),
    x = "Percent Identity Score",
    y = "Frequency")+
  
  scale_x_continuous(
    limits = c(0,100),
    breaks = seq(0, 100, 10)
  )+
  scale_y_continuous(
    limits = c(0,26800),
    breaks = seq(0, 26800, 5000)
  )+
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        legend.position = "top")

ggsave("theoretical_sequences_percent_identity_score_histogram.tiff", 
       plot= theoretical_histo, height = 11, width = 12, dpi=600)







