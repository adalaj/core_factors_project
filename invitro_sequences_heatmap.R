
#!/usr/bin/env Rscript

setwd("/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/Nathan")

library(tidyverse)
library(data.table)
library(rBLAST)
library(Biostrings)



# To see if ncbi blastn is there or not : Sys.which("blastn")
# This is the word document provided by bruce: CF Selex Supplementary Tables.docx
# I separated each tables and saved them in separate csv files based on their table title
# Task 1: To compare similarity within inviro sequences 

invitro<- fread("S3_R7_invitro.csv", sep = ",", header = TRUE)


#make a new data table that will have identifier and its sequence next to each other. 
#Transform invitro table 


reshaped_invitro <- data.frame(
  Identifier = invitro[seq(1, nrow(invitro), 2), 1], # Extract odd rows as identifiers
  Sequence = invitro[seq(2, nrow(invitro), 2), 1]   # Extract even rows as sequence
)

colnames(reshaped_invitro)<- c("invitro_identifier", "invitro_sequences")

reshaped_invitro$length<- nchar(reshaped_invitro$invitro_sequences)

reshaped_invitro$filename<- gsub(">", "", reshaped_invitro$invitro_identifier)

reshaped_invitro$filename<- gsub("[\\|/\\(\\)]", "-", reshaped_invitro$filename)

reshaped_invitro$filename<- sub("-$", "", reshaped_invitro$filename)

fwrite(reshaped_invitro, "reshaped_S3_R7_invitro.csv")


#test<- head(reshaped_invitro, 10)

# Define the temporary directory for FASTA files
temp_dir <- "/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/bruce_lab/project/Nathan"
#temp_dir<- tempdir()

# Initialize the overall BLAST results data frame
overall_blast <- data.table()


#make balstdb for all the sequences

for (j in 1:nrow(reshaped_invitro)){
  
  # Define subject sequence and identifier
  subject_id <- reshaped_invitro$filename[j]
  subject_sequence <- reshaped_invitro$invitro_sequences[j]
  
  # Create the subject DNAStringSet and assign the identifier
  subject_string_set <- DNAStringSet(subject_sequence)
  names(subject_string_set) <- subject_id
  
  # Create the subject FASTA file path
  subject_fasta <- file.path(temp_dir, paste(subject_id, ".fasta", sep = ""))
  writeXStringSet(subject_string_set, subject_fasta)
  
  
  makeblastdb (subject_fasta, 
               db_name = subject_fasta)
  
}

query_check<- head(reshaped_invitro, 1)

# Loop through each row for query sequence and subject sequence



for (i in 1:nrow(query_check)){
  # Define query sequence and identifier
  query_id<- query_check$filename[i]
 
  #read query using readDNAstringSet from biostrings
  query_fasta<- readDNAStringSet(file.path(temp_dir, paste(query_id, ".fasta", sep = ""))) # from that temp directory
  
  #Loop through the rows for subject sequences
  
  for (j in 1:nrow(reshaped_invitro)){
    #if (i !=j){ #would make sense if i dont want compare identical query and subject
    # Define subject sequence and identifier
    subject_id <- reshaped_invitro$filename[j]
    
    #read subject using readDNAstringSet from biostrings
    subject_fasta<- readDNAStringSet(file.path(temp_dir, paste(subject_id, ".fasta", sep = ""))) # from that temp directory
    
    
    ##load blast database using blast function 
    bl <- blast(db= paste(temp_dir,"/", subject_id,".fasta", sep = ""))
   
  
    two_seq_blast <- predict(bl,query_fasta,BLAST_args = "-word_size 4") # 
    colnames (two_seq_blast) <- c("qsequid", "ssequid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
    
    
    #put zero if there is no overlapping
    if (nrow(two_seq_blast)==0){
      no_similarity_blast <- data.frame(matrix(0, nrow = 1, ncol= 12)) # create empty blast with 0 entries
      colnames (no_similarity_blast) <- c("qsequid", "ssequid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
      ##these are default col names when bl was made
      no_similarity_blast [1,1] <- names(query_fasta) # which is i or query
      
      no_similarity_blast[1,2] <- names(subject_fasta)
      
      two_seq_blast <- rbind(two_seq_blast, no_similarity_blast)
    }
    
   
    
    ##next, was to resolve the issue why there are not 20 query for each subject, i realised that some reshaped_invitro like IGHMBP2
    ## the subject query is 111, meaning all small sections of query is checked for subject, and i only need 1 which has 100% precent identity and
    ## match the same length of input query.
    if (nrow(two_seq_blast) > 1) {
      colnames (two_seq_blast) <- c("qsequid", "ssequid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
      if (query_id == subject_id){
        two_seq_blast <- two_seq_blast %>% 
          filter(length == width(query_fasta))
        
      } else{
        two_seq_blast <- two_seq_blast %>%
          filter(pident == max(pident)) #& length == max(length)) }
    }
    
    
    #else was needed bcoz we are still getting 505 not 400, 
    #meaning some longer fragment is overlapping with shorter ones multiple times resulting in multiple matches
    #in such cases, we only want to retain the match which showed max percent similarity between given two sequence
    overall_blast <- rbind(overall_blast,two_seq_blast) 
    
    # files are automatically removed from temporary directory
  
    
    
  } 

  print(query_id)
  
  }
}

fwrite(overall_blast, paste("invitro_sequences_blast_word_count", z, sep = ""))
  


#####Issues with these results
#1) qseqid        sseqid        pident length mismatch gapopen qstart qend sstart send   evalue bitscore
#1 R7_1_8_20_2_2 R7_1_8_20_2_2  100.0     18        0       0      1   18      1   18   1.02e-08    34.4
#2 R7_1_8_20_2_2 R7_1_8_20_2_2   87.5      8        1       0     10   17      2    9   1.70e-01    10.4
#3 R7_1_8_20_2_2 R7_1_8_20_2_2   87.5      8        1       0      2    9     10    17   1.70e-01   10.4
#4 R7_1_8_20_2_2 R7_1_13_20_1_4  100.0     4        0       0      6    9      7    10  0.62        8.5     
#5 R7_1_8_20_2_2 R7_1_13_20_1_4  100.0     4        0       0     15   18      7    10  0.62       8.5 
#6 R7_1_8_20_2_2 R7_1_13_20_1_4  100.0     4        0       0      8   11      7    10  0.62       8.5 

#here we can see that, when it is one to one compariosn we get 100% identity percent value of length 18 (row1). I was expecting to see 100% only if we 18 length match. 
#However, when in case of row 4 R7_1_8_20_2_2 R7_1_13_20_1_4, we saw 100% with length 4
# so this is the issue with word count 4, it is not giving me desired result.

# in few cases sstart is bigger than send, which means it also reverses the sequence
#for more refer NCBI global alignment vs local alignment.docx 



###Task is to try global alignment instead of local alignment 
##all the fasta files, make db files that were generated in local alignment case are now deleted. 


#biostrings have pairwise alignment function that perform global alignment:

seq1<-"GTGTGAGGATGTGGGGAG"

seq2<- "GGGGAAGTGTGGAGTTGA" #(R7-2-11-20-21-5)

sigma<- pwalign::nucleotideSubstitutionMatrix(match = 2, mismatch = -3, baseOnly = TRUE)
#refer https://rpubs.com/Abhatt/835781

#this sigma have match and unmatch scores as per NCBI blast global alignment default parameters.

#match/mismatch is 2, -3
#existence is 5 (which is gap opening)
#gap extension is extension is 2)


alignment1v4<- pwalign::pairwiseAlignment(seq1, seq2,substitutionMatrix = sigma,type="global", gapOpening = -5, gapExtension = -2, degap = TRUE)
#Global PairwiseAlignmentsSingleSubject (1 of 1)
#pattern: GTG-TGAGGATGTGGGGAG
#subject: GCGGTGAGAGGTTGGGGA-
#  score: -5 

attempt1<- data.frame((matrix(nrow = 0, ncol=13)))
colnames(attempt1)<- c("query", "subject", "similarity_perc", "NW_score_NCBI", "aligned_length", "match",
                       "identity_perc_NCBI","mismatch", "gapOpening","gapExtension","gaps_perc",
                       "aligned_query","aligned_subject")


for (i in 1:nrow(reshaped_invitro)){
  # Define query sequence and identifier
  query_id<- reshaped_invitro$filename[i]
  query_seq <- reshaped_invitro$invitro_sequences[i]
  
  #Loop through the rows for subject sequences
  
  for (j in 1:nrow(reshaped_invitro)){
    #if (i !=j){ #would make sense if i dont want compare identical query and subject
    # Define subject sequence and identifier
    subject_id <- reshaped_invitro$filename[j]
    subject_seq <- reshaped_invitro$invitro_sequences[j]
    
    alignment<- pwalign::pairwiseAlignment(query_seq, subject_seq,substitutionMatrix = sigma,type="global", gapOpening = -5, gapExtension = -2)
    query<-  query_id
    subject<-  subject_id
    similarity_perc<- round(pwalign::pid(alignment),2) #similarity percent
    NW_score_NCBI<- pwalign::score(alignment) #NW score
    aligned_length<- width(pwalign::alignedPattern(alignment)) #The 2 objects returned by alignedPattern(x) and alignedSubject(x) are guaranteed to have the same shape (i.e. same length() and width())
    match<- pwalign::nmatch(alignment)
    identity_perc_NCBI<- round((match/aligned_length)*100, 2)#as per ncbi
    mismatch<- pwalign::nmismatch(alignment)
    gapOpening<- alignment@gapOpening
    gapExtension<- alignment@gapExtension
    gaps_perc<- round((gapExtension/aligned_length)*100, 2)
    aligned_query<-as.character(pwalign::alignedPattern(alignment)) #aligned query
    aligned_subject<- as.character(pwalign::alignedSubject(alignment)) #aligned subject
    
    
    
    final<- data.frame(query, subject, similarity_perc, NW_score_NCBI, aligned_length, match,
                       identity_perc_NCBI,mismatch, gapOpening,gapExtension,gaps_perc,
                       aligned_query,aligned_subject)
    
    attempt1 <- rbind(attempt1, final)
  }
  
}

fwrite(attempt1, "invitro_sequences_global_blast.csv")


globalblast_4c<- attempt1 %>% select(query, subject, identity_perc_NCBI) #4c - 4 columns
globalblast_4c$identifier <- paste(globalblast_4c$query, globalblast_4c$subject, sep = "&")

identifier_list <- list(reshaped_invitro$filename) #if list is removed it will save as values with 516 strings.
#it is a list of list

similarity_matrix <- data.frame(matrix(NA, nrow=0, ncol=516)) #additional column is for identifier column. This column will contain query name. 
colnames(similarity_matrix)<- t(reshaped_invitro[,4]) #4 is filename

test<- data.frame(matrix(NA, nrow = 0, ncol = 1))
colnames(test)<- "invitro_sequences"

similarity_matrix <- cbind(test, similarity_matrix) #colnames now increased to 517.


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
  similarity_matrix[nrow(similarity_matrix)+1,] <- my_row #rbind will not work because my myrow gives as a column 
  
} 


fwrite(similarity_matrix,
       "invitro_sequences_similarity_matrix.csv")
# if you plan to read the similarity table again then read as data as data frame 


similarity_matrix<- fread("invitro_sequences_similarity_matrix.csv", header= TRUE, sep = ",")
similarity_matrix <- as.data.frame(similarity_matrix)
similarity_matrix_new<- similarity_matrix[,-c(1)]
row.names(similarity_matrix_new)<- t(colnames(similarity_matrix_new)) #if as.data.frame is done then only this code will work.

similarity_matrix_new<- similarity_matrix_new %>% mutate(across(1:516, as.numeric))

#rownames does appear in excel as well in R. Deleted all files.
#fwrite(similarity_matrix_new,
#file= "/Users/jyotiadala/Library/CloudStorage/OneDrive-SUNYUpstateMedicalUniversity/project/Lieber_RLFS_Validation/Class_switch_recombination/Cluster_analysis/graph/input/template_RLFS_CSR_no_pseudogene_similarity_matrix_heatmap.csv", sep = ",")

#no point in saving as row names are deleted by default leaving no identifier
#whenever you want to plot heatmap perform the above steps, as in read similarity matrix, remove first column, transform column name and make them numeric
pheatmap(similarity_matrix_new, 
         main = "Comparative Analysis of Sequence Similarity in Invitro_sequences (R7) ", 
         treeheight_row = 100,treeheight_col=100, display_numbers = F, number_color= "white", fontsize_number = 6, border_color = "grey60")

invitro_sequence_heatmap<- pheatmap(similarity_matrix_new, 
                                    main = "Comparative Analysis of Sequence Similarity in Invitro_sequences (R7) ", 
                                    treeheight_row = 100,treeheight_col=100, display_numbers = F, number_color= "white", fontsize_number = 6, border_color = "grey60")


ggsave( "invitro_sequences(R7)_heatmap.png", 
        plot = invitro_sequence_heatmap, height= 30, width = 30, dpi = 300)



#The results are not accurate here because the scores are calculated by considering gaps.
#Gaps are introduced, extended, and then aligned, which increases the alignment length and impacts the false similarity score.
#Since these sequences are pre-aligned, we do not want to introduce any gaps. One way to address this is by increasing the gap-opening and gap-extension penalties to very high values.
#The high cost of introducing a new gap ensures that gaps are avoided in the alignment


attempt2<- data.frame((matrix(nrow = 0, ncol=10)))
colnames(attempt2)<- c("query", "subject", "similarity_perc", "NW_score_NCBI", "aligned_length", "match",
                       "identity_perc_NCBI","mismatch",
                       "aligned_query","aligned_subject")

for (i in 1:nrow(reshaped_invitro)){
  # Define query sequence and identifier
  query_id<- reshaped_invitro$filename[i]
  query_seq <- reshaped_invitro$invitro_sequences[i]
  
  #Loop through the rows for subject sequences
  
  for (j in 1:nrow(reshaped_invitro)){
    #if (i !=j){ #would make sense if i dont want compare identical query and subject
    # Define subject sequence and identifier
    subject_id <- reshaped_invitro$filename[j]
    subject_seq <- reshaped_invitro$invitro_sequences[j]
    
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




fwrite(attempt2, "invitro_sequences_global_blast_withoutgaps.csv")


##to make similarity heatmap: ##task 1: rotate the overall_blast in matrix form, where we have 516 col and 516 rows showing 266256 combinations

globalblast_4c<- attempt2 %>% select(query, subject, identity_perc_NCBI) #4c - 4 columns
globalblast_4c$identifier <- paste(globalblast_4c$query, globalblast_4c$subject, sep = "&")

identifier_list <- list(reshaped_invitro$filename) #if list is removed it will save as values with 516 strings.
#it is a list of list

similarity_matrix <- data.frame(matrix(NA, nrow=0, ncol=516)) #additional column is for identifier column. This column will  contain query name. 
colnames(similarity_matrix)<- t(reshaped_invitro[,4]) #4 is filename

test<- data.frame(matrix(NA, nrow = 0, ncol = 1))
colnames(test)<- "invitro_sequences"

similarity_matrix <- cbind(test, similarity_matrix) #colnames now increased to 517.


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
       "invitro_sequences_similarity_matrix_without_gaps.csv")
# if you plan to read the similarity table again then read as data as data frame 


similarity_matrix<- fread("invitro_sequences_similarity_matrix_without_gaps.csv", header= TRUE, sep = ",")
similarity_matrix <- as.data.frame(similarity_matrix)
similarity_matrix_new<- similarity_matrix[,-c(1)]
row.names(similarity_matrix_new)<- t(colnames(similarity_matrix_new)) #if as.data.frame is done then only this code will work.

similarity_matrix_new<- similarity_matrix_new %>% mutate(across(1:516, as.numeric))
#whenever you want to plot heatmap perform the above steps, as in read similarity matrix, remove first column, transform column name and make them numeric

#rownames does appear in excel as well in R. Deleted all files.
fwrite(similarity_matrix_new,
       "invitro_sequences_without_gaps_similarity_matrix_heatmap_input.csv", sep = ",")




#pheatmap(similarity_matrix_new, cluster_rows = TRUE, cluster_cols = FALSE)
# here the clustering of column is happening randomly. not sure of order. Also didnt see diagonal line





# Convert similarity matrix to distance matrix
# A similarity of 100% (perfect match) gives a distance of 0
# A similarity of 0% (no match) gives a distance of 1
distance_matrix <- 1 - (similarity_matrix_new / 100)
fwrite(distance_matrix, "invitro_sequences_without_gaps_similarity_matrix_heatmap_distance_input.csv")


# Convert distance matrix to dist object
dist_object <- as.dist(distance_matrix)

#Perform hierarchical clustering
row_clustering <- hclust(dist_object, method = "average")  # Choose method: average, complete, single, etc


#For two clusters R and S, the single linkage returns the minimum distance between two points i and j such that i belongs to R and j belongs to S.
#For two clusters R and S, the complete linkage returns the maximum distance between two points i and j such that i belongs to R and j belongs to S.
#For two clusters R and S, first for the distance between any data-point i in R and any data-point j in S and then the arithmetic mean of these distances are calculated.
#Average Linkage returns this value of the arithmetic mean.


invitro_sequence_heatmap_without_gaps<- pheatmap(
  similarity_matrix_new,                  # Original similarity matrix
  cluster_rows = row_clustering,     # Use custom clustering for rows
  cluster_cols = row_clustering,     # (Optional) Same clustering for columns
  scale = "none",                    # Don't scale data
  main = "Comparative Analysis of Sequence Similarity in Invitro_sequences (R7)"
)

invitro_sequence_heatmap_without_gaps_using_distance<- pheatmap(
   distance_matrix,   #it will be same but colors will be changed because similar sequence have short distance, least as 0               
   cluster_rows = row_clustering,     
   cluster_cols = row_clustering,     
   scale = "none",                    
  main = "Comparative Analysis of Sequence Similarity distance in Invitro_sequences (R7)")



ggsave( "invitro_sequences(R7)_without_gaps_heatmap.png", 
        plot = invitro_sequence_heatmap_without_gaps, height= 30, width = 30, dpi = 300)

ggsave( "invitro_sequences(R7)_without_gaps_heatmap_distance.png", 
        plot = invitro_sequence_heatmap_without_gaps_using_distance, height= 30, width = 30, dpi = 300)

#no point in saving as row names are deleted by default leaving no identifier

clusters <- cutree(row_clustering, k = 5)
table(clusters)

#clusters
#1   2   3   4   5 
#387  48   6  74   1 

annotation <- data.frame(Cluster = factor(clusters))
rownames(annotation) <- rownames(similarity_matrix_new)
# Add annotation to the heatmap
invitro_sequence_heatmap_without_gaps_with_clusters<- pheatmap(
  similarity_matrix_new,
  cluster_rows = row_clustering,
  cluster_cols = row_clustering,
  annotation_row = annotation,
  main = "invitro_sequence_heatmap_without_gaps_with_clusters"
)

ggsave("invitro_sequence_heatmap_without_gaps_with_clusters.png", 
       plot= invitro_sequence_heatmap_without_gaps_with_clusters, height = 30, width = 30, dpi=300)

cluster_assignments_4<- cluster_assignments %>% filter(Cluster ==4)

#This heatmap suggests that there is a large red cluster at the top of the figure, while the rest is filled with yellow and blue colors, 
#supporting the hypothesis that core factor binding is sequence-independent.
#The question is: what makes these red sequences different from the others? Could it be that they are more GC-rich?



#invitro_sequence_heatmap_without_gaps is a list of four lists, each containing seven values. 
#It has an 'order' list that contains the row order under tree_row. 
#The row order is a numeric value that provides the original row number in a serial fashion. 
#For example, if the first row number is 50, it means the identifier used in row 50 of the original similarity matrix is now the first row and column of the heatmap (output graph) 
#There is another list called 'labels' in the same tree_row, which provides the names of the rows in the original similarity matrix."



#clustered_row_order<- as.data.frame(invitro_sequence_heatmap_without_gaps$tree_row$order)
#this only provide numeric value therefore not that helpful


sorted_data<- similarity_matrix[invitro_sequence_heatmap_without_gaps$tree_row$order,] # This reorders the rows of similarity_matrix according to the clustering order

#The first column will provide the sequence identifier name based on first column of sorted data 
heatmap_order<- sorted_data %>% select("invitro_sequences") #captures first row

sorted_data4 <-similarity_matrix[invitro_sequence_heatmap_without_gaps$tree_row$order, #This reorders the rows of similarity_matrix according to the clustering order
                                 sorted_data$invitro_sequences]#reorder the columns based on invitro_sequences

heatmap_order <- cbind(heatmap_order, sorted_data4)
heatmap_order$original_row_order <- row.names(heatmap_order) #get original row number based on similarity_matrix
#designated to column 518

heatmap_order$serial_row_number <- 1:516 #new row number according to heatmap
#designated to column 519

heatmap_order2 <- heatmap_order %>% select(original_row_order, serial_row_number, invitro_sequences, 2:517)


fwrite(heatmap_order2, "invitro_sequences_without_gaps_similarity_matrix_heatmap_order_output.csv")

which(colnames(heatmap_order2)=="R7-1-18-20-27-5") # "R7-1-18-20-27-5" this is the column that has value of 44.44 next to 94.44
#70
#row: 67 identifier: R7-1-18-20-27-5

heatmap_order2_high_clustered<- heatmap_order2[c(1:67),c(1:70)]
fwrite(heatmap_order2_high_clustered, "invitro_sequences_without_gaps_similarity_matrix_heatmap_order_high_similarity_score_cluster.csv")




hc_similarity_matrix <- as.data.frame(heatmap_order2_high_clustered) #hc means high_clustered
hc_similarity_matrix_new<- hc_similarity_matrix[,-c(1,2,3)]
hc_row.names(hc_similarity_matrix_new)<- t(hc_colnames(similarity_matrix_new)) #if as.data.frame is done then only this code will work.

hc_similarity_matrix_new<- hc_similarity_matrix_new %>% mutate(across(1:67, as.numeric))


invitro_sequence_heatmap_without_gaps_high_similarity_cluster<- pheatmap(hc_similarity_matrix_new, cluster_rows = FALSE, cluster_cols = FALSE, SCALE= "none")

ggsave("invitro_sequence_heatmap_without_gaps_high_similarity_cluster.png", 
       plot= invitro_sequence_heatmap_without_gaps_high_similarity_cluster, height = 30, width = 30, dpi=300)

