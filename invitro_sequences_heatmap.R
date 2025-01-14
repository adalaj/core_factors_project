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