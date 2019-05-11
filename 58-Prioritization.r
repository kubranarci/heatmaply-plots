
setwd("~/Desktop/PHD/SBG_RNA-seq")
data.path1 <- paste0(getwd(), "/12-2-100Randomization-aftertik/random_nodeattributes/MV")
randomization.files <- list.files(data.path1, ".tsv$", recursive = T, full.names = T)
print(randomization.files)
data.path2 <- paste0(getwd(), "/13-2-centrality_measures_after_tik/MV")
centrality.files <- list.files(data.path2, ".txt$", recursive = T, full.names = T)
print(centrality.files)

#merge_final = data.frame(Node="", Prize="", BetwennessCentrality="",FractionOfOptimalForestsContaining="", TerminalType="", Degree.Centrality="",Betwenness.Centrality="",EigenVector.Centrality="", Set="")
merge_final =list()
#n <- ncol(merge_final)

for (random in randomization.files)  {
  random_temp <- read.delim(random, header = TRUE, stringsAsFactors = F, sep = "\t")
  random_name <- gsub("_randomTerminals_nodeattributes.tsv", "", basename(random))
  random_temp <- data.frame("Node"=random_temp$Protein, random_temp[-1], stringsAsFactors = F )
  
  for (central in centrality.files){
    central_temp <- read.delim(central, header = TRUE, stringsAsFactors = F, sep = "\t")
    central_name <- gsub("_centrality.txt", "", basename(central))
    central_temp <- data.frame("Node"=central_temp$Node.Name, "Set"=central_name, central_temp[-1], stringsAsFactors = F )
   
     if (random_name == central_name){
      merged_set <- merge(random_temp, central_temp, by = "Node") 
      merged_set <- (merged_set)[merged_set$FractionOfOptimalForestsContaining == 0.01,] 
      merged_set <- (merged_set)[merged_set$BetweennessCentrality < 0.0005,] 
      merged_set <- (merged_set)[merged_set$Degree.Centrality > 0.008,] 
      merged_set <- (merged_set)[merged_set$Betweenness.Centrality > 0.008,] 
     merged_set <- (merged_set)[merged_set$EigenVector.Centrality > 0.008,] 
      merge_final <- rbind(merge_final, merged_set, stringsAsFactors = F)
     }
  }
}

for(i in 1:length(merge_final$TerminalType)) {
  if(merge_final$TerminalType[i] == "") {
    print ("yakaladık")    
    merge_final$TerminalType[i]  <- "Steiner"
  }
  else {print ("yakalayamadık")
    merge_final$TerminalType[i]  <- "Proteomic"
    
  }
} 
data <- data.frame("Inhibitor Type" = merge_final$Set, Node = merge_final$Node, Random_Betweenness = merge_final$BetweennessCentrality, Degree = merge_final$Degree.Centrality, Eigen = merge_final$EigenVector.Centrality, "Terminal Type"= merge_final$TerminalType, HCC= "Mahlavu")
#library(tidyr)   


ggplot(data = data) + 
  geom_point(aes(sort(Eigen),Node,color= Set, shape= Steiner)) +
  theme(text = element_text(size=9))+
  ggtitle("Prioritized Nodes") +
  theme(axis.text.x = element_text(face="bold", size=7, angle=45))

#huh7
data.path1 <- paste0(getwd(), "/12-2-100Randomization-aftertik/random_nodeattributes/Huh7")
randomization.files <- list.files(data.path1, ".tsv$", recursive = T, full.names = T)
print(randomization.files)
data.path2 <- paste0(getwd(), "/13-2-centrality_measures_after_tik/Huh7")
centrality.files <- list.files(data.path2, ".txt$", recursive = T, full.names = T)
print(centrality.files)

#merge_final = data.frame(Node="", Prize="", BetwennessCentrality="",FractionOfOptimalForestsContaining="", TerminalType="", Degree.Centrality="",Betwenness.Centrality="",EigenVector.Centrality="", Set="")
merge_final =list()
#n <- ncol(merge_final)

for (random in randomization.files)  {
  random_temp <- read.delim(random, header = TRUE, stringsAsFactors = F, sep = "\t")
  random_name <- gsub("_randomTerminals_nodeattributes.tsv", "", basename(random))
  random_temp <- data.frame("Node"=random_temp$Protein, random_temp[-1], stringsAsFactors = F )
  
  for (central in centrality.files){
    central_temp <- read.delim(central, header = TRUE, stringsAsFactors = F, sep = "\t")
    central_name <- gsub("_centrality.txt", "", basename(central))
    central_temp <- data.frame("Node"=central_temp$Node.Name, "Set"=central_name, central_temp[-1], stringsAsFactors = F )
    
    if (random_name == central_name){
      merged_set <- merge(random_temp, central_temp, by = "Node") 
      merged_set <- (merged_set)[merged_set$FractionOfOptimalForestsContaining == 0.01,] 
      merged_set <- (merged_set)[merged_set$BetweennessCentrality < 0.0005,] 
      merged_set <- (merged_set)[merged_set$Degree.Centrality > 0.008,] 
      merged_set <- (merged_set)[merged_set$Betweenness.Centrality > 0.008,] 
      merged_set <- (merged_set)[merged_set$EigenVector.Centrality > 0.008,] 
      merge_final <- rbind(merge_final, merged_set, stringsAsFactors = F)
    }
  }
}
for(i in 1:length(merge_final$TerminalType)) {
  if(merge_final$TerminalType[i] == "") {
  print ("yakaladık")    
    merge_final$TerminalType[i]  <- "Steiner"
    }
  else {print ("yakalayamadık")
    merge_final$TerminalType[i]  <- "Proteomic"
    
    }
}  
data2 <- data.frame("Inhibitor Type" = merge_final$Set, Node = merge_final$Node, Random_Betweenness = merge_final$BetweennessCentrality, Degree = merge_final$Degree.Centrality, Eigen = merge_final$EigenVector.Centrality, "Terminal Type"= merge_final$TerminalType, HCC= "Huh7")

  
library(tidyr)   
install.packages("tidyverse")
library(tidyverse)
data3 <- rbind(data, data2, stringsAsFactors = F)

data3$Node = with(data3, reorder(Node,Random_Betweenness))

  ggplot(data = data3) + 
    geom_point(aes(Random_Betweenness,Node,color= Inhibitor.Type, shape= Terminal.Type)) +
    facet_grid(.~HCC) +
    xlab("Betweenness Centralities of Random Networks") + ylab("Gene Nodes") +
    theme(text = element_text(size=10), axis.text.x = element_text( size=7, angle=45)) 

