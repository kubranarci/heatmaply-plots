library(GO.db)
## List the possible values for columns

setwd("~/Desktop/PHD/SBG_RNA-seq/2019_March_BINGO_after_tik/")

GO_Diff <- read.table(file = "MBETA/mbeta_Annotated_GO_matrix_p0.0005.txt", header = TRUE, sep = "")
GO_Diff <- GO_Diff[-1,]

#deflist <- select(GO.db, keys= as.character(GO_Diff$GO_ID), columns=c("TERM","ONTOLOGY"), keytype="GOID")

#GO_evidence <- read.csv(file= "~/Desktop/DEGs_FC/goa_human.gaf", header = FALSE, sep = "\t", skip = 12)
GO_evidence_2 <- data.frame("GOID" = GO_evidence$V5, "Evidence" = GO_evidence$V7 )

merged_GO <- merge(GO_Diff, GO_evidence_2, by= "GOID", all= FALSE)
filtered_GO <- dplyr::filter(merged_GO, grepl('EXP|IDA|IPI|IMP|IGI|IEP', Evidence))
uniqu_GO <- filtered_GO[!duplicated(filtered_GO$GOID), ]

#FILTER BY p> 0.001 for huh7
merged1= (uniqu_GO)[uniqu_GO$HALPHA < 0.001 | uniqu_GO$HBETA < 0.001 | uniqu_GO$HSALPHA < 0.001 | uniqu_GO$HSBET < 0.001 | uniqu_GO$HSOR < 0.001,] 

merged2 <- merge(merged1, deflist, by = "GOID", all= FALSE) 

write.table(uniqu_GO, file= "MBETA_cluster_bingo_0.005.txt", sep = "\t" )

#FILTER BY p> 0.001 mahlavu
merged1= (uniqu_GO)[uniqu_GO$MALPHA < 0.001 | uniqu_GO$MSALPHA < 0.001 | uniqu_GO$MSBETA < 0.001 | uniqu_GO$MBETA < 0.001 | uniqu_GO$MSOR < 0.001,] 
merged2 <- merge(merged1, deflist, by = "GOID", all= FALSE) 

write.table(merged2, file= "MV_UP_Evidence_GO_.txt", sep = " " )
