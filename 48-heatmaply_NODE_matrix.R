# Generate Node/edge Matrix for heatmap analysis
install.packages("assertthat","gridExtra",
                 'devtools',"fpc","TSP","registry","gclus","gplots","RColorBrewer",
                 "stringr","labeling","yaml")

# You'll need devtools
install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
install.packages.2('devtools')
# make sure you have Rtools installed first! if not, then run:
#install.packages('installr'); install.Rtools()

devtools::install_github("ropensci/plotly") # you will probably benefit from the latest version of plotly
devtools::install_github('talgalili/heatmaply')
library('heatmaply')
library(RColorBrewer)
# display.brewer.pal(11, "BrBG")
BrBG <- colorRampPalette(brewer.pal(11, "BrBG"))
Spectral <- colorRampPalette(brewer.pal(11, "Spectral"))


setwd("~/Desktop/PHD/SBG_RNA-seq/5.2-HEATMAPs_PCST_trees_after_tik")

data.path <- paste0(getwd(), "/Nodes/MV")
count.files <- list.files(data.path, ".tsv$", recursive = T, full.names = T)
print(count.files)
snames <- gsub("_nodeattributes.tsv", "", basename(count.files))
print(snames)

merged = data.frame("Protein"="")
n <- ncol(merged)

for (countf in count.files) {
  temp <- read.delim(countf, header = TRUE, stringsAsFactors = F)
  name <- gsub("_nodeattributes.tsv", "", basename(countf))
  temp=data.frame("Protein"= temp$Protein, temp$FractionOfOptimalForestsContaining, stringsAsFactors = F )
  merged=merge(merged,temp, by="Protein", all=TRUE)
  n=n+1
  colnames(merged)[n] <-  paste("", name, sep="") 
}

rownames(merged) <- merged$Protein
merged <- merged[,-1]
merged[is.na(merged)] <- 0
#merged[merged != 0] <- 1
merged <- merged[-1,]
#merged2 <- data.frame(merged, TOTAL = merged$HALPHA +merged$HBETA +merged$HSALPHA + merged$HSBETA + merged$HSOR + merged$MALPHA + merged$MBETA + merged$MSALPHA + merged$MSBETA + merged$MSOR
#                      , stringsAsFactors = F)

#merged= (merged)[merged$HALPHA != 0 | merged$HBETA != 0 | merged$HSALPHA != 0 | merged$HSBETA != 0 | merged$HSOR != 0 |
#                   merged$MALPHA != 0 | merged$MBETA != 0 | merged$MSALPHA != 0 | merged$MSBETA != 0 | merged$MSOR != 0, ] 


write.table(merged, file= "TOTAL_node_matrixMV.txt", sep = " " )

merged3 <- merged[-2]
merged3 <- merged3[-6]
merged3 <- merged3[-8] 

heatmaply(merged,
          margins = c(80, 30, 30, 1), 
          dedogram= "raw",
          xlab = "samples"  ,ylab= "Network Nodes",
          main = "Mahlavu All Nodes from PCST Networks",
          scale_fill_gradient_fun = scale_fill_gradient2(low = "white", high = "black"),
          k_col = 3,  grid_gap = 0.1,
          showticklabels = c(TRUE, TRUE),
          file= "MV_nodes_pic_forsupport.html"
)

#for html file
heatmaply(merged,
          margins = c(10, 400, 1, 1), 
          dedogram= "raw",
          scale_fill_gradient_fun = scale_fill_gradient2(low = "white", high = "black"),
          k_col = 3, grid_gap = 0.1,ylab= "Network Nodes",
          fontsize_row = 5,dendrogram = TRUE,
          fontsize_col = 8,
          showticklabels = c(TRUE, FALSE),
          file= "MV_nodes_pic.html"
)

