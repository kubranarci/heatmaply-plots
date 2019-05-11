# Generate Node/edge Matrix for heatmap analysis
install.packages('heatmaply')

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

library(stringr)
# display.brewer.pal(11, "BrBG")
BrBG <- colorRampPalette(brewer.pal(11, "BrBG"))
Spectral <- colorRampPalette(brewer.pal(11, "Spectral"))


setwd("~/Desktop/PHD/SBG_RNA-seq/9.2-BINGO_after_tik/")

data.path <- paste0(getwd(), "HALPHA/")
count.files <- list.files(data.path, ".bgo", recursive = T, full.names = T)
print(count.files)
snames <- gsub(".bgo", "", basename(count.files))
print(snames)

merged = data.frame("GOID"="")
n <- ncol(merged)

for (countf in count.files) {
  temp <- read.delim(countf, header = TRUE, stringsAsFactors = F)
  name <- gsub("_all.bgo", "", basename(countf))
  
  temp <- temp[temp$pvalue < 0.0005, ]
  temp=data.frame("GO_ID"= temp$GOID, temp$pvalue, stringsAsFactors = F )
  
  temp$GO_ID <- str_pad(temp$GO_ID, 7, pad = "0")
  temp$GO_ID <- paste0("GO:", temp$GO_ID)
  
    merged=merge(merged,temp, by="GO_ID", all=TRUE)
    n=n+1
    colnames(merged)[n] <-  paste("", name, sep="")    
  

}

rownames(merged) <- merged$GO_ID
merged <- merged[,-1]
merged[is.na(merged)] <- 0.0005
merged[merged != 0] <- 1
merged <- merged[-1,]
merged2 <- data.frame(merged, TOTAL = merged$HALPHA +merged$HBETA +merged$HSALPHA + merged$HSBETA + merged$HSOR + merged$MALPHA + merged$MBETA + merged$MSALPHA + merged$MSBETA + merged$MSOR
                      , stringsAsFactors = F)

write.table(merged, file= "ALL_GO_matrix_p0.0005.txt", sep = " " )
merged3 <- merged[-2]
merged3 <- merged3[-6]
heatmaply(merged,
          margins = c(80, 30, 30, 1), 
          dedogram= "raw",
          xlab = "samples"  ,ylab= "GO IDs",
          main = "GOs from PCST Networks",col = cool_warm, labCol= colnames(merged),
          k_col = 1, k_row = 1, grid_gap = 0.5,
          showticklabels = c(TRUE, FALSE)
)


#for html file
heatmaply(merged,
          margins = c(120, 120, 30, 1), 
          dedogram= "raw",
          xlab = "samples"  ,ylab= "GO IDs",
          main ="GOs from PCST Networks",col = cool_warm,labCol= colnames(merged),
          k_col = 1, k_row = 1, grid_gap = 0.5,
          showticklabels = c(TRUE, TRUE)
)

