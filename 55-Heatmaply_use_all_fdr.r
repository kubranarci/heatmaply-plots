#plots heatmap using heatmaply and reports number of up/down regulated genes and the number of untranslated transcripts


# You'll need devtools
install.packages.2 <- function (pkg) if (!require(pkg)) install.packages(pkg);
install.packages.2('devtools')
# make sure you have Rtools installed first! if not, then run:
#install.packages('installr'); install.Rtools()

devtools::install_github("ropensci/plotly") # you will probably benefit from the latest version of plotly
devtools::install_github('talgalili/heatmaply')
library('heatmaply')

setwd("~/Desktop/PHD/SBG_RNA-seq/15-heatmaply-matching_to_all_edger_after_tik")
data.path <- paste0(getwd(), "/NO_FİLT/MV")
count.files <- list.files(data.path, ".txt$", recursive = T, full.names = T)
print(count.files)
snames <- gsub(".txt", "", basename(count.files))
print(snames)

merged = data.frame(ensembleid="")
n <- ncol(merged)

for (countf in count.files) {
  temp <- read.delim(countf, header = TRUE,sep = ",", stringsAsFactors = F)
  temp <- temp[temp$logFC != 0, ]
  name <- gsub(".txt", "", basename(countf))
  temp=data.frame("ensembleid"= temp$ensembleid, temp$logFC, stringsAsFactors = F )
  merged=merge(merged,temp, by="ensembleid", all=TRUE)
  n=n+1
  colnames(merged)[n] <-  paste("", name, sep="") 
}

rownames(merged) <- merged$ensembleid
merged <- merged[-1,]
merged[is.na(merged)] <- 0
merged <- data.frame(merged[,-1])

#FILTER BY -/+2 for huh7
merged1= (merged)[merged$HALPHA >= 2 | merged$HBETA >= 2 | merged$HSALPHA >= 2 | merged$HSBETA >= 2 | merged$HSOR >= 2 | merged$HALPHA <= -2 | merged$HBETA <= -2 | merged$HSALPHA <= -2  | merged$HSBETA <= -2  | merged$HSOR <= -2  | merged$HALPHA <= -2, ] 
write.table(merged, file= "ALL_DEGS_huh7.txt", sep = " " )
write.table(merged1, file= "ALL_DEGS_huh7_filtered.txt", sep = " " )

#FILTER BY -/+2 mahlavu
merged1= (merged)[merged$MALPHA >= 2 | merged$MBETA >= 2 | merged$MSALPHA >= 2 | merged$MSBETA >= 2 | merged$MSOR >= 2 | merged$MALPHA <= -2 | merged$MBETA <= -2 | merged$MSALPHA <= -2 | merged$MSBETA <= -2  | merged$MSOR <= -2  | merged$MALPHA <= -2, ] 
write.table(merged, file= "ALL_DEGS_MAHLAVU.txt", sep = " " )
write.table(merged1, file= "ALL_DEGS_MAHLAVU_filtered.txt", sep = " " )



#for html file
heatmaply(merged1,
          margins = c(10, 400, 1, 1), 
          dedogram= "raw",
          scale_fill_gradient_fun = scale_fill_gradient2(low = "green", high = "red", midpoint = 0, limits = c(-5, 5)),
          k_col = 3, k_row = 2, grid_gap = 0.1,
          fontsize_row = 5,dendrogram = TRUE,
          fontsize_col = 8,
          showticklabels = c(TRUE, FALSE),
          file= "MV_degs_pic.html"
)

