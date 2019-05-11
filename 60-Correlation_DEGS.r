install.packages("ggpubr")
library("ggpubr")


ALL_DEGS_huh7 <- read.csv(file = "~/Desktop/PHD/SBG_RNA-seq/15-heatmaply-matching_to_all_edger_after_tik/ALL_DEGS_huh7.txt", header = TRUE, sep = "")

corr_Huh7<- cor(ALL_DEGS_huh7[-1], ALL_DEGS_huh7[-1], method = c("pearson", "kendall", "spearman"))
corr_MV<- cor(ALL_DEGS_MAHLAVU[-1], ALL_DEGS_MAHLAVU[-1], method = c("pearson", "kendall", "spearman"))
##cor(ALL_DEGS_huh7[-1], ALL_DEGS_MAHLAVU[-1], method = c("pearson", "kendall", "spearman"))

install.packages("corrplot")
library(corrplot)

corr<-cor(ALL_DEGS_huh7[-1])

corrplot(corr, method="color", type = "upper")

corr2 <- cor(ALL_DEGS_MAHLAVU[-1])

corrplot(corr2, method="color", type = "lower")


ggscatter(ALL_DEGS_huh7[-1], x = "HALPHA.", y = c("HSALPHA", "HBETA", "HSBETA", "HSOR") , combine= TRUE,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "Samples", ylab = "HALPHA") 

ggscatter(ALL_DEGS_huh7[-1], x = "HSALPHA", y = c("HALPHA.", "HBETA", "HSBETA", "HSOR") , combine= TRUE,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson", 
          xlab = "Samples", ylab = "HSALPHA") 

ggscatter(ALL_DEGS_huh7[-1], x = "HBETA", y = c("HALPHA.","HSALPHA", "HSBETA", "HSOR"), 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",add(new=TRUE),
          xlab = "Samples", ylab = "HBETA")

ggscatter(ALL_DEGS_huh7[-1], x = "HSBETA", y = c("HALPHA.","HSALPHA", "HBETA", "HSOR"), 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",add(new=TRUE),
          xlab = "Samples", ylab = "HSBETA")

ggscatter(ALL_DEGS_huh7[-1], x = "HSOR", y = c("HALPHA.","HSALPHA", "HBETA", "HSBETA"), 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",add(new=TRUE),
          xlab = "Samples", ylab = "HSOR")

ggscatter(ALL_DEGS_MAHLAVU[-1], x = "MALPHA", y = c("MSALPHA", "MBETA","MSBETA", "MSOR"), 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",add(new=TRUE),
          xlab = "Samples", ylab = "MALPHA")

ggscatter(ALL_DEGS_MAHLAVU[-1], x = "MBETA", y = c("MALPHA", "MSALPHA","MSBETA", "MSOR"), 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",add(new=TRUE),
          xlab = "Samples", ylab = "MBETA")

ggscatter(ALL_DEGS_MAHLAVU[-1], x = "MSALPHA", y = c("MALPHA", "MBETA","MSBETA", "MSOR"), 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",add(new=TRUE),
          xlab = "Samples", ylab = "MSALPHA")

ggscatter(ALL_DEGS_MAHLAVU[-1], x = "MSBETA", y =  c("MALPHA", "MSALPHA","MBETA", "MSOR"), 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",add(new=TRUE),
          xlab = "Samples", ylab = "MSBETA")

ggscatter(ALL_DEGS_MAHLAVU[-1], x = "MSOR", y = c("MALPHA", "MSALPHA","MBETA", "MSBETA"), 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",add(new=TRUE),
          xlab = "Samples", ylab = "MSOR")

