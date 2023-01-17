setwd("C:/Users/andre/OneDrive/Documents/UiO/Other_projects/KetamineProject/DataCollection/QuestionnaireData")
library(lme4)
library(readxl)
library(nlme)
library(ggplot2)
library(dplyr)
library(broom)
library(ggpubr)
library(lsr)
library(effsize)
library(ggsignif)
library(gridExtra)
library(ggthemes)
library(viridis)
library(devtools)

df <- read_excel("QuestionnaireScores.xlsx")

df1 <- data.frame(t(df))
colnames(df1) <- df1[1,]
df1 <- df1[-c(1),]
attach(df1)


# Loop converts all chr to num 
for (y in 1:ncol(df1)){
  if (is.character(df1[,y]) && any(is.na(df1[,y])) == FALSE){
    if (any(grepl("^[0-9]+$", df1[,y])) == TRUE){
      df1[,y] <- sapply(df1[,y], as.numeric)
    }
  }
}


lm1 <- lm(EDI_ED_Avg ~ ST_GN_Avg + SET_Exp_Avg) 

summary(lm1)

str(df1)

print(cor.test(df1$EDI_ED_Avg, df1$SET_Exp_Avg, method = "pearson"))




