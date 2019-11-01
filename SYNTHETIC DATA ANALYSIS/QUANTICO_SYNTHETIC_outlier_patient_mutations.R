rm(list=ls())
library(rmarkdown)
library(knitr)
library(gdata)
library(readxl)
library(WriteXLS)
setwd("D:/QUANTICO/QUANTICO REPRODUCIBLE CODES/SYNTHETIC DATA ANALYSIS")

#################################
### READ TCR Merged Data
#################################


SYNTHETIC_Mutation  <- read.csv('SYNTHETIC_Mutation.csv', header=FALSE)

SYNTHETIC_Mutation_mean <- colMeans(SYNTHETIC_Mutation)

##############

########### Outlier plots 


patient_4 <- SYNTHETIC_Mutation[4,]
patient_40 <- SYNTHETIC_Mutation[40,]
patient_45 <- SYNTHETIC_Mutation[45,]
patient_51 <- SYNTHETIC_Mutation[51,]
patient_56 <- SYNTHETIC_Mutation[56,]
patient_76 <- SYNTHETIC_Mutation[76,]

## T.prod.cells = col 11, TTN = 1, MUC16 = 2, RYR2 = 3, CSMD3 = 4, TP53 = 5, USH2A = 6, TMB = 7



which_var <- 2
var.name <- c("MUT.MUC16")
colors_here <- c("darkgray","darkolivegreen3", "deeppink1","blue1", "chartreuse3", "black", "red")
values <- as.numeric(c(SYNTHETIC_Mutation_mean[which_var],patient_4[which_var],patient_40[which_var],patient_45[which_var],
            patient_51[which_var],patient_56[which_var],patient_76[which_var]))
par(mar=c(5.1,13.1,2.5,2.1))
barplot(values,
        main = var.name,
        xlab = expression(bold(Count)),
        #ylab = "Patient numbers",
        names.arg = c("Mean", "Patient 4  ", "Patient 40","Patient 45","Patient 51","Patient 56","Patient 76"),
        col = colors_here,
        las=1,
        horiz = TRUE, cex.names = 3, font = 2, cex.lab = 2)


which_var <- 5
var.name <- c("MUT.TP53")
colors_here <- c("darkgray","darkolivegreen3", "deeppink1","blue1", "chartreuse3", "black", "red")
values <- as.numeric(c(SYNTHETIC_Mutation_mean[which_var],patient_4[which_var],patient_40[which_var],patient_45[which_var],
                       patient_51[which_var],patient_56[which_var],patient_76[which_var]))
par(mar=c(5.1,13.1,2.5,2.1))
barplot(values,
        main = var.name,
        xlab = expression(bold(Count)),
        #ylab = "Patient numbers",
        names.arg = c("Mean", "Patient 4  ", "Patient 40","Patient 45","Patient 51","Patient 56","Patient 76"),
        col = colors_here,
        las=1,
        horiz = TRUE, cex.names = 3, font = 2, cex.lab = 2)


which_var <- 7
var.name <- c("MUT.TMB")
colors_here <- c("darkgray","darkolivegreen3", "deeppink1","blue1", "chartreuse3", "black", "red")
values <- as.numeric(c(SYNTHETIC_Mutation_mean[which_var],patient_4[which_var],patient_40[which_var],patient_45[which_var],
                       patient_51[which_var],patient_56[which_var],patient_76[which_var]))
par(mar=c(5.1,13.1,2.5,2.1))
barplot(values,
        main = var.name,
        xlab = expression(bold(Count)),
        #ylab = "Patient numbers",
        names.arg = c("Mean", "Patient 4  ", "Patient 40","Patient 45","Patient 51","Patient 56","Patient 76"),
        col = colors_here,
        las=1,
        horiz = TRUE, cex.names = 3, font = 2, cex.lab = 2)



