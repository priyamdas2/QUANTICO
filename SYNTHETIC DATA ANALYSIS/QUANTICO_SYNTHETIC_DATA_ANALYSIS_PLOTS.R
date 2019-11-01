rm(list=ls())
# install.packages("fda")
# install.packages("StatMethRank")
# install.packages("LaplacesDemon")
# install.packages("ggplot2")
# install.packages("TeachingDemos")
# install.packages("ggmap")
# install.packages("gplots")
# install.packages("pheatmap")


library(fda)
library(StatMethRank)
library(ggplot2)
library(ggmap)
library(gplots)
library(pheatmap)
library(LaplacesDemon)
library(TeachingDemos)
library(grid)
library(RSelenium)
library(ggplot2)
library(plyr)


source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/boxplot.with.outlier.label.r")
setwd("D:/QUANTICO/QUANTICO REPRODUCIBLE CODES/SYNTHETIC DATA ANALYSIS")

iters <- 50000

file <- paste("AAA_to_R_mean_probs_patient_T_prod_cells_iters_",iters,".csv", sep="")
mean_probs_patient_T_prod_cells <- read.csv(file, header=FALSE, sep=",")

file <- paste("AAA_to_R_mean_coef_patient_wise_prod_cells_iters_",iters,".csv", sep="")
mean_coef_patient_wise_prod_cells <- read.csv(file, header=FALSE, sep=",")

file <- paste("AAA_to_R_Mean_probs_G_TCR_prod_cells_iters_",iters,".csv", sep="")
Mean_probs_G_TCR_prod_cells <- read.csv(file, header=FALSE, sep=",")

file <- paste("AAA_to_R_Mean_probs_P_iters_",iters,".csv", sep="")
Mean_probs_P <- read.csv(file, header=FALSE, sep=",")

file <- paste("SYNTHETIC_clinical.csv", sep="")
TCR_clinical <- read.csv(file, header=TRUE, sep=",")



########## Heatmap of patients post. coef values #############
mat_here_temp <- t(as.matrix(mean_coef_patient_wise_prod_cells))

mat_here <- mat_here_temp
#mat_here <- logit(mat_here)



#row.names(TCR_clinical) <- c("MDA.Case.ID", "Diagnosis", "Gender", "Death", "Recurrence", "Tobacco", "Smoker")

annotation <- data.frame(#Tobacco = factor(TCR_clinical[,6]),
                         Death = factor(TCR_clinical[,4]),
                         Recurrence = factor(TCR_clinical[,5]),
                         Smoker = factor(TCR_clinical[,7]),
                         Stage = factor(TCR_clinical[,8]))
rownames(annotation) <- 1:87

## Change 1 patient TYPE.SMOKER to never second hand
ann_colors <- list(
  #Tobacco = c(Y = "orange", N =  "paleturquoise3"),
  Death = c("TRUE" = "black", "FALSE" =  "darkorchid1"),
  Recurrence = c(Y = "darkorange", N =  "azure4"),
  #Smoker = c("Current" = "red", "Former" = "yellow", "Never" = "white"),
  Smoker = c("Current" = "white", "Former" = "cyan3", "Never" = "darkblue"),
  Stage = c("I" = "green", "II" = "yellow", "III" = "red"))

colnames(mat_here) <- c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9" )
rownames(mat_here) <-  1:87


#setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=1, height=1, name="vp", just=c("right","top"))), action="prepend")

pheatmap(mat_here[,-9], annotation_row = annotation, annotation_colors = ann_colors,angle_col = "90",fontsize_col = 16,
         treeheight_col = 0, show_rownames = FALSE, show_colnames = TRUE, cluster_cols = FALSE, 
         #main = "Heatmap - Coeff. of T.prod.cells on CD8",
         xlab = "Patients", ylab = "quantiles")
setHook("grid.newpage", NULL, "replace")
grid.text("Quantiles", y=0.03, gp=gpar(fontsize=20))


# heatmap.2(as.matrix(mean_probs_patient_T_prod_cells), Rowv = FALSE, colv = TRUE,
#           xlab = "Patients", ylab = "quantiles")

#######################################################





#### Boxplot of coefficients of T.Prod.Cells
quantiles <- seq(0.1,.9,by = .1)

n <- max(dim(mean_coef_patient_wise_prod_cells))
array_coefs <- as.vector(t(as.matrix(mean_coef_patient_wise_prod_cells)))
array_quantiles <- as.vector(t(as.matrix(quantiles*matrix(1,9,n))))

medians <- rep(0,9)

for(i in 1:9)
{medians[i] <-  median(t(mean_coef_patient_wise_prod_cells[i,]))}


boxplot.with.outlier.label(array_coefs ~ array_quantiles,xlab="Quantiles", label_name = 1:87,ylim = c(0,9),
        main = "Effect of T.prod.cells on CD8 at different quantile levels",
        ylab="Coefficient",
        col="orange",
        border="brown", cex.lab=1.5)
lines(c(1:9),medians)


##############################################


############ G Selection Plot ####################################

Mean_probs_G_TCR_prod_cells_ARRANGED <- rbind(rep(0,9),Mean_probs_G_TCR_prod_cells)
Mean_probs_G_TCR_prod_cells_ARRANGED[1,] <- Mean_probs_G_TCR_prod_cells[7,]
Mean_probs_G_TCR_prod_cells_ARRANGED <- Mean_probs_G_TCR_prod_cells_ARRANGED[-8,]


set.seed(45)
df <- data.frame(x=rep(quantiles, 7), val=as.vector(t(as.matrix(Mean_probs_G_TCR_prod_cells_ARRANGED))), 
                 variable=rep(c("TMB","TTN","MUC16", "RYR2", "CSMD3", "TP53", "USH2A"), each=9))

# variable=c("TMB","TTN","MUC16", "RYR2", "CSMD3", "TP53", "USH2A")

plot_here <- ggplot(data = df, aes(x=x, y=val), 
                    ylab = "Post. prob. of selction") + 
  geom_line(data = df, aes(colour=variable), linetype = "longdash", size = 1.1) +
  geom_point(aes(colour=variable),size = 3)



plot_here_2 <- plot_here +xlab("Quantiles")+ylab("Post. prob. of selction")


plot_here_2+geom_hline(yintercept=0.5, linetype="dotted", color = "red")+
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1))+
  # theme(axis.text=element_text(size=50,face="bold"),
  #       axis.title=element_text(size=50,face="bold"))+
  theme(legend.text=element_text(size=50))+
  theme(legend.title = element_blank())+ theme_bw(base_size=25)

##################################################################

###### P Selection Plot ############################################


Mean_probs_P_ARRANGED <- t(Mean_probs_P[,-1])


set.seed(45)
df <- data.frame(x=rep(quantiles, 4), val=as.vector(t(as.matrix(Mean_probs_P_ARRANGED))), 
                 variable=rep(c("T.clonality","T.entropy","T.prod.cells", "T.richness"), each=9))

# variable=c("TMB","TTN","MUC16", "RYR2", "CSMD3", "TP53", "USH2A")

plot_here <- ggplot(data = df, aes(x=x, y=val), 
                    ylab = "Post. prob. of selction") + 
                    geom_line(data = df, aes(colour=variable), linetype = "longdash", size = 1.1) +
                    geom_point(aes(colour=variable),size = 3)



plot_here_2 <- plot_here +xlab("Quantiles")+ylab("Post. prob. of selction")


plot_here_2+geom_hline(yintercept=0.5, linetype="dotted", color = "red")+
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1))+
# theme(axis.text=element_text(size=50,face="bold"),
#       axis.title=element_text(size=50,face="bold"))+
  theme(legend.text=element_text(size=50))+
  theme(legend.title = element_blank())+ theme_bw(base_size=25)


###################  Posterior Probability based on recurrnece, death etc categorization ################
# Recurrence

head(TCR_clinical)

cat_1_recur <- which(TCR_clinical[,5] == "Y")
cat_2_recur <- which(TCR_clinical[,5] == "N")

rowMeans(mean_probs_patient_T_prod_cells[,cat_1_recur])
rowMeans(mean_probs_patient_T_prod_cells[,cat_2_recur])

# Death 

cat_1_vital <- which(TCR_clinical[,4] == "TRUE")
cat_2_vital <- which(TCR_clinical[,4] == "FALSE")

rowMeans(mean_probs_patient_T_prod_cells[,cat_1_vital])
rowMeans(mean_probs_patient_T_prod_cells[,cat_2_vital])

# Tobacco 

cat_1_tobac <- which(TCR_clinical[,6] == "Y")
cat_2_tobac <- which(TCR_clinical[,6] == "N")

rowMeans(mean_probs_patient_T_prod_cells[,cat_1_tobac])
rowMeans(mean_probs_patient_T_prod_cells[,cat_2_tobac])

# Smoker

cat_1_smoke <- which(TCR_clinical[,7] == "Current")
cat_2_smoke <- which(TCR_clinical[,7] == "Former")
cat_3_smoke <- which(TCR_clinical[,7] == "Never")

rowMeans(mean_probs_patient_T_prod_cells[,cat_1_smoke])
rowMeans(mean_probs_patient_T_prod_cells[,cat_2_smoke])
rowMeans(mean_probs_patient_T_prod_cells[,cat_3_smoke])


MAT_ALL <- rbind(rowMeans(mean_probs_patient_T_prod_cells[,cat_1_recur]),
      rowMeans(mean_probs_patient_T_prod_cells[,cat_2_recur]),
rowMeans(mean_probs_patient_T_prod_cells[,cat_1_vital]),
rowMeans(mean_probs_patient_T_prod_cells[,cat_2_vital]),
rowMeans(mean_probs_patient_T_prod_cells[,cat_1_tobac]),
rowMeans(mean_probs_patient_T_prod_cells[,cat_2_tobac]),
rowMeans(mean_probs_patient_T_prod_cells[,cat_1_smoke]),
rowMeans(mean_probs_patient_T_prod_cells[,cat_2_smoke]),
rowMeans(mean_probs_patient_T_prod_cells[,cat_3_smoke]))

rowMeans(MAT_ALL)

MAT_ALL_final <- round(cbind(MAT_ALL,rowMeans(MAT_ALL)),3)

# write.table( MAT_ALL_final,"Category_wise_post_prob.csv", sep=",", row.names = FALSE, col.names=FALSE)


###################  Reg coeff based on recurrnece, death etc categorization ################
# Recurrence

head(TCR_clinical)

cat_1_recur <- which(TCR_clinical[,5] == "Y")
cat_2_recur <- which(TCR_clinical[,5] == "N")

rowMeans(mean_coef_patient_wise_prod_cells[,cat_1_recur])
rowMeans(mean_coef_patient_wise_prod_cells[,cat_2_recur])

# Death 

cat_1_vital <- which(TCR_clinical[,4] == "TRUE")
cat_2_vital <- which(TCR_clinical[,4] == "FALSE")

rowMeans(mean_coef_patient_wise_prod_cells[,cat_1_vital])
rowMeans(mean_coef_patient_wise_prod_cells[,cat_2_vital])

# Stage 

cat_1_stage <- which(TCR_clinical[,8] == "I")
cat_2_stage <- which(TCR_clinical[,8] == "II")
cat_3_stage <- which(TCR_clinical[,8] == "III")


rowMeans(mean_coef_patient_wise_prod_cells[,cat_1_stage])
rowMeans(mean_coef_patient_wise_prod_cells[,cat_2_stage])
rowMeans(mean_coef_patient_wise_prod_cells[,cat_3_stage])

# Smoker

cat_1_smoke <- which(TCR_clinical[,7] == "Current")
cat_2_smoke <- which(TCR_clinical[,7] == "Former")
cat_3_smoke <- which(TCR_clinical[,7] == "Never")

rowMeans(mean_coef_patient_wise_prod_cells[,cat_1_smoke])
rowMeans(mean_coef_patient_wise_prod_cells[,cat_2_smoke])
rowMeans(mean_coef_patient_wise_prod_cells[,cat_3_smoke])


MAT_ALL <- rbind(rowMeans(mean_coef_patient_wise_prod_cells[,cat_1_recur]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_2_recur]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_1_vital]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_2_vital]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_1_stage]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_2_stage]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_3_stage]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_1_smoke]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_2_smoke]),
                 rowMeans(mean_coef_patient_wise_prod_cells[,cat_3_smoke]))

rowMeans(MAT_ALL)

MAT_ALL_final <- round(cbind(MAT_ALL,rowMeans(MAT_ALL)),3)

# write.table( MAT_ALL_final,"Category_wise_reg_coeff.csv", sep=",", row.names = FALSE, col.names=FALSE)

########## DENSITY PLOTS ###################


# Smoking DENSITY PLOTS

for(quantile in 1:8)
{cat_1 <- as.numeric(mean_coef_patient_wise_prod_cells[quantile,cat_1_smoke])
cat_2 <- as.numeric(mean_coef_patient_wise_prod_cells[quantile,cat_2_smoke])
cat_3 <- as.numeric(mean_coef_patient_wise_prod_cells[quantile,cat_3_smoke])


df1 <- data.frame(x = cat_1, label=rep('Current', length(cat_1)))
df2 <- data.frame(x = cat_2, label=rep('Former', length(cat_2)))
df3 <- data.frame(x = cat_3, label=rep('Never', length(cat_3)))
df=rbind(df1, df2)
df=rbind(df, df3)
temp_plot <- ggplot(df, aes(x, y=..scaled.., fill=label)) + 
  geom_density(color = "black", alpha = 0.7) + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9","green"))+
  scale_x_continuous(name = "Effect size")+
  scale_y_continuous(name = "Density")+
  theme_bw()+
  theme_grey(base_size = 22)
print(temp_plot)
ggsave(temp_plot, file=paste0("figures/plot_SMOKE_", quantile,".png"), width = 15, height = 10, units = "cm")}



