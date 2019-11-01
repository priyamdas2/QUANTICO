rm(list=ls())
#install.packages("fda")
#install.packages("StatMethRank")
# install.packages("LaplacesDemon")
#install.packages("ggplot2")
# install.packages("ggmap")
# install.packages("TeachingDemos")
# install.packages("RPMG")
# install.packages("rmarkdown")
# install.packages("knitr")
# install.packages("gdata")
# install.packages("readxl")
# install.packages("WriteXLS")
# install.packages("tikzDevice")
# install.packages("latex2exp")
#library(rmarkdown)
#library(knitr)
#library(gdata)
#library(readxl)
#library(WriteXLS)
#library(tikzDevice)
#library(latex2exp)
#library(ggmap)
#library(gplots)
#library(pheatmap)
#library(StatMethRank)

require(graphics)
library(fda)
library(ggplot2)
library(LaplacesDemon)
library(TeachingDemos)
library(grid)
library(RPMG)

# options(tikzMetricPackages = c("\\usepackage{amsmath}",
#                                "\\usepackage{xcolor}",
#                                "\\usepackage{tikz}",
#                                "\\usetikzlibrary{calc}"))
# tikz("test-tikz.tex",standAlone = TRUE,packages=c("\\usepackage{amsmath}",
#                                                   "\\usepackage{tikz}",
#                                                   "\\usepackage{xcolor}",
#                                                   "\\usetikzlibrary{calc}",
#                                                   "\\usepackage[active,tightpage,psfixbb]{preview}",
#                                                   "\\PreviewEnvironment{pgfpicture}"))
# 


setwd("Y:/QUANTICO/QUANTICO REPRODUCIBLE CODES/SIMULATION_PLOTS")

#################################
### READ TCR Merged Data
#################################

sample_size <- 200


All_P_probs  <- read.table('for_R_plot_P.csv', header=FALSE, sep = ",")

G1_probs  <- read.table('for_R_plot_G1.csv', header=FALSE, sep = ",")

G2_probs  <- read.table('for_R_plot_G2.csv', header=FALSE, sep = ",")

QR_Y_sim_ESTs  <- read.table('QR_Y_sim_ESTs_200.csv', header=FALSE, sep = ",")

QR_Y_sim_TRUEs  <- read.table('QR_Y_sim1_TRUEs_200.csv', header=FALSE, sep = ",")

QR_Y_sim_ESTs_high  <- read.table('QR_Y_sim_ESTs_high_200.csv', header=FALSE, sep = ",")

QR_Y_sim_ESTs_low <- read.table('QR_Y_sim_ESTs_low_200.csv', header=FALSE, sep = ",")



#### PLOT 1

# colors_here <- c("brown1", "blue", "darkgreen", "darkorange", "darkviolet")
# colors_here <- rev(rainbow.colors(5))
# par(mar=c(5.1,5.1,4.1,2.1))
# matplot(All_P_probs[,2:6], type = c("b"),pch=3,col = colors_here, lty = 1,font.lab=2,
#         xlab = "Quantiles", ylab = "Post. selection prob.",axes=F, cex.lab=1.5)
# quantiles <- seq(from = 0.1, to = 0.9, by = 0.1)
# axis(side=1,at=1:9,labels=quantiles)
# axis(2)
# legend(x = 6.5,y = 0.8, 
#        legend = c(expression(`$P_1$ (TRUE)` = paste("", "P", phantom()[{paste("1")}],
#                                                     "", " ", "(", "TRUE", ")", "")),
#                   expression(`$P_2$ (TRUE)` = paste("", "P", phantom()[{paste("2")}],
#                                                     "", " ", "(", "TRUE", ")", "")),
#                   expression(`$P_3$ (NOISE/TRUE)` = paste("", "P", phantom()[{paste("3")}],
#                                                           "", " ", "(", "NOISE/TRUE", ")", "")),
#                   expression(`$P_4$ (NOISE)` = paste("", "P", phantom()[{paste("4")}],
#                                                      "", " ", "(", "NOISE", ")", "")),
#                   expression(`$P_5$ (NOISE)` = paste("", "P", phantom()[{paste("5")}],
#                                                      "", " ", "(", "NOISE", ")", ""))),
#        col = colors_here, pch=3, lty = 1, cex=1.2)


#### GGPLOT 
par(mar=c(5.1,5.1,4.1,2.1))
set.seed(45)
mat_now <- t(as.matrix(All_P_probs[,2:6]))

quantiles <- seq(0.1,.9,by = .1)
df <- data.frame(x=rep(quantiles, 5), val=as.vector(t(as.matrix(mat_now))), 
                 variable=rep(c("P_1 (TRUE)","P_2 (TRUE)","P_3 (NOISE/TRUE)","P_4 (NOISE)","P_5 (NOISE)"),each=9)
                   # c(paste("", "P", phantom()[{paste("1")}],
                   #                                                "", " ", "(", "TRUE", ")", "")),
                   #              paste("", "P", phantom()[{paste("2")}],
                   #                                                "", " ", "(", "TRUE", ")", "")),
                   #              paste("", "P", phantom()[{paste("3")}],
                   #                                                      "", " ", "(", "NOISE/TRUE", ")", ""),
                   #              paste("", "P", phantom()[{paste("4")}],
                   #                                                 "", " ", "(", "NOISE", ")", ""),
                   #              paste("", "P", phantom()[{paste("5")}],
                   #                                                 "", " ", "(", "NOISE", ")", ""),
)

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
  # theme(legend.text=element_text(size=50))+
  # theme(legend.title = element_blank())+ 
  theme_bw(base_size=25)+theme(legend.position = "none") 


#######################################################
#### PLOT 2 (G selection for P_1)


# colors_here <- rev(rainbow.colors(6))
# par(mar=c(5.1,5.1,4.1,2.1))
# matplot(t(G1_probs), type = c("b"),pch=3,col = colors_here, lty = 1,
#         xlab = "Quantiles", ylab = "Post. selection prob.",axes=F, font.lab=2, cex.lab=1.5)
# quantiles <- seq(from = 0.1, to = 0.9, by = 0.1)
# axis(side=1,at=1:9,labels=quantiles)
# axis(2)
# legend(x = 7,y = 0.9, legend = c(expression(`$G_1$ (TRUE)` = paste("", "G", phantom()[{paste("1")}],
#                                                                    "", " ", "(", "TRUE", ")", "")),
#                                  expression(`$G_2$ (NOISE)` = paste("", "G", phantom()[{paste("2")}],
#                                                                    "", " ", "(", "NOISE", ")", "")),
#                                  expression(`$G_3$ (NOISE)` = paste("", "G", phantom()[{paste("3")}],
#                                                                          "", " ", "(", "NOISE", ")", "")),
#                                  expression(`$G_4$ (NOISE)` = paste("", "G", phantom()[{paste("4")}],
#                                                                     "", " ", "(", "NOISE", ")", "")),
#                                  expression(`$G_5$ (NOISE)` = paste("", "G", phantom()[{paste("5")}],
#                                                                     "", " ", "(", "NOISE", ")", "")),
#                                  expression(`$G_6$ (NOISE)` = paste("", "G", phantom()[{paste("6")}],
#                                                                     "", " ", "(", "NOISE", ")", ""))),
#        col = colors_here, pch=3, lty = 1, cex=1.2)


#### GGPLOT 
par(mar=c(5.1,5.1,4.1,2.1))
set.seed(45)
mat_now <- G1_probs

quantiles <- seq(0.1,.9,by = .1)
df <- data.frame(x=rep(quantiles, 6), val=as.vector(t(as.matrix(mat_now))), 
                 variable=rep(c("G_1 (TRUE)","G_2 (NOISE)","G_3 (NOISE)","G_4 (NOISE)","G_5 (NOISE)"
                                ,"G_6 (NOISE)"),each=9))
                 


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
  # theme(legend.text=element_text(size=50))+
  # theme(legend.title = element_blank())+ 
  theme_bw(base_size=25)+theme(legend.position = "none") 



#######################################################
#### PLOT 3 (G selection for P_4)


# colors_here <- rev(rainbow.colors(6))
# par(mar=c(5.1,5.1,4.1,2.1))
# matplot(t(G2_probs), type = c("b"),pch=3,col = colors_here, lty = 1,
#         xlab = "Quantiles", ylab = "Post. selection prob.",axes=F, font.lab=2, cex.lab=1.5)
# quantiles <- seq(from = 0.1, to = 0.9, by = 0.1)
# axis(side=1,at=1:9,labels=quantiles)
# axis(2)
# legend(x = 7,y = 0.9, legend = c(expression(`$G_1$ (NOISE)` = paste("", "G", phantom()[{paste("1")}],
#                                                                    "", " ", "(", "NOISE", ")", "")),
#                                  expression(`$G_2$ (TRUE)` = paste("", "G", phantom()[{paste("2")}],
#                                                                     "", " ", "(", "TRUE", ")", "")),
#                                  expression(`$G_3$ (NOISE)` = paste("", "G", phantom()[{paste("3")}],
#                                                                     "", " ", "(", "NOISE", ")", "")),
#                                  expression(`$G_4$ (NOISE)` = paste("", "G", phantom()[{paste("4")}],
#                                                                     "", " ", "(", "NOISE", ")", "")),
#                                  expression(`$G_5$ (NOISE)` = paste("", "G", phantom()[{paste("5")}],
#                                                                     "", " ", "(", "NOISE", ")", "")),
#                                  expression(`$G_6$ (NOISE)` = paste("", "G", phantom()[{paste("6")}],
#                                                                     "", " ", "(", "NOISE", ")", ""))),
#        col = colors_here, pch=3, lty = 1, cex=1.2)
# 

#### GGPLOT 
par(mar=c(5.1,5.1,4.1,2.1))
set.seed(45)
mat_now <- G2_probs

quantiles <- seq(0.1,.9,by = .1)
df <- data.frame(x=rep(quantiles, 6), val=as.vector(t(as.matrix(mat_now))), 
                 variable=rep(c("G_1 (TRUE)","G_2 (NOISE)","G_3 (NOISE)","G_4 (NOISE)","G_5 (NOISE)"
                                ,"G_6 (NOISE)"),each=9))



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
  # theme(legend.text=element_text(size=50))+
  # theme(legend.title = element_blank())+ 
  theme_bw(base_size=25)+theme(legend.position = "none") 



#######################################################
#### PLOT 4 Patient_wise confidence bands


# for(i in 1:20)
#   {which_patient <- i
# 
# plot(0,0,xlim = c(0.1,0.9),ylim = c(-5,25),type="l", axes=F,
#      xlab = "quantiles", ylab = "Y", main = paste0("patient_no ",which_patient))
# axis(side=1,at=quantiles,labels=quantiles)
# axis(2)
# lines(quantiles,QR_Y_sim_TRUEs[which_patient,], col = "red", lty = "solid" )
# lines(quantiles,QR_Y_sim_ESTs[which_patient,], col = "blue", lty = "solid")
# lines(quantiles,QR_Y_sim_ESTs_high[which_patient,], col = "blue", lty = "dashed" )
# lines(quantiles,QR_Y_sim_ESTs_low[which_patient,], col = "blue", lty = "dashed")}




random_patients <- c(85, 88, 197)

#TeX("Q(\\tau|x)")


for(i in random_patients)
{which_patient <- i
par(mar=c(5.1,5.1,1.1,1.1))
#par(mar=c(5.1,5.1,4.1,2.1))
plot(0,0,xlim = c(0.1,0.9),ylim = c(-3,33),type="l", axes=F,
     xlab = "Quantiles", ylab = expression(`$Q(\tau|X)$` = bold(paste("Q", "(", bold(tau), "|", X, ")", "", ""))),
     font.lab=2,cex.lab=1.5)
axis(side=1,at=quantiles,labels=quantiles)
axis(2)
lines(quantiles,QR_Y_sim_TRUEs[which_patient,], col = "red", lty = "solid", lwd=2.5)
lines(quantiles,QR_Y_sim_ESTs[which_patient,], col = "blue",lty = 1,cex=3, lwd=2.5)
points(quantiles,QR_Y_sim_ESTs[which_patient,], col = "blue",pch=19,lty = 1, cex=1.5)
lines(quantiles,QR_Y_sim_ESTs_high[which_patient,], col = "blue", lty = "dashed", pch=3, lwd=2.5)
lines(quantiles,QR_Y_sim_ESTs_low[which_patient,], col = "blue", lty = "dashed",pch=3, lwd=2.5)

legend("topleft",cex=1.2,legend = c("True quantiles", "Est. quantiles", "95% credible bounds"),col = c("red","blue","blue"), lty = c(1,1,2), lwd = c(2.5,2.5,2.5))

}

# 
# which_patient <- 86
# par(mar=c(5.1,5.1,4.1,2.1))
# plot(0,0,xlim = c(0.1,0.9),ylim = c(-5,35),type="l", axes=F,
#      xlab = "Quantiles", ylab = expression(`$Q(\tau|X)$` = bold(paste("Q", "(", bold(tau), "|", X, ")", "", ""))),
#      font.lab=2,cex.lab=1.5)
# axis(side=1,at=quantiles,labels=quantiles)
# axis(2)
# lines(quantiles,QR_Y_sim_TRUEs[which_patient,], col = "red", lty = "solid" )
# lines(quantiles,QR_Y_sim_ESTs[which_patient,], col = "blue",pch=3,lty = 1)
# lines(quantiles,QR_Y_sim_ESTs_high[which_patient,], col = "blue", lty = "dashed", pch=3)
# lines(quantiles,QR_Y_sim_ESTs_low[which_patient,], col = "blue", lty = "dashed",pch=3)
# 
# legend("topleft",cex=1.2,legend = c("True quantiles", "Est. quantiles", "95% credible bounds"),col = c("red","blue","blue"), lty = c(1,1,2))
# 


