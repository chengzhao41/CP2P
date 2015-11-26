library("sva")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("gridBase")
source('~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Script/preparing_data_helper.R')
source('~/Dropbox/SNF_DRUG_PROJECT/Script/comGENE.R')

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Main Figure Scripts")

##########
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/epirubicin_data.RData")
source('~/Dropbox/SNF_DRUG_PROJECT/Script/comGENE.R')
temp.data <- comGENE(scale(epirubicin$heiser), scale(epirubicin$patient))
dim(temp.data[[1]])
dim(temp.data[[2]])
mean(temp.data[[1]])
mean(temp.data[[2]])

temp.label <- epirubicin.labels$AUC_combined.source
temp.label[temp.label == "heiser"] <- "Heiser Cell Lines"
temp.label[temp.label == "patient"] <- "Patient"

p1 <- show_pca(input_data = rbind(temp.data[[1]], temp.data[[2]]), label = temp.label, pca_line_plot = FALSE, print = FALSE)
p1 <- p1 + theme(legend.direction = 'horizontal', 
                 legend.position = c(0.5, 1.2),#'top', 
                 plot.margin = unit(c(0,3,-10,2)/2, "lines"), 
                 text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2)
                 )  + scale_color_discrete(name = 'Dataset')
p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17))))
p1 <- p1 + geom_point(aes(shape = factor(temp.label), colour = factor(temp.label)), show_guide = FALSE)
p1

p2 <- show_pca(input_data = epirubicin$combined_AUC.sva, label = temp.label, pca_line_plot = FALSE, print = FALSE)
p2 <- p2 + theme(legend.direction = 'horizontal', 
                 legend.position = c(0.5, 1.2),#'top', 
                 plot.margin = unit(c(0,3,-10,2)/2, "lines"),
                 text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2)                 
                 ) + scale_color_discrete(name = 'Dataset')
p2 <- p2 + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17))))
p2 <- p2 + geom_point(aes(shape = factor(temp.label), colour = factor(temp.label)), show_guide = FALSE)
p2

#########
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Epirubicin/WS/epi_best_plot.RData")

library("ggplot2")
library("grid")
library("RColorBrewer")

cpp.m.slope <- data.frame(cbind(cpp.best$x_values, cpp.best$slope.mean, cpp.best$slope.sdev))
colnames(cpp.m.slope) <- c("xval", "yval", "se")
cpp.m.slope <- cpp.m.slope[2:10, ]
cpp.m.auc <- data.frame(cbind(cpp.best$x_values, cpp.best$auc.mean, cpp.best$auc.sdev))
colnames(cpp.m.auc) <- c("xval", "yval", "se")
cpp.m.auc <- cpp.m.auc[2:10, ]

pp.m <- cbind(pp.best$x_values, y_values = pp.best$mean, sdev = pp.best$sdev)
colnames(pp.m) <- c("xval", "yval", "se")

all.m <- rbind(cpp.m.slope, cpp.m.auc, pp.m)
all.m$Approach <- c(rep(cpp.best$slope.name, dim(cpp.m.slope)[1]), rep(cpp.best$auc.name, dim(cpp.m.auc)[1]), rep(pp.best$name, dim(pp.m)[1]))

all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
all.m$Approach[nrow(all.m)] <- cp.best$auc.name
all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
all.m$Approach[nrow(all.m)] <- cp.best$slope.name

all.m$Approach <- factor(all.m$Approach, levels = c(cp.best$auc.name, cp.best$slope.name,
                                                    pp.best$name, 
                                                    cpp.best$auc.name, cpp.best$slope.name))

pd <- position_dodge(2) # move them .05 to the left and right

p4 <- ggplot(all.m, aes(x=xval, y=yval, colour = Approach, ymax = 0.7)) + theme_bw() +   
  geom_errorbar(aes(ymin= yval - se, ymax = yval + se), width=5, position=pd) + 
  geom_line(position=pd) + 
#   geom_line() +
  geom_point(aes(shape=Approach, colour = Approach), size = 4, na.rm = TRUE) +
  geom_hline(aes(yintercept = cp.best$slope, colour = cp.best$slope.name), show_guide = FALSE) +   
  geom_hline(aes(yintercept = cp.best$auc, colour = cp.best$auc.name), show_guide = FALSE) +   
  scale_color_manual(values =  brewer.pal(n = 5, name = 'Dark2')) + 
  scale_shape_manual(values = c("C2P-AUC SNF all" = NA, 
                                "C2P-Slope SNF L1000" = NA,                                 
                                "P2P SVM lin mRMR1000" = 16,                                
                                "CP2P-AUC SNF mRMR1000" = 17, 
                                "CP2P-Slope SNF mRMR1000" = 18)) +
  scale_y_continuous(breaks = seq(0.3, 1, 0.05), "Test AUROC") +
  scale_x_continuous(breaks = seq(10, 100, by = 10), "# Number of Patient Samples in Training")  
p4

p4 <- p4 + theme(legend.direction = 'horizontal', 
                 legend.position = c(0.5, 1.2),#'top', 
                 #legend.position = 'top', 
                 plot.margin = unit(c(5.1, 7, 2, 3.5)/2, "lines"), 
                 text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2)) 
p4 <- p4 + guides(col=guide_legend(ncol=2))
p4

##############

pdf("epi_plot2.pdf", width = 21, height = 7, onefile = FALSE)
par(mfrow=c(1, 3))
plot.new()

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p4)

g <- gtable:::cbind_gtable(g1, g2, "first")
g <- gtable:::cbind_gtable(g, g3, "first")

grid.arrange(g)

#mtext(expression(bold('a')), adj=-1.8, side=3, outer=F, cex=1.5, line = 2)  
#mtext(expression(bold('b')), adj=-1.8, side=3, outer=F, cex=1.5, line = -34)  
mtext(expression(bold('A')), adj=-0.05, side=3, outer=F, cex=1.5, line = 1)
plot.new()
mtext(expression(bold('B')), adj=-0.1, side=3, outer=F, cex=1.5, line = 1)
plot.new()
mtext(expression(bold('C')), adj=-0.1, side=3, outer=F, cex=1.5, line = 1)

dev.off()