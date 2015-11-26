library("sva")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("gridBase")
source('~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Script/preparing_data_helper.R')
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/docetaxel_data.RData")

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Main Figure Scripts")

# AUC
#show_pca(input_data = docetaxel$breast_AUC.ComBat, label = docetaxel.labels$AUC_breast)
temp.labels <- docetaxel.labels$AUC_breast.source
temp.labels[temp.labels == "cgp"] = "CGP breast cell lines"
temp.labels[temp.labels == "patient"] = "Patient"
p1 <- show_pca(input_data = docetaxel$breast_AUC.ComBat, label = temp.labels, pca_line_plot = FALSE, print = FALSE)
p1 <- p1 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top', 
                 plot.margin = unit(c(5.1,7,4.5,3.5)/2, "lines"), text = element_text(size=15), axis.title.x=element_text(vjust=-1.5)) + 
  scale_color_discrete(name = 'Dataset') 
p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17))))
p1 <- p1 + geom_point(aes(shape = factor(temp.labels), colour = factor(temp.labels)), show_guide = FALSE)
p1

g1 <- ggplotGrob(p1)

temp.labels <- docetaxel.labels$AUC_combined.source
temp.labels[temp.labels == "cgp"] = "CGP"
temp.labels[temp.labels == "patient"] = "Patient"
p2 <- show_pca(input_data = docetaxel$combined_AUC.ComBat, label = temp.labels, pca_line_plot = FALSE, print = FALSE)
p2 <- p2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top', 
                 plot.margin = unit(c(5.1,7,4.5,3.5)/2, "lines"), text = element_text(size=15), axis.title.x=element_text(vjust=-1.5)) + 
  scale_color_discrete(name = 'Dataset') 
p2 <- p2 + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17))))
p2 <- p2 + geom_point(aes(shape = factor(temp.labels), colour = factor(temp.labels)), show_guide = FALSE)
p2

g2 <- ggplotGrob(p2)


####### bar graph
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/doc_best_plot.RData")
gg.cpp.auc.mean <- mean(cpp.best$AUC.test_auc[14:23])
gg.cpp.auc.se <- sd(cpp.best$AUC.test_auc[14:23])

gg.pp.auc.mean <- mean(pp.best$test_auc)
gg.pp.auc.se <- sd(pp.best$test_auc)

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cp_boxplot_data.RData")
gg.cp.auc.mean <- mean(cp.boxplot_data[, 20])
gg.cp.auc.se <- sd(cp.boxplot_data[, 20])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cpp_var_auc_breast.RData")
temp.cpp_breast.best <- which(max(apply(cpp.varying_training_matrix, 2, mean)) == apply(cpp.varying_training_matrix, 2, mean))
gg.cpp_breast.auc.mean <- mean(cpp.varying_training_matrix[, temp.cpp_breast.best])
gg.cpp_breast.auc.se <- sd(cpp.varying_training_matrix[, temp.cpp_breast.best])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/cp_breast_only_boxplot.RData")
temp.cp_breast.best <- which(max(apply(cp.breast_only.boxplot_data , 2, mean)) == apply(cp.breast_only.boxplot_data , 2, mean))
gg.cp_breast.auc.mean <- mean(cp.breast_only.boxplot_data[, temp.cp_breast.best])
gg.cp_breast.auc.se <- sd(cp.breast_only.boxplot_data[, temp.cp_breast.best])

gg.all <- data.frame(mean = c(gg.cpp.auc.mean, gg.pp.auc.mean, gg.cp.auc.mean, gg.cpp_breast.auc.mean, gg.cp_breast.auc.mean),
                     se = c(gg.cpp.auc.se, gg.pp.auc.se, gg.cp.auc.se, gg.cpp_breast.auc.se, gg.cp_breast.auc.se), 
                     Approach = c("CP2P", "P2P", "C2P", "CP2P Breast", "C2P Breast"))
gg.all$Approach <- factor(gg.all$Approach, levels = c("C2P Breast", "C2P", "P2P", "CP2P Breast", "CP2P"))
#cp.breast_only.boxplot_data
#levels(gg.all$Approach) <- c("C2P Breast", "C2P", "P2P", "CP2P Breast", "CP2P")

cols <- brewer.pal(n = 5, name = 'Set1')
p3 <- ggplot(gg.all, aes(x = Approach, y = mean, fill = Approach)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + theme_bw() +  
  scale_y_continuous(breaks = seq(0, 0.95, 0.05), "Test AUROC", limit = c(0, 0.95))
p3 <- p3 + theme(legend.direction = 'horizontal', 
                 #legend.position = 'top', 
                 #legend.position = c(0.5, 1.075),#'top', 
                 legend.position="none",
                 #plot.margin = unit(c(5.1, 0, 4.5, 1)/2, "lines"), 
                 plot.margin = unit(c(10, 7, 2, 3.5)/2, "lines"), 
                 #plot.margin = unit(c(15, 7, 2, 3.5)/2, "lines"), 
                   text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2)) + 
  scale_fill_manual(name="Approach", 
                    labels = c("C2P Breast", "C2P", "P2P", "CP2P Breast", "CP2P"),                    
                    values = cols) + 
  scale_x_discrete(name ="Approach", labels=c("C2P Breast", "C2P", "P2P", "CP2P Breast", "CP2P"))
p3

########################### compare IC50, AUC, and slope ################################
library("ggplot2")
library("grid")
library("RColorBrewer")

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Docetaxel/WS/doc_best_plot.RData")
cl <- brewer.pal(n = 7, name = 'Dark2')
#cl <- cl[c(2, 5, 7, 1, 3, 6, 4)]

cpp.ic50.m <- data.frame(xval = cpp.best$x_values, yval = cpp.best$IC50.test_auc)
cpp.ic50.m <- cpp.ic50.m[14:23, ]
cpp.auc.m <- data.frame(xval = cpp.best$x_values, yval = cpp.best$AUC.test_auc)
cpp.auc.m <- cpp.auc.m[14:23, ]
cpp.slope.m <- data.frame(xval = cpp.best$x_values, yval = cpp.best$slope.test_auc)
cpp.slope.m <- cpp.slope.m[14:23, ]
pp.m <- data.frame(xval = pp.best$x_values, yval = pp.best$test_auc)

all.m <- rbind(cpp.ic50.m, cpp.auc.m, cpp.slope.m, pp.m)
all.m$Approach <- c(rep(cpp.best$IC50.name, dim(cpp.ic50.m)[1]),
                    rep(cpp.best$AUC.name, dim(cpp.auc.m)[1]),
                    rep(cpp.best$slope.name, dim(cpp.slope.m)[1]),
                    rep(pp.best$name, dim(pp.m)[1]))

all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
all.m$Approach[nrow(all.m)] <- cp.best$IC50.name
all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
all.m$Approach[nrow(all.m)] <- cp.best$AUC.name
all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
all.m$Approach[nrow(all.m)] <- cp.best$slope.name

all.m$Approach <- factor(all.m$Approach, levels = c(cp.best$IC50.name, cp.best$AUC.name, cp.best$slope.name,
                                                   pp.best$name, 
                                                   cpp.best$IC50.name, cpp.best$AUC.name, cpp.best$slope.name))

p4 <- ggplot(all.m, aes(x=xval, y=yval, colour = Approach, ymax = 0.95)) + theme_bw() + 
  geom_line() + 
  geom_point(aes(shape=Approach, colour = Approach), size = 4, na.rm = TRUE) +   
  geom_hline(aes(yintercept = cp.best$IC50, colour = cp.best$IC50.name)) + 
  geom_hline(aes(yintercept = cp.best$AUC, colour = cp.best$AUC.name)) +   
  geom_hline(aes(yintercept = cp.best$slope, colour = cp.best$slope.name)) +   
  scale_color_manual(values = cl) +
  scale_shape_manual(values = c("C2P-IC50 Elastic Net mRMR1000"=NA, 
                                "C2P-AUC SVM rbf mRMR1000"=NA, 
                                "C2P-Slope RF mRMR1000"=NA, 
                                "P2P SVM rbf mRMR1000"=16, 
                                "CP2P-IC50 SVM rbf mRMR1000"=17, 
                                "CP2P-AUC SVM rbf mRMR1000"=18, 
                                "CP2P-Slope SVM rbf mRMR1000"=19)) +
  scale_y_continuous(breaks = seq(0.65, 0.95, 0.05), "Test AUROC") +
  scale_x_continuous(breaks = seq(14, 23, by = 2), "# Number of Patient Samples in Training")
#p4
p4 <- p4 + theme(legend.direction = 'horizontal',                  
                 legend.position = c(0.5, 1.165),#'top', 
                 #plot.margin = unit(c(5.1, 0, 4.5, 1)/2, "lines"),                  
                 plot.margin = unit(c(15, 7, 2, 3.5)/2, "lines"),                  
                 text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2))
p4 <- p4 + guides(col=guide_legend(ncol=2))

############
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Main Figure Scripts")

require(gtable)
require(gridExtra)
# #Extract Grobs
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)

#Extract Grobs
# #Bind the tables
gT <- gtable:::cbind_gtable(g1, g2, "first")
gT <- gtable_add_rows(gT, unit(0.5,"cm"), pos=nrow(gT))
gB <- gtable:::cbind_gtable(g4, g3, "first")
#g <- gtable_add_cols(g, unit(5,"cm"), pos=nrow(g1))

pdf("doc_plot2.pdf", width = 14, height = 14, onefile=FALSE)

par(mfrow=c(2, 2))
plot.new()
plot.new()
plot.new()

gAll <- gtable:::rbind_gtable(gT, gB, "first")
#gAll <- gtable_add_rows(g, unit(5,"cm"), pos=nrow(g1))
grid.arrange(gAll)

#grid.arrange(gT, gB, nrow = 2, ncol = 1)

mtext(expression(bold('A')), adj=-0.1, side=3, outer=F, cex=1.5, line = 37)  
mtext(expression(bold('C')), adj=-0.1, side=3, outer=F, cex=1.5, line = 11)  

plot.new()
mtext(expression(bold('B')), adj=-0.1, side=3, outer=F, cex=1.5, line = 37)  
mtext(expression(bold('D')), adj=-0.1, side=3, outer=F, cex=1.5, line = 11)  

dev.off()
