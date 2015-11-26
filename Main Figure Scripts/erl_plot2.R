library("sva")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("gridBase")
source('~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Script/preparing_data_helper.R')

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erlotinib_data.RData")
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Main Figure Scripts")

#####
source('~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Script/preparing_data_helper.R')

temp.label <- erlotinib.labels$lung_all.IC50.source
temp.label[temp.label == "GEO"] = "BATTLE cell lines"
temp.label[temp.label == "patient"] = "Patient"
temp.label[temp.label == "CGP"] = "CGP lung cell lines"

p1 <- show_pca(input_data = erlotinib$lung_all.IC50.ComBat, label = temp.label, pca_line_plot = FALSE, print = FALSE)
p1 <- p1 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top', 
                 plot.margin = unit(c(5.1,7,4.5,3.5)/2, "lines"), text = element_text(size=15), axis.title.x=element_text(vjust=-1.5)) + 
  scale_color_discrete(name = 'Dataset') 
p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17, 18))))
p1 <- p1 + geom_point(aes(shape = factor(temp.label), colour = factor(temp.label)), show_guide = FALSE)
p1

g1 <- ggplotGrob(p1)


temp.label <- erlotinib.labels$all.source.slope
temp.label[temp.label == "GEO"] = "BATTLE cell lines"
temp.label[temp.label == "patient"] = "Patient"

p2 <- show_pca(input_data = erlotinib$all.ComBat.slope, label = temp.label, pca_line_plot = FALSE, print = FALSE)
p2 <- p2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top', 
                 plot.margin = unit(c(5.1,7,4.5,3.5)/2, "lines"), text = element_text(size=15), axis.title.x=element_text(vjust=-1.5)) + 
  scale_color_discrete(name = 'Dataset') 
p2 <- p2 + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17, 18))))
p2 <- p2 + geom_point(aes(shape = factor(temp.label), colour = factor(temp.label)), show_guide = FALSE)
p2

g2 <- ggplotGrob(p2)


####### bar graph
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erl_best_plot.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erl_cp_lung_only.RData")
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erl_cpp_lung_only.RData")

gg.cpp.auc.mean <- mean(cpp.best$AUC.test_auc[13:24])
gg.cpp.auc.se <- sd(cpp.best$AUC.test_auc[13:24])

gg.pp.auc.mean <- mean(pp.best$test_auc)
gg.pp.auc.se <- sd(pp.best$test_auc)

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_all_boxplot_slope.RData")
gg.cp_all.auc.best <- which(max(apply(cp_all.boxplot.slope, 2, mean)) == apply(cp_all.boxplot.slope, 2, mean))
gg.cp_all.auc.mean <- mean(cp_all.boxplot.slope[, gg.cp_all.auc.best ])
gg.cp_all.auc.se <- sd(cp_all.boxplot.slope[, gg.cp_all.auc.best ])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cpp_lung_var_lung_auc.RData")
gg.cpp_lung.auc.best <- which(max(apply(cpp.varying_training_matrix2, 2, mean)) == apply(cpp.varying_training_matrix2, 2, mean))
gg.cpp_lung.auc.mean <- mean(cpp.varying_training_matrix2[, gg.cpp_lung.auc.best])
gg.cpp_lung.auc.se <- sd(cpp.varying_training_matrix2[, gg.cpp_lung.auc.best])

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/cp_lung_boxplot_slope.RData")
gg.cp_lung.best <- which(max(apply(cp.lung_only.boxplot.slope , 2, mean)) == apply(cp.lung_only.boxplot.slope , 2, mean))
gg.cp_lung.auc.mean <- mean(cp.lung_only.boxplot.slope[, gg.cp_lung.best])
gg.cp_lung.auc.se <- sd(cp.lung_only.boxplot.slope[, gg.cp_lung.best])

gg.all <- data.frame(mean = c(gg.cpp.auc.mean, gg.pp.auc.mean, gg.cp_all.auc.mean, gg.cpp_lung.auc.mean, gg.cp_lung.auc.mean),
                     se = c(gg.cpp.auc.se, gg.pp.auc.se, gg.cp_all.auc.se, gg.cpp_lung.auc.se, gg.cp_lung.auc.se), 
                     Approach = c("CP2P", "P2P", "C2P", "CP2P Lung", "C2P Lung"),
                     Order = c("E", "C", "B", "D", "A"))

cols <- brewer.pal(n = 5, name = 'Set1')
p3 <- ggplot(gg.all, aes(x = Approach, y = mean, fill = Approach)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.9)) + theme_bw() +  
  scale_y_continuous(breaks = seq(0, 0.9, 0.05), "Test AUROC", limit = c(0, 0.9))
p3 <- p3 + theme(legend.direction = 'horizontal', 
                 #legend.position = 'top', 
                 #legend.position = c(0.5, 1.075),#'top', 
                 legend.position="none",
                 #plot.margin = unit(c(5.1, 0, 4.5, 1)/2, "lines"), 
                 plot.margin = unit(c(10, 7, 2, 3.5)/2, "lines"), 
                 #plot.margin = unit(c(15, 7, 2, 3.5)/2, "lines"), 
                 text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2)) + 
  scale_fill_manual(name="Approach", 
                    labels = c("C2P Lung", "C2P", "P2P", "CP2P Lung", "CP2P"),                    
                    values = cols) + 
  scale_x_discrete(name ="Approach", labels=c("C2P Lung", "C2P", "P2P", "CP2P Lung", "CP2P"))
p3

##### compare best
library("ggplot2")
library("grid")
library("RColorBrewer")

load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Erlotinib/WS/erl_best_plot.RData")
cl <- brewer.pal(n = 5, name = 'Dark2')

#cpp.ic50.m <- data.frame(xval = cpp.best$x_values, yval = cpp.best$IC50.test_auc)
cpp.auc.m <- data.frame(xval = cpp.best$x_values, yval = cpp.best$AUC.test_auc)
cpp.auc.m <- cpp.auc.m[13:24, ]
cpp.slope.m <- data.frame(xval = cpp.best$x_values, yval = cpp.best$slope.test_auc)
cpp.slope.m <- cpp.slope.m[13:24, ]
pp.m <- data.frame(xval = pp.best$x_values, yval = pp.best$test_auc)

all.m <- rbind(cpp.auc.m, cpp.slope.m, pp.m)
all.m$Approach <- c(
                    rep(cpp.best$AUC.name, dim(cpp.auc.m)[1]),
                    rep(cpp.best$slope.name, dim(cpp.slope.m)[1]),
                    rep(pp.best$name, dim(pp.m)[1]))

all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
all.m$Approach[nrow(all.m)] <- cp.best$AUC.name
all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
all.m$Approach[nrow(all.m)] <- cp.best$slope.name

all.m$Approach <- factor(all.m$Approach, levels = c(cp.best$AUC.name, cp.best$slope.name,
                                                    pp.best$name, 
                                                    cpp.best$AUC.name, cpp.best$slope.name))

p4 <- ggplot(all.m, aes(x=xval, y=yval, colour = Approach, ymax = 0.95)) + theme_bw() + 
  geom_line() + 
  geom_point(aes(shape=Approach, colour = Approach), size = 4, na.rm = TRUE) +
  geom_hline(aes(yintercept = cp.best$slope, colour = cp.best$slope.name), show_guide = FALSE) +   
  geom_hline(aes(yintercept = cp.best$AUC, colour = cp.best$AUC.name), show_guide = FALSE) +   
  scale_color_manual(values = cl) + 
  scale_shape_manual(values = c("C2P-AUC SNF all"=NA, 
                                "C2P-Slope SVM lin mRMR1000"=NA,                                 
                                "P2P SVM rbf mRMR1000"=16,                                
                                "CP2P-AUC SVM rbf all"=17, 
                                "CP2P-Slope SVM lin all"=18)) +
  scale_y_continuous(breaks = seq(0.5, 0.95, 0.05), "Test AUROC") +
  scale_x_continuous(breaks = seq(13, 24, by = 1), "# Number of Patient Samples in Training")
p4
p4 <- p4 + theme(legend.direction = 'horizontal', 
                 legend.position = c(0.5, 1.12),#'top', 
                 #plot.margin = unit(c(5.1, 0, 4.5, 1)/2, "lines"),                  
                 plot.margin = unit(c(12, 7, 2, 3.5)/2, "lines"), 
                 text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2))
p4 <- p4 + guides(col=guide_legend(ncol=2))
p4



######
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

pdf("erl_plot2.pdf", width = 14, height = 14, onefile=FALSE)

par(mfrow=c(2, 2))
plot.new()
plot.new()
plot.new()

gAll <- gtable:::rbind_gtable(gT, gB, "first")
#gAll <- gtable_add_rows(g, unit(5,"cm"), pos=nrow(g1))
grid.arrange(gAll)

#grid.arrange(gT, gB, nrow = 2, ncol = 1)

mtext(expression(bold('A')), adj=-0.1, side=3, outer=F, cex=1.5, line = 40)  
mtext(expression(bold('C')), adj=-0.1, side=3, outer=F, cex=1.5, line = 5)  

plot.new()
mtext(expression(bold('B')), adj=-0.1, side=3, outer=F, cex=1.5, line = 40)  
mtext(expression(bold('D')), adj=-0.1, side=3, outer=F, cex=1.5, line = 5)  

dev.off()

