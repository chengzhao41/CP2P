library("sva")
library("RColorBrewer")
library("ggplot2")
library("doParallel")
library("gridExtra")
library("gridBase")
source('~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Script/preparing_data_helper.R')
source('~/Dropbox/SNF_DRUG_PROJECT/Script/comGENE.R')

getwd()
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Main Figure Scripts")
#pdf("bor_plot2.pdf", width = 14, height = 14)#, onefile=FALSE)
plot.new()

########## Figure A and B
load("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/WS/bortezomib_data.RData")
temp.data <- comGENE(scale(bortezomib$cgp_slope), scale(bortezomib$patient.combat))
mean(temp.data[[1]])
mean(temp.data[[2]])

# source('~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Script/preparing_data_helper.R')
# bortezomib$combined_slope.ComBat <- ComBat_combine(batch = bortezomib.labels$slope_combined.source,
#                                              label = bortezomib.labels$slope_combined,
#                                              input_data = rbind(temp.data[[1]], temp.data[[2]]), prior.plots = FALSE)

# 
p1 <- show_pca(input_data = rbind(temp.data[[1]], temp.data[[2]]) , label = bortezomib.labels$slope_combined.source, pca_line_plot = FALSE, print = FALSE)
# #p1 <- p1 + coord_fixed(ratio = 1) + geom_point(size = 0.01)
p1 <- p1 + theme(legend.direction = 'horizontal', 
               legend.position = 'top', plot.margin = unit(c(5.1,7,4.5,3.5)/2, "lines"), 
               text = element_text(size=15), axis.title.x=element_text(vjust=-1.5)) + 
  scale_color_discrete(name = 'Dataset')
p1 <- p1 + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17))))
p1 <- p1 + geom_point(aes(shape = factor(bortezomib.labels$slope_combined.source), colour = factor(bortezomib.labels$slope_combined.source)), show_guide = FALSE)

#print(p1)
p2 <- show_pca(input_data = bortezomib$combined_slope.sva, label = bortezomib.labels$slope_combined.source, pca_line_plot = FALSE, print = FALSE)
#p2 <- p2 + coord_fixed(ratio = 1) + geom_point(size = 0.01)
p2 <- p2 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top', plot.margin = unit(c(5.1,7,4.5,3.5)/2, "lines"), text = element_text(size=15), axis.title.x=element_text(vjust=-1.5)) + 
  scale_color_discrete(name = 'Dataset')
p2 <- p2 + guides(colour = guide_legend(override.aes = list(size=2.5, linetype=0, shape = c(16, 17))))
p2 <- p2 + geom_point(aes(shape = factor(bortezomib.labels$slope_combined.source), colour = factor(bortezomib.labels$slope_combined.source)), show_guide = FALSE)


#### Figure C
load(paste0("~/output_temp/Bortezomib/cpp_var/slope/bortezomib_cpp_1to100_slope_var13.RData")) 

# Boxplot Data ---------------------------------
cpp.boxplot_matrix <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
  c(
    cpp.snf.single.l1000[[temp.ind]]$test.auc
    , cpp.other_model.l1000[[temp.ind]]$elasticNet$test_auc
    , cpp.other_model.l1000[[temp.ind]]$lasso$test_auc
    , cpp.other_model.l1000[[temp.ind]]$ridge$test_auc
    , cpp.other_model.l1000[[temp.ind]]$rf$test_auc
    , cpp.other_model.l1000[[temp.ind]]$svmLinear$test_auc
    , cpp.other_model.l1000[[temp.ind]]$svmRadial$test_auc    
    
    , cpp.snf.single.all[[temp.ind]]$test.auc
    , cpp.other_model.all[[temp.ind]]$elasticNet$test_auc
    , cpp.other_model.all[[temp.ind]]$lasso$test_auc
    , cpp.other_model.all[[temp.ind]]$ridge$test_auc
    , cpp.other_model.all[[temp.ind]]$rf$test_auc
    , cpp.other_model.all[[temp.ind]]$svmLinear$test_auc
    , cpp.other_model.all[[temp.ind]]$svmRadial$test_auc
    
    , cpp.snf.single.mRMR1000[[temp.ind]]$test.auc
    , cpp.other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$lasso$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$ridge$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$rf$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
    , cpp.other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc    
  )  
}

cpp.boxplot_data <- data.frame(
  cpp.boxplot_matrix)

temp.names <- c(
  "SNF L1000"
  , "Elastic Net L1000" 
  , "Lasso L1000" 
  , "Ridge L1000" 
  , "RF L1000" 
  , "SVM lin L1000" 
  , "SVM rbf L1000"             
  
  , "SNF all"
  , "Elastic Net all" 
  , "Lasso all"           
  , "Ridge all"
  , "RF all" 
  , "SVM lin all"
  , "SVM rbf all"     
  
  , "SNF mRMR1000"
  , "Elastic Net mRMR1000"
  , "Lasso mRMR1000"
  , "Ridge mRMR1000"
  , "RF mRMR1000"
  , "SVM lin mRMR1000"
  , "SVM rbf mRMR1000"
)

temp.order.cp2p <- order(apply(cpp.boxplot_data, 2, median), decreasing = TRUE)

#####
load(paste0("~/output_temp/Bortezomib/pp_var/bortezomib_pp_1to100_13_var13.RData"))   
load(paste0("~/output_temp/Bortezomib/pp_var/bortezomib_pp_1to100_13_var13_l1000.RData"))   

# Boxplot Data ---------------------------------
pp.boxplot_matrix <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
  c(
    pp.snf.single.l1000[[temp.ind]]$test.auc    
    , pp.other_model.l1000[[temp.ind]]$elasticNet$test_auc
    , pp.other_model.l1000[[temp.ind]]$lasso$test_auc
    , pp.other_model.l1000[[temp.ind]]$ridge$test_auc
    , pp.other_model.l1000[[temp.ind]]$rf$test_auc
    , pp.other_model.l1000[[temp.ind]]$svmLinear$test_auc
    , pp.other_model.l1000[[temp.ind]]$svmRadial$test_auc
    
    , pp.snf.single.all[[temp.ind]]$test.auc
    , pp.other_model.all[[temp.ind]]$elasticNet$test_auc
    , pp.other_model.all[[temp.ind]]$lasso$test_auc
    , pp.other_model.all[[temp.ind]]$ridge$test_auc
    , pp.other_model.all[[temp.ind]]$rf$test_auc
    , pp.other_model.all[[temp.ind]]$svmLinear$test_auc
    , pp.other_model.all[[temp.ind]]$svmRadial$test_auc    
    
    , pp.snf.single.mRMR1000[[temp.ind]]$test.auc    
    , pp.other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
    , pp.other_model.mRMR1000[[temp.ind]]$lasso$test_auc
    , pp.other_model.mRMR1000[[temp.ind]]$ridge$test_auc
    , pp.other_model.mRMR1000[[temp.ind]]$rf$test_auc
    , pp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
    , pp.other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc    
  )  
}
pp.boxplot_data <- data.frame(
  pp.boxplot_matrix)

temp.order.p2p <- order(apply(pp.boxplot_data, 2, median), decreasing = TRUE)

#######################
load(paste0("~/output_temp/Bortezomib/cp_var/slope/bortezomib_cp_1to100_slope_var30.RData"))

# Boxplot Data ---------------------------------
cp.boxplot_matrix <- foreach (temp.ind = 1:100, .combine=rbind) %do% {
  c(
    cp.snf.single.l1000[[temp.ind]]$test.auc    
    , cp.other_model.l1000[[temp.ind]]$elasticNet$test_auc
    , cp.other_model.l1000[[temp.ind]]$lasso$test_auc
    , cp.other_model.l1000[[temp.ind]]$ridge$test_auc
    , cp.other_model.l1000[[temp.ind]]$rf$test_auc
    , cp.other_model.l1000[[temp.ind]]$svmLinear$test_auc
    , cp.other_model.l1000[[temp.ind]]$svmRadial$test_auc
    
    , cp.snf.single.all[[temp.ind]]$test.auc
    , cp.other_model.all[[temp.ind]]$elasticNet$test_auc
    , cp.other_model.all[[temp.ind]]$lasso$test_auc
    , cp.other_model.all[[temp.ind]]$ridge$test_auc
    , cp.other_model.all[[temp.ind]]$rf$test_auc
    , cp.other_model.all[[temp.ind]]$svmLinear$test_auc
    , cp.other_model.all[[temp.ind]]$svmRadial$test_auc    
    
    , cp.snf.single.mRMR1000[[temp.ind]]$test.auc    
    , cp.other_model.mRMR1000[[temp.ind]]$elasticNet$test_auc
    , cp.other_model.mRMR1000[[temp.ind]]$lasso$test_auc
    , cp.other_model.mRMR1000[[temp.ind]]$ridge$test_auc
    , cp.other_model.mRMR1000[[temp.ind]]$rf$test_auc
    , cp.other_model.mRMR1000[[temp.ind]]$svmLinear$test_auc
    , cp.other_model.mRMR1000[[temp.ind]]$svmRadial$test_auc    
  )  
}
cp.boxplot_data <- data.frame(
  cp.boxplot_matrix)

temp.order.c2p <- order(apply(cp.boxplot_data, 2, median), decreasing = TRUE)

all.order <- c(temp.order.cp2p[1:3], temp.order.p2p[1:3], temp.order.c2p[1:3])
all.order 
all.order <- unique(all.order)
all.order <- rev(all.order)
temp.names[all.order]

# for significant values
#wilcox.test(cpp.boxplot_matrix[, all.order], cp.boxplot_matrix[, all.order], paired = TRUE, alternative = "greater")
#wilcox.test(cpp.boxplot_matrix[, all.order], pp.boxplot_matrix[, all.order], paired = TRUE, alternative = "greater")

library(reshape2)
library(ggplot2)
cl <- brewer.pal(n = 3, name = 'Dark2')

temp.pp <- pp.boxplot_data[, all.order]
colnames(temp.pp) <- temp.names[all.order]
temp.pp <- melt(temp.pp)
#
temp.cp <- cp.boxplot_data[, all.order]
colnames(temp.cp) <- temp.names[all.order]
temp.cp <- melt(temp.cp)
#
temp.cpp <- cpp.boxplot_data[, all.order]
colnames(temp.cpp) <- temp.names[all.order]
temp.cpp <- melt(temp.cpp)

best.boxplot_data.m <- rbind(temp.pp, temp.cp, temp.cpp)

colnames(best.boxplot_data.m)[1] <- "model"
best.boxplot_data.m$model <- as.vector(best.boxplot_data.m$model)

#best.boxplot_data.m$Approach <- factor(c(rep("P2P", 500), rep("C2P", 500), rep("CP2P", 500))))
best.boxplot_data.m <- cbind(best.boxplot_data.m, Approach = factor(c(rep("B", 500), rep("A", 500), rep("C", 500))))

p3 <- ggplot(best.boxplot_data.m, aes(x = factor(model, levels = unique(model)), y = value))
#p + geom_boxplot()
p3 <- p3 + geom_boxplot(aes(fill = Approach)) + theme_bw() + theme(legend.direction = 'horizontal', 
                                                                   legend.position = 'top', 
                                                                   plot.margin = unit(c(5.1, 7, 4.5, 3.5)/2, "lines"),                      
                                                                   text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.text.x=element_text(angle = -15, hjust = 0)
                                                                   , axis.title.y=element_text(vjust=1.5) 
) + xlab("Models") + ylab("Test AUROC") + scale_fill_discrete(name="Approach", breaks = c("A", "B", "C"), labels = c("C2P", "P2P", "CP2P")) + 
  scale_fill_manual(values = brewer.pal(n = 3, name = 'Dark2'), breaks = c("A", "B", "C"), labels = c("C2P", "P2P", "CP2P")) + 
  scale_y_continuous(breaks = seq(0.45, 1, 0.05), "Test AUROC")

#p3
################
setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Bortezomib/WS")
load("bor_best_plot.RData")

####################### SLOPE
library("ggplot2")
library("grid")
library("RColorBrewer")

cpp.m <- data.frame(cbind(cpp.best$x_values, cpp.best$slope.mean, cpp.best$slope.sdev))
colnames(cpp.m) <- c("xval", "yval", "se")
cpp.m <- cpp.m[2:15, ]

pp.m <- cbind(pp.best$x_values, y_values = pp.best$mean, sdev = pp.best$sdev)
colnames(pp.m) <- c("xval", "yval", "se")

all.m <- rbind(cpp.m, pp.m)
all.m$Approach <- c(rep("CP2P", dim(cpp.m)[1]), rep("P2P", dim(pp.m)[1]))
all.m$Approach <- factor(all.m$Approach, levels = c("C2P", "P2P", "CP2P"))

all.m[nrow(all.m) + 1, ] <- all.m[nrow(all.m) + 1, ] # Creates a new row filled with NAs
all.m$Approach[nrow(all.m)] <- "C2P"

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(2) # move them .05 to the left and right

cols <- brewer.pal(n = 3, name = 'Dark2')

p4 <- ggplot(all.m, aes(x=xval, y=yval, colour = Approach, ymax = 0.95)) + theme_bw() + 
  geom_errorbar(aes(ymin= yval - se, ymax = yval + se), width=5, position=pd) + 
  geom_line(position=pd) + 
  geom_point(aes(shape=Approach, colour = Approach), size = 3, na.rm = TRUE) + 
  geom_hline(aes(yintercept = cp.best$slope, colour = "C2P")) + 
  
  scale_color_manual(values = c("C2P" = cols[1], "P2P" = cols[2], "CP2P" = cols[3])) + 
  scale_shape_manual(values = c("C2P" = NA, "P2P" = 16, "CP2P" = 17)) +
  scale_y_continuous(breaks = seq(0.4, 0.95, 0.05), "Test AUROC") +
  scale_x_continuous(breaks = seq(20, 150, by = 20), "# Number of Patient Samples in Training")
p4 <- p4 + theme(legend.direction = 'horizontal', 
                 legend.position = 'top', 
                 plot.margin = unit(c(5.1, 7, 4.5, 3.5)/2, "lines"), 
                 text = element_text(size=15), axis.title.x=element_text(vjust=-1.5), axis.title.y=element_text(vjust=2))   
p4

#
require(gtable)
require(gridExtra)
# #Extract Grobs
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)


# 
# #Bind the tables
gT <- gtable:::cbind_gtable(g1, g2, "first")
gB <- gtable:::cbind_gtable(g3, g4, "first")
#g <- gtable_add_cols(g, unit(5,"cm"), pos=nrow(g1))

setwd("~/Dropbox/Cell Line Drug Response Prediction Project/Ensemble Of Similarity Networks/Cheng/Main Figure Scripts")

pdf("bor_plot2.pdf", width = 14, height = 14, onefile=FALSE)


par(mfrow=c(2, 2))
plot.new()
plot.new()
plot.new()

gAll <- gtable:::rbind_gtable(gT, gB, "first")
grid.arrange(gAll)


mtext(expression(bold('A')), adj=-0.1, side=3, outer=F, cex=1.5, line = 42)  
mtext(expression(bold('C')), adj=-0.1, side=3, outer=F, cex=1.5, line = 2)  

plot.new()
mtext(expression(bold('B')), adj=-0.1, side=3, outer=F, cex=1.5, line = 42)  
mtext(expression(bold('D')), adj=-0.1, side=3, outer=F, cex=1.5, line = 2)  

dev.off()

# #g <- gtable_add_cols(g, unit(-1,"cm"), pos=nrow(g1))
# #draw

#grid.newpage()
#grid.arrange(gT, gB, ncol = 1, nrow = 2)
#grid.arrange(gB, ncol = 2, nrow = 1)



#dev.off()