################################################################################
# Written by James M Roe, Ph.D.
# Center for Lifespan Changes in Brain and Cognition, Department of Psychology
# University of Oslo, Oslo, Norway
# November 12, 2020
#-------------------------------------------------------------------------------
# Adapted by Loïc Labache, Ph.D.
# Holmes Lab Department of Psychology - Yale University
# December 6, 2023
################################################################################

# open libraries
#..............#
packages <- c("dplyr", "stringr", "magrittr", "forcats")
lapply(packages, require, character.only = T)


# load trajectories within signficant clusters
#............................................#
resdir = "path"
fit_val = read.csv(file.path(resdir,"3_1_fit_trajectories.csv"))[,-1]
# asym
fit_asym = fit_val[,c(seq(1,dim(fit_val)[2],by=6)-1)[-1]]
colnames(fit_asym) = gsub("...AsymFit", "", colnames(fit_asym))
# left
fit_L = fit_val[,c(seq(1,dim(fit_val)[2],by=6)+1)[-38]]
colnames(fit_L) = gsub("...Lfit", "", colnames(fit_L))
# right
fit_R = fit_val[,c(seq(1,dim(fit_val)[2],by=6)+3)[-38]]
colnames(fit_R) = gsub("...Rfit", "", colnames(fit_R))
fit_val_tmp = fit_val = fit_asym



significance = read.csv(file.path(resdir,"3_2_gradient_significance.csv"))
roi_sig = significance[significance$significance!=0,]$X
length(roi_sig)

fit_val = fit_val[, colnames(fit_val) %in% roi_sig]
dim(fit_val)
fit_L = fit_L[, colnames(fit_L) %in% roi_sig]
dim(fit_L)
fit_R = fit_R[, colnames(fit_R) %in% roi_sig]
dim(fit_R)

end = dim(fit_val)[2]
mat.dist.fit = matrix(rep(0,end*end), nrow=end)

#compute dissimilaritiy matrix
# ...........................#
pb = txtProgressBar(min=1, max=end, style=3)
for(k in 1:end) {
  setTxtProgressBar(pb,k)
  for(j in 1:end){
    mat.dist.fit[k,j] <- sum( (fit_val[,k] - fit_val[,j]) ^2 ) #dissimilarity matrix based on sum of squares between difference trajectories
  }
}

################################################################################
################################################################################
hemieffect = read.csv(file.path(resdir, "3_3_mapHCoef.csv")) #hemisphere effect
names(hemieffect) = colnames(fit_val_tmp)
hemieffect = hemieffect[names(hemieffect) %in% roi_sig]

age <- fit_val$age

# load libraries
# .............#
packages <- c("dplyr", "stringr", "magrittr", "cluster","ggplot2", "basicPlotteR")
sapply(packages, require, character.only = T)

# compute silhouettes #
# ...................#
cl.sil = 2:7
sil.coef.dist = NULL
for (i in cl.sil) {
  print(i)
  
  # PAM clustering (2-7 solutions)
  cl.dist = pam(mat.dist.fit, k=i, diss=T) 
  
  # silhouette
  sil.dist = silhouette(cl.dist) 
  
  # extract silhouette coefficients 
  sil.coef.dist[i] = summary(sil.dist)$avg.width
}

sols = data.frame("cl"=c(2:max(cl.sil)),
                  "dist"=sil.coef.dist[-1])


#silhouette coefficients
ps=ggplot(sols) +
  geom_vline(xintercept = sols$cl[sols$dist==max(sols$dist)], 
             linetype=2,col="dark grey") +
  geom_point(aes(x=cl,y=dist),size=2,col="#00886e") +
  geom_line(aes(x=cl,y=dist),linewidth=1,col="#00886e") +
  ylim(c(0,0.8)) +
  xlim(c(2,max(cl.sil))) +
  theme_classic() +
  labs(y="Silhouette coefficient",
       x="N clusters") +
  theme(text = element_text(size=12),
        axis.title.x = element_text(vjust=1),
        plot.title = element_text(hjust=0.5))
ps

## plot cluster solutions ##
# .........................#
LM.dist = sols$cl[sols$dist==max(sols$dist)] 
saveord=0 # 1
col.cl.main = c('lightblue1','orange1')
col.cl.mean = c('lightblue4','orange4')
ymin = min(fit_val[,]+as.vector(t(hemieffect[])))- 0.01
ymax = max(fit_val[,]+as.vector(t(hemieffect[]))) + 0.01

for (jj in 1:length(LM.dist)) {
  ncl <- LM.dist[jj]
  print(paste(ncl,"cluster solution"))
  cluster = pam(mat.dist.fit, k=ncl, diss=T)
  cl.ord = cluster$clustering
  for (cl in 1:ncl) {
    plot(0, type='n', xlab='Age', ylab='Gradient 1 asymmetry (LH-RH)', 
         main=paste('N =', toString(length(which(cl.ord == cl)))), 
         ylim=c(ymin,ymax),
         xlim=c(min(age),max(age)))
    #plot vertex trajectories with added Hemisphere effect
    for (n in which(cl.ord == cl)) {
      lines(age,fit_val[,n] + as.vector(t(hemieffect[n])), 
            type = 'l', lty = 1, lwd = 6, col = col.cl.main[cl])
    }
    #mean trajectory with added Hemisphere effect
    lines(age,(rowMeans(fit_val[,which(cl.ord == cl)])) + 
            mean(as.vector(t(hemieffect[,which(cl.ord == cl)]))), 
          type = 'l', lty = 1, lwd = 6, col = col.cl.mean[cl])
    #symmetry line
    lines(age, rep(0,length(age)), 
          type = "l", lty = 2, lwd = 3, col="black")
    for (n in which(cl.ord == cl)) {
      position = runif(1, min = 1, max = 100)
      text(x=age[position],
           y=(fit_val[,n] + as.vector(t(hemieffect[n])))[position], 
           pos=2,
           labels=names(hemieffect)[n])
    }
  }
  if (saveord == 1) {
    write.csv(setNames(cl.ord, names(hemieffect)),
              file = file.path(resdir,
                               paste('s.order_clusters.fit', toString(ncl),'.txt', 
                                     sep = "")))
  }
}

################################################################################
################################################################################
# parametre PLot 
lab.size = 12
title.size = 14
line.thickness = 2
theme_perso = theme(text = element_text(family = "Arial"), 
                    plot.title = element_text(hjust=0.5),
                    axis.text.y = element_text(size=lab.size, colour = "black"),
                    axis.title.y = element_text(size=title.size,colour = "black"),
                    axis.text.x = element_text(size=lab.size,colour = "black"),
                    axis.title.x = element_text(size=title.size,colour = "black"),
                    axis.line = element_line(colour = "black", size = 0.5))
  
# Plot variation Asym
# .............#
#cluster1
meancurve1 = rowMeans(fit_val[,which(cl.ord == 1)])+mean(as.vector(t(hemieffect[,which(cl.ord == 1)])))
sdd1 = apply(fit_val[,which(cl.ord == 1)],1,sd)
res1 = data.frame(meancurve1, sdd1, age)
#cluster2
meancurve2 = rowMeans(fit_val[,which(cl.ord == 2)])+mean(as.vector(t(hemieffect[,which(cl.ord == 2)])))
sdd2 = apply(fit_val[,which(cl.ord == 2)],1,sd)
res2 = data.frame(meancurve2, sdd2, age)

intersect(meancurve1, meancurve2)
#plot means and SD's
res = data.frame(res1,res2) %>% mutate(zer=0)
pvar=ggplot(res) +
  geom_line(aes(x = age, y = meancurve2), col = col.cl.main[2],
            linewidth = line.thickness, alpha=0.4) +
  geom_line(aes(x = age, y = meancurve1), col = col.cl.main[1], 
            linewidth = line.thickness) +
  geom_ribbon(aes(x = age, ymin = meancurve2-sdd2, ymax = meancurve2+sdd2),
              alpha = 0.2, fill = col.cl.main[2]) +
  geom_ribbon(aes(x = age, ymin = meancurve1-sdd1, ymax = meancurve1+sdd1),
              alpha = 0.6, fill = col.cl.main[1]) +
  geom_line(aes(x = age, y = zer), col = "black",
            linewidth = (line.thickness/2), linetype = 3) +
  xlab("Age") +
  ylab("Mean (SD)") +
  theme_classic() + 
  theme_perso
# find the intersection point:
f1 <- approxfun(ggplot_build(pvar)$data[[1]]$x, 
                ggplot_build(pvar)$data[[1]]$y)
f2 <- approxfun(ggplot_build(pvar)$data[[2]]$x, 
                ggplot_build(pvar)$data[[2]]$y)
intersection_courbes = optimize(function(t0) abs(f1(t0) - f2(t0)), 
                                interval = range(ggplot_build(pvar)$data[[1]]$x))
print(intersection_courbes)
pvar = pvar + geom_vline(xintercept = intersection_courbes$minimum,
                         color = "darkgrey",
                         linewidth = (line.thickness/2), 
                         linetype = "solid")
pvar
#...............................................................................
# Plot of proportion change:
res_prop <- res %>%
  mutate(IncrementalChange_meancurve1 = (((meancurve1 - meancurve1[1]) / abs(meancurve1[1])) * 100),
         IncrementalChange_meancurve2 = (((meancurve2 - meancurve2[1]) / abs(meancurve2[1])) * 100))
ggplot(res_prop) +
  geom_line(aes(x = age,
                y = IncrementalChange_meancurve2),
            col = col.cl.main[2], 
            linewidth = line.thickness, 
            alpha=0.4) +
  geom_line(aes(x = age, 
                y = IncrementalChange_meancurve1),
            col = col.cl.main[1],
            linewidth = line.thickness) +
  xlab("Age") +
  ylab("Percentage change") +
  theme_classic() + 
  theme_perso +
  geom_vline(xintercept = intersection_courbes$minimum,
             color = "darkgrey",
             size = (line.thickness/2),
             linetype = "solid")


################################################################################
################################################################################
significance = read.csv(file.path(resdir,"3_4_LMN_hROIs.csv"))
significance = significance[,c(4,7,10)]
significance_proportion = significance
significance = significance[significance$Cluster_assignment!=0,]
dim(significance)

cluster_number = 1 # 1 or 2 
all_net = TRUE # TRUE or FALSE 
LMN_part = "L" # L or M
LroiFit = significance[which(significance$Cluster_assignment == cluster_number),]
if (all_net==FALSE){
  LroiFit = LroiFit[which(LroiFit$LMNpart == LMN_part),]
}

#===============================
# Plot histogramme proportion LMN par clsuter
significance_proportion$functNet = as.factor(ifelse(significance_proportion$LMNpart !="L", "LM", "LCORE"))
significance_proportion$cluFact = as.factor(significance_proportion$Cluster_assignment)
prop_col = col.cl.main
prop_col[3] = "#555555" # gray
prop_col = prop_col[c(3,1,2)]
significance_proportion_cond = significance_proportion %>%
  count(functNet, cluFact) %>%
  group_by(functNet) %>%
  mutate(lab = round(prop.table(n) * 100, 2))
# graph
prop_lmn = ggplot(significance_proportion_cond) +
  geom_col(alpha=0.5,
           aes(x = functNet, 
               y = lab,
               fill = cluFact),
           position = "dodge") +
  scale_fill_manual(values = prop_col) +
  geom_text(aes(functNet, lab, 
                label = sprintf("%2.1f", lab), 
                group = cluFact), 
            position = position_dodge(width = 0.9)) +
  theme_classic() + theme_perso + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),)
prop_lmn


significance_proportion_byCluster = significance_proportion %>%
  count(functNet, cluFact) %>%
  group_by(cluFact) %>%
  mutate(lab = round(prop.table(n) * 100, 2))
# graph 2
significance_proportion_byCluster$cluFact = factor(significance_proportion_byCluster$cluFact,
                                                   levels = c("1", "2", "0"))
prop_lmn_byCluster = ggplot(significance_proportion_byCluster) +
  geom_col(alpha=0.5,
           aes(x = functNet, 
               y = lab,
               fill = cluFact),
           position = "dodge") +
  scale_fill_manual(values = prop_col[c(2,3,1)]) +
  geom_text(aes(functNet, lab, 
                label = sprintf("%2.1f", lab), 
                group = cluFact), 
            position = position_dodge(width = 0.9)) +
  theme_classic() + theme_perso + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),)
prop_lmn_byCluster


#===============================
# Plot variation by hemisphere
# 1 = L, 0 = R 
if (dim(LroiFit)[1]==1){ 
  #cluster1
  meancurve1 = (fit_L[,colnames(fit_L) %in% LroiFit$Abbreviation])+(as.vector(t(hemieffect[,names(hemieffect) %in% LroiFit$Abbreviation])))
  sdd1 = rep(0, 100)
  res1 = data.frame(meancurve1, sdd1, age)
  #cluster2
  meancurve2 = (fit_R[,colnames(fit_R) %in% LroiFit$Abbreviation])
  sdd2 = sdd1
  res2 = data.frame(meancurve2, sdd2, age)
  #Asym
  meancurve3 = (fit_val[,colnames(fit_val) %in% LroiFit$Abbreviation])+(as.vector(t(hemieffect[,names(hemieffect) %in% LroiFit$Abbreviation])))
  sdd3 = sdd1
  res3 = data.frame(meancurve2, sdd2, age)
}else{
  #cluster1
  meancurve1 = rowMeans(fit_L[,colnames(fit_L) %in% LroiFit$Abbreviation])+mean(as.vector(t(hemieffect[,names(hemieffect) %in% LroiFit$Abbreviation])))
  sdd1 = apply(fit_L[,colnames(fit_L) %in% LroiFit$Abbreviation],1,sd)
  res1 = data.frame(meancurve1, sdd1, age)
  #cluster2
  meancurve2 = rowMeans(fit_R[,colnames(fit_R) %in% LroiFit$Abbreviation])
  sdd2 = apply(fit_R[,colnames(fit_R) %in% LroiFit$Abbreviation],1,sd)
  res2 = data.frame(meancurve2, sdd2, age)
  #Asym
  meancurve3 = rowMeans(fit_val[,colnames(fit_val) %in% LroiFit$Abbreviation])+mean(as.vector(t(hemieffect[,names(hemieffect) %in% LroiFit$Abbreviation])))
  sdd3 = apply(fit_val[,colnames(fit_val) %in% LroiFit$Abbreviation],1,sd)
  res3 = data.frame(meancurve2, sdd2, age)
}

#plot means and SD's
res = data.frame(res1,res2) %>% mutate(zer=0)
if(cluster_number==1){
  perso_col = c("#C1D8DF", "#6C7A80") # blue gradient LH / RH
  col.cl.main
}else{
  perso_col = c("#DFB568","#604B24") # orange gradient LH / RH
}
pvar_clust=ggplot(res) +
  geom_line(aes(x = age, y = meancurve2), col = perso_col[2], 
            linewidth = line.thickness, alpha=0.4) +
  geom_line(aes(x = age, y = meancurve1), col = perso_col[1],
            linewidth = line.thickness) +
  geom_ribbon(aes(x = age, ymin = meancurve2-sdd2, ymax = meancurve2+sdd2),
              alpha = 0.8, fill = perso_col[2]) +
  geom_ribbon(aes(x = age, ymin = meancurve1-sdd1, ymax = meancurve1+sdd1),
              alpha = 0.8, fill = perso_col[1]) +
  geom_line(aes(x = age, y = zer), col = "black", 
            linewidth = (line.thickness/2), linetype = 3) +
  xlab("Age") +
  ylab("Mean (SD)") +
  theme_classic() + 
  theme_perso +
  geom_vline(xintercept = intersection_courbes$minimum,
             color = "darkgrey", 
             linetype = "solid")
pvar_clust
