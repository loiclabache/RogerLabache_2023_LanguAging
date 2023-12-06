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

base="path"
database = file.path(paste0(base, "1_1_gradient.xlsx"))

library(readxl)
donnees = read_xlsx(database)

nb_roi = c(15:51)

packages = c("dplyr", "stringr", "numDeriv","gamm4","magrittr",
             "ggplot2","scales", "gridExtra")
sapply(packages, require, character.only = T)

ROI = colnames(donnees)[nb_roi]
subset.size = length(ROI)
nvtx = subset.size
knots = 6

vtxmat = donnees[,nb_roi] %>% as.data.frame()

db = donnees[,c(1,6,14,8,4)]
colnames(db) = c("fsid_base", "Age", "hemi", "Sex", "Site_Name")

if(length(table(donnees$Side)) == 2){db = db %>%
  mutate(scanner_demean = Site_Name,
         sex_demean=Sex-mean(Sex),
         scanner_demean = (scanner_demean-mean(scanner_demean)) /sd(scanner_demean)) %>%
  select(fsid_base,Age,hemi,sex_demean,scanner_demean)} else{
    db = db %>%
      mutate(scanner_demean = Site_Name,
             sex_demean=Sex-mean(Sex),
             scanner_demean = (scanner_demean-mean(scanner_demean)) /sd(scanner_demean)) %>%
      select(fsid_base,Age,sex_demean,scanner_demean)
  }

i=1
Yorig = as.matrix(vtxmat)
N = ceiling(nvtx/subset.size)
print(N)
if (i == N) {
  end = dim(Yorig)[2]
} else {
  end = subset.size
}
Y = Yorig[,1:end]

nn = 100
Opt = list()
AGE = db$Age
Opt$Age = seq(min(AGE), max(AGE), length.out = nn) 
if(length(table(donnees$Side)) == 2){
  Opt$fake.frame = data.frame("Age" = Opt$Age,
                              "hemi" =rep(1,nn),
                              "sex_demean" = rep(0,nn),
                              "scanner_demean" = rep(0,nn))
}else{
  Opt$fake.frame = data.frame("Age" = Opt$Age,
                              "sex_demean" = rep(0,nn),
                              "scanner_demean" = rep(0,nn))
}

if(length(table(donnees$Side)) == 2){
  pdat = rbind(Opt$fake.frame,
               Opt$fake.frame)
  pdat$hemi[1:100] = 0
} else{
  pdat = Opt$fake.frame
}

res_BetaEffect = matrix(NA, ncol = 4, nrow = end)
colnames(res_BetaEffect) = c("Age_Beta", "Age_pValue",
                             "ageXhemi_Beta", "ageXhemi_pValue")
rownames(res_BetaEffect) = colnames(Yorig)

res_fit = as.data.frame(matrix(data=NA, nrow=100, ncol = ((6*end)+1)))
colnames(res_fit) = c("age", 
                      paste0(rep(colnames(Yorig), each=6), c(" - Lfit", 
                                                             " - Lse",
                                                             " - Rfit",
                                                             " - Rse",
                                                             " - AsymFit",
                                                             " - AsymSe")))
compteur = c(2:7)

for (j in 1:end) {
  if (j == 1) {
    gamm.trajectories = list()
    ogamm.trajectories = list()
    RR = list()
    pb = txtProgressBar(min=1, max=end, style=3)
    ph = list()
  }
  setTxtProgressBar(pb,j)
  
  Y = data.frame((Yorig)[,j])
  names(Y) = "Y"
  dat = data.frame(Y,db)
  
  if(length(table(donnees$Side)) == 2){  
    gamm.trajectories[[j]] = gamm4(Y ~ s(Age, by = as.factor(hemi), k = knots) + 
                                     as.factor(hemi) + sex_demean + scanner_demean, 
                                   data = dat, random = ~ (1 |fsid_base))
  }else{
    gamm.trajectories[[j]] = gamm4(Y ~ s(Age, k = knots) + 
                                     sex_demean + scanner_demean, 
                                   data = dat)
  }
  
  gamm.sum = summary(gamm.trajectories[[j]]$gam)
  g = gamm.trajectories[[j]]$gam
  
  Xp = predict(g, newdata = pdat, type = "lpmatrix") 
  
  if(length(table(donnees$Side)) == 2){
    c1 = grepl("hemi\\)1", colnames(Xp))
    c2 = grepl("hemi\\)0", colnames(Xp))
    
    r1 = with(pdat, hemi == 1)
    r2 = with(pdat, hemi == 0)

    X = Xp[r1, ] - Xp[r2, ]
    X[, !(c1 | c2)] = 0
    X[, !grepl("s\\(", colnames(Xp))] = 0
  }else{
    X =Xp
  }
  
  dif = X %*% coef(g)
  se = sqrt(rowSums((X %*% vcov(g, unconditional =T)) * X))
  comp = data.frame("OptAge" = Opt$Age, dif, se)
  
  if(length(table(donnees$Side)) == 2){
    dat = mutate(dat,
                 ohemi = ifelse(dat$hemi == 1, "left", "right"),
                 ohemi = factor(ohemi, levels = c("left","right"),ordered = T))
    ogamm.trajectories[[j]] = gamm4(Y ~ as.factor(hemi) + s(Age) + s(Age, by = ohemi, k = knots) + sex_demean + scanner_demean, 
                                    data = dat, random = ~ (1 |fsid_base))
  }else{
    ogamm.trajectories[[j]] = gamm4(Y ~  s(Age, k = knots) + 
                                      sex_demean + scanner_demean, 
                                    data = dat)
  }
  ogamm.sum = summary(ogamm.trajectories[[j]]$gam)
  
  plotData <- list()
  trace(mgcv:::plot.gam, at = list(c(27, 1)),
        quote({
          message("assigning into globalenv()'s plotData...")
          plotData <<- pd
        }))
  
  mgcv::plot.gam(gamm.trajectories[[j]]$gam, seWithMean = TRUE, pages = 1)
  
  if(length(table(donnees$Side)) == 2){
    plotdf=data.frame(Age=plotData[[1]]$x,
                      Lfit=plotData[[2]]$fit,
                      Rfit=plotData[[1]]$fit,
                      xl=plotData[[1]]$xlim,
                      Lse=plotData[[2]]$se,
                      Rse=plotData[[1]]$se)
  }else{
    plotdf=data.frame(Age=plotData[[1]]$x,
                      Lfit=plotData[[1]]$fit,
                      xl=plotData[[1]]$xlim,
                      Lse=plotData[[1]]$se)
  }
  
  ph[[j]] = ggplot(plotdf) +
    {if(length(table(donnees$Side)) == 2) 
      geom_line(aes(Age,Lfit),col="#F94144",size=0.7)
      else
        geom_line(aes(Age,Lfit),col="#000000",size=0.7)
    } +
    {if(length(table(donnees$Side)) == 2) 
      geom_line(aes(Age,Rfit),col="#90BE6D",size=0.7) 
    } +
    {if(length(table(donnees$Side)) == 2) 
      geom_ribbon(aes(Age,ymin=Lfit-Lse,ymax = Lfit+Lse), alpha = 0.2,fill="#F94144")
      else
        geom_ribbon(aes(Age,ymin=Lfit-Lse,ymax = Lfit+Lse), alpha = 0.2,fill="#000000") 
    } +
    {if(length(table(donnees$Side)) == 2) geom_ribbon(aes(Age,ymin=Rfit-Rse,
                                                          ymax = Rfit+Rse), 
                                                      alpha = 0.2,fill="#90BE6D")} +
    {if(length(table(donnees$Side)) == 2) geom_line(data=comp,aes(Opt$Age,dif),
                                                    col="#277DA1",size=0.7)} +
    {if(length(table(donnees$Side)) == 2) geom_ribbon(data=comp,
                                                      aes(x=Opt$Age,
                                                          ymin = dif-se, 
                                                          ymax = dif+se),
                                                      alpha = 0.2,fill="#008571")} +
    geom_hline(aes(yintercept=0), linetype='dotted') +
    theme_classic() +
    scale_x_continuous(limits = c(min(plotdf$Age), max(plotdf$Age)),
                       breaks = scales::pretty_breaks(n = 10),
                       expand = c(1e-2, 1e-2)) +
    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 5),
      expand = c(0, 0)) +
    ggtitle(ROI[j]) +
    ylab(NULL) + xlab(NULL) 
  
  
  res_fit[,1] = plotdf$Age
  res_fit[,compteur[1]] = plotdf$Lfit
  res_fit[,compteur[2]] = plotdf$Lse
  res_fit[,compteur[3]] = plotdf$Rfit
  res_fit[,compteur[4]] = plotdf$Rse
  res_fit[,compteur[5]] = comp$dif
  res_fit[,compteur[6]] = comp$se
  compteur = compteur+6
  
  RR$edf = RR$edf %>% cbind(., ogamm.sum$edf[2])
  if(length(table(donnees$Side)) == 2){
    RR$Fval = RR$Fval %>% cbind(., ogamm.sum$s.table[[2,3]])
  }else{
    RR$Fval = RR$Fval %>% cbind(., ogamm.sum$s.table[[3]]) 
  }
  RR$p.spl = RR$p.spl %>% cbind(., -log10(ogamm.sum$s.pv[2]))
  RR$pp.spl = RR$pp.spl %>% cbind(., ogamm.sum$s.pv[2])
  
  RR$fit_valL = RR$fit_valL %>% cbind(.,
                                      predict.gam(gamm.trajectories[[j]]$gam, 
                                                  newdata = pdat[101:200,]))
  if(length(table(donnees$Side)) == 2){RR$s = RR$fit_valR %>% cbind(.,
                                                                    predict.gam(gamm.trajectories[[j]]$gam, 
                                                                                newdata = pdat[1:100,]))
  }
  if(length(table(donnees$Side)) == 2){RR$fit_valdiff = RR$fit_valdiff %>% cbind(., dif)
  }else{
    RR$fit_valdiff =  RR$fit_val 
  }
  RR$se = RR$se %>% cbind(., se)
  
  RR$hT = RR$hT %>% cbind(., gamm.sum$p.t[[2]])
  RR$hCoef = RR$hCoef %>% cbind(., gamm.sum$p.coeff[[2]])
  RR$hP = RR$hP %>% cbind(., gamm.sum$p.pv[[2]])
  RR$hPlog = RR$hPlog %>% cbind(., -log10(gamm.sum$p.pv[[2]]))
  
  if(length(table(donnees$Side)) == 2){
    res_BetaEffect[j, 1] = coef(summary(ogamm.trajectories[[j]]$mer))[5,1]
    res_BetaEffect[j, 2] = ogamm.sum$s.pv[1]
    res_BetaEffect[j, 3] = coef(summary(ogamm.trajectories[[j]]$mer))[6,1]
    res_BetaEffect[j, 4] = ogamm.sum$s.pv[2]
  }else{
    res_BetaEffect[j, 1] = coef(summary(ogamm.trajectories[[j]]$mer))[4,1]
    res_BetaEffect[j, 2] = ogamm.sum$s.pv
  }
  
  
}
names(ph) = ROI

resultats = list()
resultats[[1]] = res_BetaEffect
resultats[[2]] = ph

names(resultats) = c("BetaEffect", "graphes")


head(resultats$BetaEffect)
resultats$graphes[[1]]
RR[[(length(RR)+1)]] = dimnames(Yorig)[[2]]
names(RR)[(length(RR))] = "roi_order"


global_x_min <- Inf
global_x_max <- -Inf
global_y_min <- Inf
global_y_max <- -Inf
for (plot in resultats$graphes) {
  x_limits <- ggplot_build(plot)$layout$panel_params[[1]]$x.range
  y_limits <- ggplot_build(plot)$layout$panel_params[[1]]$y.range
  
  global_x_min <- min(global_x_min, x_limits[1])
  global_x_max <- max(global_x_max, x_limits[2])
  global_y_min <- min(global_y_min, y_limits[1])
  global_y_max <- max(global_y_max, y_limits[2])
}
cat("Global x-axis range: [", global_x_min, ", ", global_x_max, "]\n", sep = "")
cat("Global y-axis range: [", global_y_min, ", ", global_y_max, "]", sep = "")
updated_graphes <- list()
for (i in seq_along(resultats$graphes)) {
  plot <- resultats$graphes[[i]]
  updated_plot <- plot +
    scale_x_continuous(limits = c(global_x_min, global_x_max)) +
    scale_y_continuous(limits = c(global_y_min, global_y_max))
  updated_graphes[[i]] <- updated_plot
}
names(updated_graphes) = names(resultats$graphes)

roi_order = read.csv(paste0(base, "1_2_gradient_cluster_traj.txt"))
matched_indices <- match(names(updated_graphes), roi_order$Abbreviation)
sorted_graphes <- updated_graphes[order(matched_indices)]

customize_axes <- function(plot, row, col, ncol) {
  plot <- plot + theme(text = element_text(family = "Arial", size=8),
                       plot.title = element_text(hjust=0.5),
                       axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.title.y = element_blank(),
                       axis.line = element_line(colour = "black", size = 0.5))
  return(plot)
}
ncol <- 7
nplots <- length(sorted_graphes)
for (i in seq_along(sorted_graphes)) {
  row <- ceiling(i / ncol)
  col <- (i - 1) %% ncol + 1
  sorted_graphes[[i]] <- customize_axes(sorted_graphes[[i]], row, col, ncol)
}

matched_indices <- match(names(sorted_graphes), roi_order$Abbreviation[order(roi_order$Abbreviation)])
sorted_graphes <- sorted_graphes[order(matched_indices)]
wraped_plot_updated_sorted <- do.call("grid.arrange",
                                      c(sorted_graphes, 
                                        ncol = ncol))

cluster_0_graphes <- sorted_graphes[match(roi_order[roi_order$Cluster_assignment == 2,]$Abbreviation, 
                                          names(sorted_graphes))]
wraped_plot_cluster_0 <- do.call("grid.arrange",
                                 c(cluster_0_graphes, 
                                   ncol = 2))


# for further reading see Gavin Simpson's smooth term comparison procedure outlined here
# https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
# https://fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/