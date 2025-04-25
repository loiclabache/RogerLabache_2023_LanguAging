################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# April 25, 2025                                                               #
################################################################################

# Libraries.....................................................................
#...............................................................................
packages <- c("readxl", "ggridges", "ggplot2", "tidyr", "pals", "dplyr",
              "acca", "corrplot", "factoextra", "PCAtools", "viridis",
              "broom", "here")
lapply(packages, require, character.only = T)

# User parameters...............................................................
#...............................................................................
cluster_choice = 1 # 1 or 2 
choix_PC = "PC1" # PC1 or PC2
asym_bool = TRUE # do we want to work on asymmetry directly? yes or no
grad_bool = TRUE # Do we want to work with gradient & VolNorm only?

# Data..........................................................................
#...............................................................................
resdir = '/MyProject'
behavior_data = read_excel(file.path(resdir, "Data/behaviorialData.xlsx"))
# Keep only relevant behavioral data............................................
behavior_data = behavior_data[,c(1:4, 6, 5, 14:15, 18, 20)]
dim(behavior_data)
colnames(behavior_data)
for(i in 6:10){
  print(table(is.na(behavior_data[,i])))
  behavior_data = behavior_data[!is.na(behavior_data[,i]),]
  behavior_data[,i] = as.numeric(as.data.frame(behavior_data[,i])[,1])
  if(i==7){
    behavior_data[,i] = max(behavior_data[,i]) - behavior_data[,i]
  }
  if(i==9){
    behavior_data[,i] = max(behavior_data[,i]) - behavior_data[,i]
  }
  behavior_data[,i] = scale(behavior_data[,i])
}

# Behavioral Data Visualization.................................................
#...............................................................................
# Histograms by cohorts.........................................................
behavior_data_hist = read_excel(file.path(resdir, "Data/behaviorialData.xlsx"))
behavior_data_hist = behavior_data_hist[,c(1:4, 6, 5, 14:15, 18, 20)]
bd = read_excel(file.path(resdir,"Data/gradientData.xlsx"))
bd = bd[bd$Sujet %in% behavior_data_hist$Sujet,]
behavior_data_hist = behavior_data_hist[behavior_data_hist$Sujet %in% bd$Sujet,]
dim(behavior_data_hist)
behavior_data_hist_long <- behavior_data_hist %>%
  pivot_longer(cols = c(Age, MMSE, TOT_Summary_ToT_ratio, SynSem_nValid, SynSem_RTmax, Picture__Primming_Summary_ACC_baseline_all),
               names_to = "variable",
               values_to = "value") %>%
  drop_na(value)
behavior_data_hist_long$Acquisition_site = factor(behavior_data_hist_long$Acquisition_site,
                                                  levels = c("Grenoble",
                                                             "Omaha", 
                                                             "CamCan"))
ggplot(behavior_data_hist_long, aes(x = value, 
                                    fill = Acquisition_site)) +
  geom_histogram(position = "stack", alpha = 0.90,
                 color = "black",
                 bins = 30) +
  facet_wrap(~ variable, ncol = 2, scales = "free") +
  scale_fill_viridis(discrete = TRUE,
                     option = "mako", 
                     direction = -1) +
  theme_classic()
# Histograms by age bins........................................................
age_bins = behavior_data_hist_long %>%
  filter(variable == "Age") %>%
  select(Sujet, age = value) %>%
  mutate(AgeGroup = cut(age,
                        breaks = seq(10, 100, by = 10),
                        include.lowest = TRUE,
                        right = FALSE,
                        labels = paste0(seq(10,90,10), "–", seq(19,99,10)))) %>%
  select(-age)
plot_data = behavior_data_hist_long %>%
  filter(variable %in% c("MMSE",
                         "Picture__Primming_Summary_ACC_baseline_all",
                         "SynSem_nValid", "SynSem_RTmax",
                         "TOT_Summary_ToT_ratio")) %>%
  left_join(age_bins, by = "Sujet")
ggplot(plot_data, aes(x = value, fill = AgeGroup)) +
  geom_histogram(position = "stack", bins = 30, color = "black", alpha = 0.9) +
  facet_wrap(~ variable, ncol = 2, scales = "free") +
  scale_fill_viridis_d(option = "plasma", direction = -1) +
  theme_classic()

# PCA...........................................................................
#...............................................................................
# Reading data..................................................................
brain_data = read_excel(file.path(resdir,"Data/gradientData.xlsx"))
brain_data = brain_data[brain_data$Sujet %in% behavior_data$Sujet,]
behavior_data = behavior_data[behavior_data$Sujet %in% brain_data$Sujet,]
# "Pivot" data..................................................................
L_brain = brain_data[brain_data$Side=="L",]
L_brain = L_brain[,c(15:51)]
R_brain = brain_data[brain_data$Side=="R",]
R_brain = R_brain[,c(15:51)]
if (asym_bool == FALSE){
  colnames(L_brain) = paste0(colnames(L_brain), "_L")
  colnames(R_brain) = paste0(colnames(R_brain), "_R")
  for(i in 1:dim(L_brain)[2]){
    L_brain[,i] = scale(L_brain[,i])
    R_brain[,i] = scale(R_brain[,i])
  }
  brain_data = cbind(L_brain, R_brain)
}else {
  brain_data = L_brain-R_brain
  for(i in 1:dim(L_brain)[2]){
    brain_data[,i] = scale(brain_data[,i])
  }
  brain_data$Sujet = behavior_data$Sujet 
}
# More brain data...............................................................
volume = read_excel(file.path(resdir,"Data/VolNorm.xlsx"))
hiic = read_excel(file.path(resdir,"Data/hiic.xlsx"))
ct = read_excel(file.path(resdir,"Data/CTmean.xlsx"))
ct = ct[ct$Sujet %in% behavior_data$Sujet,]
# get the same participants for all data set....................................
volume = volume[volume$Sujet %in% ct$Sujet,]
hiic = hiic[hiic$Sujet %in% ct$Sujet,]
brain_data = brain_data[brain_data$Sujet %in% ct$Sujet,]
brain_data$Sujet = NULL
behavior_data = behavior_data[behavior_data$Sujet %in% ct$Sujet,]
L_volume = volume[volume$Side=="L",]
L_volume = L_volume[,c(15:51)]
R_volume = volume[volume$Side=="R",]
R_volume = R_volume[,c(15:51)]
L_ct = volume[volume$Side=="L",]
L_ct = L_ct[,c(15:51)]
R_ct = ct[ct$Side=="R",]
R_ct = R_ct[,c(15:51)]
hiic = hiic[,c(15:51)]
if (asym_bool == FALSE){
  colnames(L_volume) = paste0(colnames(L_volume), "_L")
  colnames(R_volume) = paste0(colnames(R_volume), "_R")
  brain_volume = cbind(L_volume, R_volume)
  colnames(L_ct) = paste0(colnames(L_ct), "_L")
  colnames(R_ct) = paste0(colnames(R_ct), "_R")
  brain_ct = cbind(L_ct, R_ct)
}else{
  brain_volume = L_volume-R_volume
  brain_ct = L_ct-R_ct
  for(i in 1:dim(brain_volume)[2]){
    brain_volume[,i] = scale(brain_volume[,i])
    brain_ct[,i] = scale(brain_ct[,i])
  }
}
# Select significant ROI........................................................
significance = read.csv(file.path(resdir,"Data/LMN_hROIs.csv"))
significance = significance[significance$Cluster_assignment!=0,]
significance$RoiSide = ifelse(significance$Cluster_assignment==1, 
                              paste0(significance$Abbreviation, "_R"),
                              paste0(significance$Abbreviation, "_L"))
# One pCCA by cluster ..........................................................
significance = significance[significance$Cluster_assignment==cluster_choice,]
if (asym_bool == FALSE){
  brain_data = brain_data[,colnames(brain_data) %in% significance$RoiSide]
  brain_volume = brain_volume[,colnames(brain_volume) %in% significance$RoiSide]
  brain_ct = brain_ct[,colnames(brain_ct) %in% significance$RoiSide]
}else{
  brain_data = brain_data[,colnames(brain_data) %in% significance$Abbreviation]
  brain_volume = brain_volume[,colnames(brain_volume) %in% significance$Abbreviation]
  brain_ct = brain_ct[,colnames(brain_ct) %in% significance$Abbreviation]
}
hiic = hiic[,colnames(hiic) %in% significance$Abbreviation]
for(i in 1:dim(hiic)[2]){
  hiic[,i] = scale(hiic[,i])
}
colnames(brain_data) = paste0(colnames(brain_data), "_G1")
colnames(brain_volume) = paste0(colnames(brain_volume), "_VolNorm")
colnames(hiic) = paste0(colnames(hiic), "_hiic")
colnames(brain_ct) = paste0(colnames(brain_ct), "_CT")
if (grad_bool == TRUE){
  colnames(brain_data) = paste0(colnames(brain_data), " * ", 
                                significance$LMNpart, " * ",
                                significance$Lobe)
  colnames(brain_volume) = paste0(colnames(brain_volume), " * ", 
                                  significance$LMNpart, " * ",
                                  significance$Lobe)
  full_brain = cbind(brain_data, brain_volume)
}else{
  full_brain = cbind(brain_data, brain_volume, hiic, brain_ct)
}
# Select the columns with Age and MMSE as covariates............................
covariates = as.data.frame(behavior_data[,c("Sex", "Age", "MMSE")])
covariates$Sex = ifelse(covariates$Sex=="male", 1, 0)
# Reduce number of brain variable...............................................
pca_brain = prcomp(full_brain, center = F, scale. = F)
# Select PCs which explain X amount of variance.................................
selectedPCs = findElbowPoint(get_eig(pca_brain)$variance.percent)
print(get_eig(pca_brain))
brainPCs = pca_brain$x[,1:selectedPCs]
# Inspect loading's to see which scales load onto the selected PCs..............
loadings = as.data.frame(pca_brain$rotation[,1:selectedPCs] )
loadings$ROI = rownames(loadings)
# plot the loadings to interpret the main components............................
loadings_long = loadings %>%
  pivot_longer(!ROI, names_to = "PCs", values_to = "eigenVal") %>% 
  group_by(PCs) %>% arrange(PCs, eigenVal) 
loadings_long_order = loadings_long %>%
  do(tibble(al=levels(reorder(interaction(.$PCs, .$ROI, drop=TRUE), 
                              .$eigenVal)))) %>% 
  pull(al)
loadings_long = loadings_long %>%
  mutate(al=factor(interaction(PCs, ROI), levels=loadings_long_order)) 
ggplot(loadings_long[loadings_long$PCs=="PC1" | loadings_long$PCs=="PC2",],
       aes(y = al, x = eigenVal)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           aes(fill = eigenVal)) +
  geom_text(size=5,
            aes(label=gsub("_G1", " G1", 
                           gsub("_VolNorm"," Vol", 
                                gsub(" \\*.*", "", loadings_long[loadings_long$PCs=="PC1" | loadings_long$PCs=="PC2",]$ROI))))) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100]) +
  facet_wrap( ~ PCs, scales = "free_y") +
  theme_classic()
# same thing in circle for 2 first comp.........................................
biplot_data = loadings[,c(1:2,dim(loadings)[2])]
biplot_data$Var = gsub(" \\*.*", "", biplot_data$ROI)
biplot_data$Var = sapply(strsplit(biplot_data$Var, "_\\s*"), tail, 1)
biplot_data$ROI = gsub(" \\*.*", "", biplot_data$ROI)
biplot_data$ROI = gsub('(.*)_\\w+', '\\1', biplot_data$ROI)
lim_ax = round(max(c(abs(max(biplot_data[,1:2])), abs(min(biplot_data[,1:2])))),2)
biplot_col = c("#64B67D", "#B6649D")
ggplot(biplot_data, aes(x = PC1, y = PC2)) +
  xlab("PC1") +
  ylab("PC2") + 
  scale_y_continuous(limits = c(-lim_ax, lim_ax)) +
  scale_x_continuous(limits = c(-lim_ax, lim_ax)) + 
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_vline(xintercept = 0,  color = "darkgrey", linetype = "solid") +
  geom_segment(aes(xend = 0, yend = 0,
                   color = Var,
                   alpha = sqrt((PC1 - 0)^2 + (PC2 - 0)^2)),
               arrow = arrow(length = unit(0.025, "inches"),
                             type = "closed", ends="first")) +
  scale_color_manual(values = biplot_col) + 
  geom_text_repel(size = 2,
                  aes(label=ROI, 
                      alpha = sqrt((PC1 - 0)^2 + (PC2 - 0)^2))) + 
  theme_classic()

# CCA...........................................................................
#...............................................................................
set.seed(13)
mod <- cc(X  = brainPCs,
          Y  = behavior_data[,c(7:10)],
          Zx = as.matrix(covariates),
          Zy = as.matrix(covariates))
mod$cor 
mod$prop_expl_var
# Permuration testing on CCA to derive p-values for each cannonical variate pair
mod_p <- cc_inference(mod, B = 1000)
mod_p$p_values
# Visualize CCA.................................................................
# Plot 1 laodings : brain mode 1 vs initial variables (=PC).....................
loadings_mode_brain = as.data.frame(mod$corr$corr.X.xscores)
loadings_mode_brain$Var = rownames(loadings_mode_brain)
ggplot(loadings_mode_brain, 
                         aes(y = reorder(Var, Cx1), x = Cx1)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           aes(fill = Cx1)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100])+
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_text(x=0, size=3,
            aes(label=loadings_mode_brain$Var,
                alpha=abs(loadings_mode_brain$Cx1))) +
  xlim(-max(abs(c(min(loadings_mode_brain$Cx1),max(loadings_mode_brain$Cx1)))),
       max(abs(c(min(loadings_mode_brain$Cx1),max(loadings_mode_brain$Cx1)))))+
  theme_classic() + 
  xlab("Loadings")
# Plot 2 laodings : behavior mode 1 vs initial variables (=scores)..............
loadings_mode_behavior = as.data.frame(mod$corr$corr.Y.yscores)
loadings_mode_behavior$Var = c("Tip of the the tongue",
                               "Accuracy",
                               "Reaction Time",
                               "Naming")
ggplot(loadings_mode_behavior, 
       aes(y = reorder(Var, Cy1), x = Cy1)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           aes(fill = Cy1)) +
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100])+
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_text(x=0, size=3,
            aes(label=loadings_mode_behavior$Var,
                alpha=abs(loadings_mode_behavior$Cy1))) +
  xlim(-max(abs(c(min(loadings_mode_behavior$Cy1),max(loadings_mode_behavior$Cy1)))),
       max(abs(c(min(loadings_mode_behavior$Cy1),max(loadings_mode_behavior$Cy1)))))+
  theme_classic() +
  xlab("Loadings") +
  scale_y_discrete(labels = abbreviate)
# ............................
# Plot 3 : correlation between both significant mode: (1 dot = 1 participant)...
# colored by the strongest association with behavioral score....................
# and size of a dot = the strongest association with brain PC = asym de G1......
individuals_mode = data.frame(x = mod$scores$xscores[,1],
                              y = mod$scores$yscores[,1])
ggplot(individuals_mode, aes(x, y)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(size=1, aes(colour=behavior_data$Age)) +
  geom_smooth(method='gam', color = "black") +
  scale_colour_viridis(option = "C") +
  theme_classic() + 
  xlab("Brain mode") +
  ylab("Behavioral mode") +
  labs(color="Age")
