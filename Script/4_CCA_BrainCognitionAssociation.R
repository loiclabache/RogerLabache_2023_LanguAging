################################################################################
# Written by Loïc Labache, Ph.D.
# Holmes Lab Department of Psychology - Yale University
# December 6, 2023
################################################################################

# open libraries
#..............#
packages <- c("readxl", "ggridges", "ggplot2", "tidyr", "pals", "dplyr",
              "acca", "corrplot", "factoextra", "PCAtools", "viridis",
              "broom", "here")
lapply(packages, require, character.only = T)

cluster_choice = 1 # 1 or 2 # <.................................................
choix_PC = "PC1" # PC1 or PC2 <.................................................

#............................................#
resdir = "path"
behavior_data = read_excel(file.path(resdir,"4_1_participants_behavior.xlsx"))
# Keep only relevant behavioral data
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

# Reshape data into long format
behavior_long <- gather(behavior_data, variable, value, MMSE:Picture__Primming_Summary_ACC_baseline_all)
# Calculate the mean value of each variable
mean_values <- behavior_long %>%
  group_by(variable) %>%
  summarize(mean = mean(value))
# Join the mean values to the data
behavior_long <- behavior_long %>%
  left_join(mean_values, by = "variable")
# Create a vector of Combas-inspired colors to use
colors <- c("#f44336", "#2196f3", "#4caf50", "#9c27b0", "#ff9800")
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
                    axis.line = element_line(colour = "black", size = 0.5),
                    legend.position="none")
# Create a ridgeline plot
distribution_behav = ggplot(behavior_long,
                            aes(x = value, 
                                y = reorder(variable, -value),  
                                color = variable)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(size = 3, shape = 20, alpha=0.2) +
  geom_density_ridges(scale = 1.5, 
                      size = line.thickness, 
                      alpha = 0.5) +
  scale_color_manual(values = colors) +
  theme_classic() + 
  theme_perso
distribution_behav


#............................
# Histograms
behavior_data_hist = read_excel(file.path(resdir,"4_1_participants_behavior.xlsx"))
behavior_data_hist = behavior_data_hist[,c(1:4, 6, 5, 14:15, 18, 20)]
bd = read_excel(file.path(resdir,"1_1_gradient.xlsx"))
bd = bd[bd$Sujet %in% behavior_data_hist$Sujet,]
behavior_data_hist = behavior_data_hist[behavior_data_hist$Sujet %in% bd$Sujet,]
dim(behavior_data_hist)
behavior_data_hist_long <- behavior_data_hist %>%
  pivot_longer(cols = c(Age, MMSE, TOT_Summary_ToT_ratio, SynSem_nValid, SynSem_RTmax, Picture__Primming_Summary_ACC_baseline_all),
               names_to = "variable",
               values_to = "value") %>%
  drop_na(value)
p <- ggplot(behavior_data_hist_long, aes(x = value)) +
  geom_histogram(position = "identity", alpha = 0.25,
                 bins = 30, color = "#555555") +
  facet_wrap(~ variable, ncol = 2, scales = "free") +
  theme_classic() + 
  theme_perso
print(p)


#............................
# Behavior depending of age
behavior_long$facet = factor(behavior_long$variable,
                             levels = c("Picture__Primming_Summary_ACC_baseline_all",
                                        "TOT_Summary_ToT_ratio",
                                        "SynSem_RTmax",
                                        "SynSem_nValid",
                                        "MMSE"))
levels(behavior_long$facet) <- c("Naming", "Tip of the the tongue",
                                 "Reaction Time", "Accuracy", "MMSE")
cor_age_behav = ggplot(behavior_long[behavior_long$variable!="MMSE",],
                       aes(x = Age, 
                           y = value,
                           color = variable,
                           fill = variable)) +
  facet_wrap( ~ facet) +
  geom_point(size = 3, shape = 20, alpha=0.1, color="black") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_vline(xintercept = 52.5, color = "darkgrey", linetype = "solid") +
  geom_smooth(method="gam") +
  theme_classic() + 
  theme_perso +
  labs(y = "Normalized scores") +
  scale_color_viridis(discrete=TRUE, option = "C") +
  scale_fill_viridis(discrete=TRUE, option = "C") 
cor_age_behav



################################################################################
################################################################################
# pCCA
# <...................................................................................................................
# <...................................................................................................................
# do we want to work on asymetry directly? yes or no
asym_bool = TRUE 
# Do we want to work with gradient & VolNorm only?
grad_bool = TRUE
# <...................................................................................................................
# <...................................................................................................................
# Reading data
brain_data = read_excel(file.path(resdir,"1_1_gradient.xlsx"))
brain_data = brain_data[brain_data$Sujet %in% behavior_data$Sujet,]
dim(brain_data)[1]/2
behavior_data = behavior_data[behavior_data$Sujet %in% brain_data$Sujet,]
dim(behavior_data)[1]
# "Pivot" data 
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
#################
# more brain data 
volume = read_excel(file.path(resdir,"4_2_VolNorm.xlsx"))
volume = volume[volume$Sujet %in% behavior_data$Sujet,]
# get the same participants for all data set
brain_data = brain_data[brain_data$Sujet %in% volume$Sujet,]
brain_data$Sujet = NULL
behavior_data = behavior_data[behavior_data$Sujet %in% volume$Sujet,]
dim(volume)[1]/2
dim(brain_data)[1]
dim(behavior_data)[1]
#-
L_volume = volume[volume$Side=="L",]
L_volume = L_volume[,c(15:51)]
R_volume = volume[volume$Side=="R",]
R_volume = R_volume[,c(15:51)]
if (asym_bool == FALSE){
  colnames(L_volume) = paste0(colnames(L_volume), "_L")
  colnames(R_volume) = paste0(colnames(R_volume), "_R")
  brain_volume = cbind(L_volume, R_volume)
}else{
  brain_volume = L_volume-R_volume
  for(i in 1:dim(brain_volume)[2]){
    brain_volume[,i] = scale(brain_volume[,i])
  }
}
# Select significant ROI
significance = read.csv(file.path(resdir,"3_4_LMN_hROIs.csv"))
significance = significance[significance$Cluster_assignment!=0,]
significance$RoiSide = ifelse(significance$Cluster_assignment==1, 
                              paste0(significance$Abbreviation, "_R"),
                              paste0(significance$Abbreviation, "_L"))
################################################################################
# One pCCA by cluster 
significance = significance[significance$Cluster_assignment==cluster_choice,]
#####
if (asym_bool == FALSE){
  brain_data = brain_data[,colnames(brain_data) %in% significance$RoiSide]
  brain_volume = brain_volume[,colnames(brain_volume) %in% significance$RoiSide]
}else{
  brain_data = brain_data[,colnames(brain_data) %in% significance$Abbreviation]
  brain_volume = brain_volume[,colnames(brain_volume) %in% significance$Abbreviation]
}
# rename columns of everything
colnames(brain_data) = paste0(colnames(brain_data), "_G1")
colnames(brain_volume) = paste0(colnames(brain_volume), "_VolNorm")
# selection of brain variables to inject into the PCA:
if (grad_bool == TRUE){
  colnames(brain_data) = paste0(colnames(brain_data), " * ", 
                                significance$LMNpart, " * ",
                                significance$Lobe)
  colnames(brain_volume) = paste0(colnames(brain_volume), " * ", 
                                  significance$LMNpart, " * ",
                                  significance$Lobe)
  full_brain = cbind(brain_data, brain_volume)
}else{
  full_brain = cbind(brain_data, brain_volume)
}
dim(full_brain)
# Select the columns with Age and MMSE as covariates
covariates = as.data.frame(behavior_data[,c("Sex", "Age", "MMSE")])
covariates$Sex = ifelse(covariates$Sex=="male", 1, 0)


#===================================
# Reduce number of brain variable 
pca_brain = prcomp(full_brain, center = F, scale. = F)
# Check scree plot (how much variance explained by each PC)
fviz_eig(pca_brain)
# Select PCs which explain X amount of variance 
selectedPCs = findElbowPoint(get_eig(pca_brain)$variance.percent)
print(get_eig(pca_brain))
brainPCs = pca_brain$x[,1:selectedPCs]
# Inspect loading's to see which scales load onto the selected PCs
loadings = as.data.frame(pca_brain$rotation[,1:selectedPCs] )
loadings$ROI = rownames(loadings)
# plot the loadings to interpret the main components
loadings_long = loadings %>%
  pivot_longer(!ROI, names_to = "PCs", values_to = "eigenVal") %>% 
  group_by(PCs) %>% arrange(PCs, eigenVal) 
loadings_long_order = loadings_long %>%
  do(tibble(al=levels(reorder(interaction(.$PCs, .$ROI, drop=TRUE), 
                              .$eigenVal)))) %>% 
  pull(al)
loadings_long = loadings_long %>%
  mutate(al=factor(interaction(PCs, ROI), levels=loadings_long_order)) 
# plot the main components
loadings_pca = ggplot(loadings_long[loadings_long$PCs=="PC1" | loadings_long$PCs=="PC2",],
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
                       high = coolwarm(100)[100])+
  facet_wrap( ~ PCs, scales = "free_y") +
  theme_classic() + theme_perso +
  theme(strip.text = element_text(size=16),
        strip.background = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) 
loadings_pca

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# same thing in circle for 2 first comp
biplot_data = loadings[,c(1:2,dim(loadings)[2])]
biplot_data$Var = gsub(" \\*.*", "", biplot_data$ROI)
biplot_data$Var = sapply(strsplit(biplot_data$Var, "_\\s*"), tail, 1)
biplot_data$ROI = gsub(" \\*.*", "", biplot_data$ROI)
biplot_data$ROI = gsub('(.*)_\\w+', '\\1', biplot_data$ROI)
lim_ax = round(max(c(abs(max(biplot_data[,1:2])), abs(min(biplot_data[,1:2])))),2)
biplot_col = c("#64B67D", "#B6649D") # vert et rose
biplot_pca = ggplot(biplot_data, aes(x = PC1, y = PC2)) +
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
  theme_classic() + theme_perso + 
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5))
biplot_pca


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Visualization of the relationship between the 2 first component
main_comp_ind = data.frame(PC1 = pca_brain$x[,1],
                           PC2 = pca_brain$x[,2])
relation_comp = ggplot(main_comp_ind, aes(x = PC1, y = PC2)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(aes(size=2,
                 colour=covariates$Age)) +
  scale_colour_viridis(option = "E") +
  theme_classic() + 
  xlab("PC1") +
  ylab("PC2") +
  labs(color="Age") +
  guides(size = "none") +
  theme(text = element_text(family = "Arial"))
relation_comp


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Visualization of the relationship between a component and original var
# <...................................................................................................................
# <...................................................................................................................
PC_by_brain = full_brain
PC_by_brain$PC1 = main_comp_ind$PC1
PC_by_brain$PC2 = main_comp_ind$PC2
PC_by_brain_long <- PC_by_brain %>%
  gather(key = "Brain_Region", value = "Value", -PC1, -PC2)
PC_by_brain_long_split <- PC_by_brain_long %>%
  separate(Brain_Region, into = c("Region", "Hemisphere", "Lobe"),
           sep = " \\* ", remove = TRUE)
PC_by_brain_long_split$Hemisphere=NULL
PC_by_brain_long_split = PC_by_brain_long_split %>%
  mutate(Metric = ifelse(grepl("G1", Region), "G1", 
                         ifelse(grepl("VolNorm", Region), "VolNorm", NA)))
PC_by_brain_long_split$Region = gsub("_G1", "", PC_by_brain_long_split$Region)
PC_by_brain_long_split$Region = gsub("_VolNorm", "", PC_by_brain_long_split$Region)
PC_by_brain_long_split$RegionMetric = paste(PC_by_brain_long_split$Region,
                                            PC_by_brain_long_split$Metric,
                                            sep=" | ")
label_df = PC_by_brain_long_split %>% group_by(Metric, Region) %>% 
  do(augment(lm(as.formula(paste("Value ~", choix_PC)), data=.))) %>% 
  top_n(1, get(choix_PC)) %>%
  select(Metric, Region, all_of(choix_PC), .fitted) %>% 
  rename(Value = .fitted, ROI = Region)
pal_col_binaire = c("#159090", ## blue
                    "#E7298A") ## pink
pc_varBrain = ggplot(PC_by_brain_long_split,
                     aes(x = get(choix_PC), y = Value,
                         color = factor(Metric))) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(aes(size=2), alpha = 0.1) + 
  geom_text_repel(data = label_df,
                  max.overlaps =13,
                  aes(label = ROI),
                  size = 3,
                  nudge_x = 1) + 
  geom_smooth(method="lm", se = FALSE,
              aes(group = RegionMetric)) +
  scale_color_manual(values = pal_col_binaire) +
  theme_classic() +
  xlab(choix_PC) +
  ylab("Brain Var") +
  guides(alpha = "none", size = "none") +
  theme(text = element_text(family = "Arial", size = 16),
        legend.title = element_blank())
pc_varBrain


#==========================================
# Compute CCA
set.seed(13)
mod <- cc(X  = brainPCs,
          Y  = behavior_data[,c(7:10)],
          Zx = as.matrix(covariates),
          Zy = as.matrix(covariates))
mod$cor #print R (canonical correlation coef between each cca mode, which is our primary test-stat)
mod$prop_expl_var # Variance explained by each mode of the CCA
# Permuration testing on CCA to derive p-values for each cannonical variate pair
mod_p <- cc_inference(mod, 
                      B = 1000) # B= number of perms 
mod_p$p_values
cormat = cor(cbind(brainPCs, 
                   behavior_data[,c(7:10)],
                   as.matrix(covariates)))
corrplot(cormat,
         method = 'color', 
         col = coolwarm(100),
         diag=FALSE, tl.cex=0.7)
#===============
# Visualize CCA  
# Plot 1 laodings : brain mode 1 vs initial variables (=PC)
# loadings are correlations between variables and the canonical variates.
loadings_mode_brain = as.data.frame(mod$corr$corr.X.xscores)
loadings_mode_brain$Var = rownames(loadings_mode_brain)
brain_mode_plot = ggplot(loadings_mode_brain, 
                         aes(y = reorder(Var, Cx1), x = Cx1)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           aes(fill = Cx1)) + # Border color
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100])+
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_text(x=0, size=3,
            aes(label=loadings_mode_brain$Var,
                alpha=abs(loadings_mode_brain$Cx1))) +
  xlim(-max(abs(c(min(loadings_mode_brain$Cx1),max(loadings_mode_brain$Cx1)))),
       max(abs(c(min(loadings_mode_brain$Cx1),max(loadings_mode_brain$Cx1)))))+
  theme_classic() + theme_perso +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Loadings")
brain_mode_plot

#............................
# Plot 2 laodings : behavior mode 1 vs initial variables (=scores)
loadings_mode_behavior = as.data.frame(mod$corr$corr.Y.yscores)
loadings_mode_behavior$Var = c("Tip of the the tongue",
                               "Accuracy",
                               "Reaction Time",
                               "Naming")
behav_mode_plot = ggplot(loadings_mode_behavior, 
                         aes(y = reorder(Var, Cy1), x = Cy1)) +
  geom_bar(stat = "identity",
           show.legend = FALSE,
           aes(fill = Cy1)) + # Border color
  scale_fill_gradient2(low = coolwarm(100)[1],
                       mid = coolwarm(100)[50],
                       high = coolwarm(100)[100])+
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_text(x=0, size=3,
            aes(label=loadings_mode_behavior$Var,
                alpha=abs(loadings_mode_behavior$Cy1))) +
  xlim(-max(abs(c(min(loadings_mode_behavior$Cy1),max(loadings_mode_behavior$Cy1)))),
       max(abs(c(min(loadings_mode_behavior$Cy1),max(loadings_mode_behavior$Cy1)))))+
  theme_classic() + theme_perso +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank()) +
  xlab("Loadings") +
  scale_y_discrete(labels = abbreviate)
behav_mode_plot

# ............................
# Plot 3 : correlation between both significant mode: (1 dot = 1 participant)
# colored by the strongest association with behavioral score 
# and size of a dot = the strongest association with brain PC = asym de G1 
individuals_mode = data.frame(x = mod$scores$xscores[,1],
                              y = mod$scores$yscores[,1])
mod_cca = ggplot(individuals_mode, aes(x, y)) +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_hline(yintercept = 0, color = "darkgrey", linetype = "solid") +
  geom_point(size=1, aes(colour=behavior_data$Age)) +
  geom_smooth(method='gam', color = "black") +
  scale_colour_viridis(option = "C") +
  theme_classic() + 
  xlab("Brain mode") +
  ylab("Behavioral mode") +
  labs(color="Age") + 
  theme_perso
mod_cca

