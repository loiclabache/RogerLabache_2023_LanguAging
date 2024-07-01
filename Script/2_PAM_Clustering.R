################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# June 21, 2024                                                                #
################################################################################

# Libraries.....................................................................
#...............................................................................
packages = c("here", "dplyr", "stringr", "cluster", "ggplot2")
lapply(packages, require, character.only = T)

# Data..........................................................................
#...............................................................................
desc_lum = read.csv(here("Atlas", "language_memory_atlas.txt"))
traj_data_LRA = read.csv(here("Data", "trajectories_hrois.csv"))[, -1]
sig_roi = read.csv(here("Data", "effect_HemiByAge.csv"))[, -1]
beta_hemi = read.csv(here("Data", "effect_Hemi.csv"))
colnames(beta_hemi)[1] = "region"

# Data Wrangling................................................................
#...............................................................................
# Selecting Asymmetrical Trajectories...........................................
index_columns = grep("Asym.*pred|pred.*Asym",
                     colnames(traj_data_LRA), 
                     value = TRUE)
traj_data_all = traj_data_LRA %>%
  select(age, all_of(index_columns))
# Selecting Region with a Significant Hemisphere by Age Effect..................
pat_sig_roi = str_c(sig_roi[sig_roi$sig_roi == TRUE, ]$ROI,
                    collapse = "|")
index_sig_roi = colnames(traj_data_all)[str_detect(colnames(traj_data_all), 
                                               pat_sig_roi)]
traj_data = traj_data_all %>%
  select(age, all_of(index_sig_roi))
beta_hemi = beta_hemi[beta_hemi$region %in%
                        sig_roi[sig_roi$sig_roi == TRUE, ]$ROI, ]

# Partition Around Medoids Classification.......................................
#...............................................................................
# Computation Sum of Squares between Trajectories...............................
ssdist = as.matrix(dist(t(traj_data[, -1]), method = "euclidean"))^2
# Determination of the Number of Clusters.......................................
cl_sil = 2:7
sil_coef_dist = numeric(length(cl_sil))
for (i in seq_along(cl_sil)) {
  k = cl_sil[i]
  cl_dist = pam(ssdist, k = k, diss = TRUE)
  sil_dist = silhouette(cl_dist)
  sil_coef_dist[i] = summary(sil_dist)$avg.width
}
sols = data.frame("cl" = cl_sil, "dist" = sil_coef_dist)
nb_clust = sols$cl[sols$dist==max(sols$dist)] 
# Clustering....................................................................
pam_clust = pam(ssdist, k = nb_clust, diss = TRUE)
clusters = pam_clust$clustering

# Visualization Trajectories: Asymmetries.......................................
#...............................................................................
# Cluster 1.....................................................................
mean_beta_clusterOne = mean(beta_hemi[beta_hemi$region %in% 
                                        gsub("...Asym...pred", "", 
                                             names(clusters[clusters == 1])), ]$coef_hemi)
mean_clusterOne = rowMeans(traj_data[, colnames(traj_data) %in% 
                                       names(clusters[clusters == 1])]) 
mean_clusterOne = mean_clusterOne + mean_beta_clusterOne
sd_clusterOne = apply(traj_data[, colnames(traj_data) %in%
                                  names(clusters[clusters == 1])], 1, sd)
res_cOne = data.frame(mean_clusterOne, sd_clusterOne, age = traj_data$age)
# Cluster 2.....................................................................
mean_beta_clusterTwo = mean(beta_hemi[beta_hemi$region %in% 
                                        gsub("...Asym...pred", "", 
                                             names(clusters[clusters == 2])), ]$coef_hemi)
mean_clusterTwo = rowMeans(traj_data[, colnames(traj_data) %in% 
                                       names(clusters[clusters == 2])]) 
mean_clusterTwo = mean_clusterTwo + mean_beta_clusterTwo
sd_clusterTwo = apply(traj_data[, colnames(traj_data) %in%
                                  names(clusters[clusters == 2])], 1, sd)
res_cTwo = data.frame(mean_clusterTwo, sd_clusterTwo)
# Plot Parameters...............................................................
lab.size = 12
title.size = 14
line.thickness = 2
theme_perso = theme(text = element_text(family = "Arial"), 
                    plot.title = element_text(hjust=0.5),
                    axis.text.y = element_text(size=lab.size, colour = "black"),
                    axis.title.y = element_text(size=title.size,colour = "black"),
                    axis.text.x = element_text(size=lab.size,colour = "black"),
                    axis.title.x = element_text(size=title.size,colour = "black"),
                    axis.line = element_line(colour = "black", linewidth = 0.5))
col_cl_asym = c('lightblue1','orange1')
# Visualization.................................................................
res = data.frame(res_cOne, res_cTwo) %>% mutate(zero = 0)
traj_clusters = ggplot(res) +
  geom_line(aes(x = age, y = mean_clusterTwo), 
            col = col_cl_asym[2],
            linewidth = line.thickness, 
            alpha = 0.4) +
  geom_line(aes(x = age, y = mean_clusterOne), 
            col = col_cl_asym[1], 
            linewidth = line.thickness) +
  geom_ribbon(aes(x = age,
                  ymin = mean_clusterTwo - sd_clusterTwo,
                  ymax = mean_clusterTwo + sd_clusterTwo),
              alpha = 0.2, 
              fill = col_cl_asym[2]) +
  geom_ribbon(aes(x = age, 
                  ymin = mean_clusterOne - sd_clusterOne,
                  ymax = mean_clusterOne + sd_clusterOne),
              alpha = 0.6, 
              fill = col_cl_asym[1]) +
  geom_line(aes(x = age, y = zero), 
            col = "black",
            linewidth = (line.thickness/2), 
            linetype = 3) +
  xlab("Age") +
  ylab("Mean (SD)") +
  theme_classic() +
  theme_perso
# Intersection between Trajectories Curves.......................................
f1 = approxfun(ggplot_build(traj_clusters)$data[[1]]$x, 
               ggplot_build(traj_clusters)$data[[1]]$y)
f2 = approxfun(ggplot_build(traj_clusters)$data[[2]]$x, 
               ggplot_build(traj_clusters)$data[[2]]$y)
intersection_courbes = optimize(function(t0) abs(f1(t0) - f2(t0)), 
                                interval = range(ggplot_build(traj_clusters)$data[[1]]$x))
traj_clusters = traj_clusters + 
  geom_vline(xintercept = intersection_courbes$minimum,
             color = "darkgrey",
             linewidth = (line.thickness/2), 
             linetype = "solid")
traj_clusters

# Visualization Proportion of Clusters by Function..............................
#...............................................................................
# Clusters Proportion by Function...............................................
clusters_roi = data.frame(region = gsub("...Asym...pred", "", names(clusters)), 
                          Cluster = as.numeric(clusters))
cluster_desc = merge(desc_lum[desc_lum$Hemisphere == "L", ], 
                     clusters_roi, 
                     by.x = "abbreviation",
                     by.y = "region", 
                     all.x = TRUE)
cluster_desc$Cluster[is.na(cluster_desc$Cluster)] = 0
cluster_desc$Function = as.factor(ifelse(cluster_desc$Function !="L",
                                         "LM", "LCORE"))
cluster_desc$cluFact = factor(cluster_desc$Cluster,
                              levels = c("1", "2", "0"))
cluster_proportion = cluster_desc %>%
  count(Function, Cluster) %>%
  group_by(Cluster) %>%
  mutate(lab = round(prop.table(n) * 100, 2))
# Visualization.................................................................
prop_col = col_cl_asym
prop_col[3] = "#555555" # gray
prop_col = prop_col[c(3,1,2)]
cluster_proportion$lab = round(cluster_proportion$lab)
cluster_proportion$Cluster = factor(cluster_proportion$Cluster,
                                    levels = c("1", "2", "0"))
prop_lmn_byCluster = ggplot(cluster_proportion) +
  geom_col(alpha=0.5,
           aes(x = Function, 
               y = lab,
               fill = Cluster),
           position = "dodge") +
  scale_fill_manual(values = prop_col[c(2,3,1)]) +
  geom_text(aes(Function, lab, 
                label = sprintf("%2.1f", lab), 
                group = Cluster), 
            position = position_dodge(width = 0.9)) +
  theme_classic() +
  theme_perso + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(),)
prop_lmn_byCluster

# Visualization Trajectories: by Hemispheres....................................
#...............................................................................
# Selecting Left Trajectories with a Significant Hemisphere by Age Effect.......
idx_col_L = grep("L.*pred|pred.*L",
                 colnames(traj_data_LRA), 
                 value = TRUE)
traj_data_L = traj_data_LRA %>%
  select(age, all_of(idx_col_L))
idx_sig_L = colnames(traj_data_L)[str_detect(colnames(traj_data_L), 
                                             pat_sig_roi)]
traj_data_L = traj_data_L %>%
  select(age, all_of(idx_sig_L))
# Selecting Right Trajectories with a Significant Hemisphere by Age Effect......
idx_col_R = grep("R.*pred|pred.*R",
                 colnames(traj_data_LRA), 
                 value = TRUE)
traj_data_R = traj_data_LRA %>%
  select(age, all_of(idx_col_R))
idx_sig_R = colnames(traj_data_R)[str_detect(colnames(traj_data_R), 
                                             pat_sig_roi)]
traj_data_R = traj_data_R %>%
  select(age, all_of(idx_sig_R))
# Average Trajectories By Hemisphere Cluster 1..................................
# Left..........................................................................
mean_clusterOne_L = rowMeans(traj_data_L[, colnames(traj_data_L) %in% 
                                           gsub("Asym", "L", names(clusters[clusters == 1]))]) 
mean_clusterOne_L = mean_clusterOne_L + mean_beta_clusterOne
sd_clusterOne_L = apply(traj_data_L[, colnames(traj_data_L) %in%
                                      gsub("Asym", "L", names(clusters[clusters == 1]))], 1, sd)
res_cOne_L = data.frame(mean_clusterOne_L, sd_clusterOne_L, age = traj_data_L$age)
# Right.........................................................................
mean_clusterOne_R = rowMeans(traj_data_R[, colnames(traj_data_R) %in% 
                                           gsub("Asym", "R", names(clusters[clusters == 1]))])
sd_clusterOne_R = apply(traj_data_R[, colnames(traj_data_R) %in%
                                      gsub("Asym", "R", names(clusters[clusters == 1]))], 1, sd)
res_cOne_R = data.frame(mean_clusterOne_R, sd_clusterOne_R, age = traj_data_R$age)
res_cOne_LR = data.frame(res_cOne_L, res_cOne_R) %>% mutate(zero = 0)
# Average Trajectories By Hemisphere Cluster 2..................................
# Left..........................................................................
mean_clusterTwo_L = rowMeans(traj_data_L[, colnames(traj_data_L) %in% 
                                           gsub("Asym", "L", names(clusters[clusters == 2]))]) 
mean_clusterTwo_L = mean_clusterTwo_L + mean_beta_clusterTwo
sd_clusterTwo_L = apply(traj_data_L[, colnames(traj_data_L) %in%
                                      gsub("Asym", "L", names(clusters[clusters == 2]))], 1, sd)
res_cTwo_L = data.frame(mean_clusterTwo_L, sd_clusterTwo_L, age = traj_data_L$age)
# Right.........................................................................
mean_clusterTwo_R = rowMeans(traj_data_R[, colnames(traj_data_R) %in% 
                                           gsub("Asym", "R", names(clusters[clusters == 2]))])
sd_clusterTwo_R = apply(traj_data_R[, colnames(traj_data_R) %in%
                                      gsub("Asym", "R", names(clusters[clusters == 2]))], 1, sd)
res_cTwo_R = data.frame(mean_clusterTwo_R, sd_clusterTwo_R, age = traj_data_R$age)
res_cTwo_LR = data.frame(res_cTwo_L, res_cTwo_R) %>% mutate(zero = 0)
# Visualization Cluster 1.......................................................
col_hemi_cOne = c("#C1D8DF", "#6C7A80")
traj_cOne = ggplot(res_cOne_LR) +
  geom_line(aes(x = age, y = mean_clusterOne_L), 
            col = col_hemi_cOne[1],
            linewidth = line.thickness) +
  geom_line(aes(x = age, y = mean_clusterOne_R), 
            col = col_hemi_cOne[2], 
            linewidth = line.thickness, 
            alpha = 0.4) +
  geom_ribbon(aes(x = age,
                  ymin = mean_clusterOne_L - sd_clusterOne_L,
                  ymax = mean_clusterOne_L + sd_clusterOne_L),
              alpha = 0.8, 
              fill = col_hemi_cOne[1]) +
  geom_ribbon(aes(x = age, 
                  ymin = mean_clusterOne_R - sd_clusterOne_R,
                  ymax = mean_clusterOne_R + sd_clusterOne_R),
              alpha = 0.8, 
              fill = col_hemi_cOne[2]) +
  geom_line(aes(x = age, y = zero), 
            col = "black",
            linewidth = (line.thickness/2), 
            linetype = 3) +
  xlab("Age") +
  ylab("Mean (SD)") +
  theme_classic() +
  theme_perso
traj_cOne
# Visualization Cluster 2.......................................................
col_hemi_cTwo = c("#DFB568","#604B24")
traj_cTwo = ggplot(res_cTwo_LR) +
  geom_line(aes(x = age, y = mean_clusterTwo_L), 
            col = col_hemi_cTwo[1],
            linewidth = line.thickness) +
  geom_line(aes(x = age, y = mean_clusterTwo_R), 
            col = col_hemi_cTwo[2], 
            linewidth = line.thickness, 
            alpha = 0.4) +
  geom_ribbon(aes(x = age,
                  ymin = mean_clusterTwo_L - sd_clusterTwo_L,
                  ymax = mean_clusterTwo_L + sd_clusterTwo_L),
              alpha = 0.8, 
              fill = col_hemi_cTwo[1]) +
  geom_ribbon(aes(x = age, 
                  ymin = mean_clusterTwo_R - sd_clusterTwo_R,
                  ymax = mean_clusterTwo_R + sd_clusterTwo_R),
              alpha = 0.8, 
              fill = col_hemi_cTwo[2]) +
  geom_line(aes(x = age, y = zero), 
            col = "black",
            linewidth = (line.thickness/2), 
            linetype = 3) +
  xlab("Age") +
  ylab("Mean (SD)") +
  theme_classic() +
  theme_perso
traj_cTwo