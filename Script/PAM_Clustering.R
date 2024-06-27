################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# June 21, 2024                                                                #
################################################################################

# Libraries.....................................................................
#...............................................................................
packages = c("here", "dplyr", "stringr", "cluster")
             
             
             
             # "readxl", "gamm4", "progress",
             # "tidyr", "", "ggplot2")
lapply(packages, require, character.only = T)

# Data..........................................................................
#...............................................................................
desc_lum = read.csv(here("Atlas", "language_memory_atlas.txt"))
traj_data = read.csv(here("Data", "trajectories_hrois.csv"))[, -1]
sig_roi = read.csv(here("Data", "effect_HemiByAge.csv"))[, -1]
# beta_hemi = read.csv(here("Data", "effect_Hemi.csv"))
# colnames(beta_hemi)[1] = "region"

# Data Wrangling................................................................
#...............................................................................
# Selecting Asymmetrical Trajectories...........................................
index_columns = grep("Asym.*pred|pred.*Asym",
                     colnames(traj_data), 
                     value = TRUE)
traj_data = traj_data %>%
  select(age, all_of(index_columns))
# Selecting Region with a Significant Hemisphere by Age Effect..................
pat_sig_roi = str_c(sig_roi[sig_roi$sig_roi == TRUE, ]$ROI,
                    collapse = "|")
index_sig_roi = colnames(traj_data)[str_detect(colnames(traj_data), 
                                               pat_sig_roi)]
traj_data = traj_data %>%
  select(age, all_of(index_sig_roi))
# beta_hemi = beta_hemi[beta_hemi$region %in% 
#                         sig_roi[sig_roi$sig_roi == TRUE, ]$ROI, ]

# Partition Around Medoids Classification.......................................
#...............................................................................
# Computation Sum of Squares between Trajectories...............................
ssdist = as.matrix(dist(t(traj_data[, -1]), method = "euclidean"))^2
# Determination of the Number of Clusters..............................................
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



# Visualization Trajectories....................................................
# Cluster 1.....................................................................
meancurve1 = rowMeans(traj_data[, which(cl.ord == 1)]) + 
  mean(as.vector(t(hemieffect[, which(cl.ord == 1)])))
sdd1 = apply(fit_val[,which(cl.ord == 1)],1,sd)
res1 = data.frame(meancurve1, sdd1, age)
# Cluster 2.....................................................................
meancurve2 = rowMeans(fit_val[,which(cl.ord == 2)])+mean(as.vector(t(hemieffect[,which(cl.ord == 2)])))
sdd2 = apply(fit_val[,which(cl.ord == 2)],1,sd)
res2 = data.frame(meancurve2, sdd2, age)








col_cl_main = c('lightblue1','orange1')
col_cl_mean = c('lightblue4','orange4')
ymin = min(traj_data[, -1]) - 0.01
ymax = max(traj_data[, -1]) + 0.01
# ymin = min(traj_data[, -1] + as.vector(t(beta_hemi[, 3]))) - 0.01
# ymax = max(traj_data[, -1] + as.vector(t(beta_hemi[, 3]))) + 0.01
