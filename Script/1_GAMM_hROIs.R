################################################################################
# Written by Loïc Labache, Ph.D.                                               #
# Holmes Lab, Department of Psychiatry - Rutgers University                    #
# June 17, 2024                                                                #
################################################################################

# Libraries.....................................................................
#...............................................................................
packages = c("here", "readxl", "dplyr", "gamm4", "progress",
             "tidyr", "stringr", "ggplot2")
lapply(packages, require, character.only = T)

# Data..........................................................................
#...............................................................................
path_folder = "/Users/loiclabache/Library/CloudStorage/Dropbox/Collaboration externe/Projet Aging - MB, ER, GD/Conference/OHBM_2023"
grad_data = read_xlsx(here(path_folder, "/Data/728sujets_LMN_firstGradient.xlsx"))

# Generalized Additive Mixed Models.............................................
#...............................................................................
# Data Wrangling................................................................
gamm_data = grad_data[, c(1,6,14,8,4, 15:51)]
colnames(gamm_data)[1:5] = c("fsid_base", "Age", "hemi", "Sex", "Site_Name")
gamm_data$scanner_zscored = (gamm_data$Site_Name - mean(gamm_data$Site_Name)) / sd(gamm_data$Site_Name)
gamm_data$sex_demean = gamm_data$Sex - mean(gamm_data$Sex)
gamm_data[, c(4:5)] = NULL
# gamm_data = gamm_data[gamm_data$scanner_zscored == "-0.378664155473584", ]

# GAMM..........................................................................
# Some GAMM Parameters..........................................................
roi_index = c(4:40)
roi_name = colnames(gamm_data)[roi_index]
knots = 6
seq_age = seq(min(gamm_data$Age), 
              max(gamm_data$Age), 
              length.out = 100)
# Data Frames Results...........................................................
beta_Hemi = setNames(as.data.frame(matrix(NA, 
                                          ncol = 3,
                                          nrow = length(roi_index))), 
                     c("tvalue_hemi", "coef_hemi",
                       "pValue_hemi"))
rownames(beta_Hemi) = roi_name
beta_HemiByAge = setNames(as.data.frame(matrix(NA, 
                                               ncol = 3,
                                               nrow = length(roi_index))), 
                          c("Fstat_ageBYhemi", "edf_ageBYhemi",
                            "pValue_ageBYhemi"))
rownames(beta_HemiByAge) = roi_name
fit_traj = setNames(as.data.frame(matrix(data = NA,
                                         nrow = length(seq_age),
                                         ncol = ((6 * length(roi_index)) + 1))),
                    c("age", paste0(rep(roi_name, each = 6),
                                    c(" - L - pred", " - L - SE", 
                                      " - R - pred", " - R - SE",
                                      " - Asym - pred", " - Asym - SE"))))
fit_traj$age = seq_age
# Computation of Predicted Trajectories for each Region.........................
pb = progress_bar$new(format = "Progress: [:bar] :percent eta: :eta",
                      total = length(roi_name))
for (r in 1:length(roi_name)){
  pb$tick()
  # GAMM Data...................................................................
  gamm_roi = gamm_data[, -roi_index[-r]]
  current_roi = colnames(gamm_data[, -roi_index[-r]])[4]
  index_roi = rownames(beta_Hemi) == current_roi
  colnames(gamm_roi)[4] = "Y"
  # GAMM Model..................................................................
  gamm_traj = gamm4(Y ~ s(Age, by = as.factor(hemi), k = knots) + 
                      as.factor(hemi) + sex_demean + scanner_zscored, 
                    data = gamm_roi,
                    random = ~ (1 | fsid_base))
  gamm_traj_sum = summary(gamm_traj$gam)
  # Main Hemisphere Effect......................................................
  beta_Hemi[index_roi, ]$tvalue_hemi = gamm_traj_sum$p.t[[2]]
  beta_Hemi[index_roi, ]$coef_hemi = gamm_traj_sum$p.coeff[[2]]
  beta_Hemi[index_roi, ]$pValue_hemi = gamm_traj_sum$p.pv[[2]]
  # GAMM Model with Interaction.................................................
  gamm_roi = mutate(gamm_roi,
                    Hemi = ifelse(gamm_roi$hemi == 1, 
                                  "left", "right"),
                    Hemi = factor(Hemi, levels = c("left", "right"), 
                                  ordered = T))
  gamm_traj_inter = gamm4(Y ~ as.factor(hemi) +
                            s(Age) + s(Age, by = Hemi, k = knots) +
                            sex_demean + scanner_zscored, 
                          data = gamm_roi,
                          random = ~ (1 |fsid_base))
  gamm_traj_inter_sum = summary(gamm_traj_inter$gam)
  # Age by Hemisphere Effect....................................................
  beta_HemiByAge[index_roi, ]$Fstat_ageBYhemi = gamm_traj_inter_sum$s.table[[2,3]]
  beta_HemiByAge[index_roi, ]$edf_ageBYhemi = gamm_traj_inter_sum$edf[2]
  beta_HemiByAge[index_roi, ]$pValue_ageBYhemi = gamm_traj_inter_sum$s.pv[[2]]
  # Computation of Predicted Trajectories.......................................
  model = gamm_traj$gam
  new_data = data.frame("Age" = c(seq_age, seq_age),
                        "hemi" = c(rep(0, 100), rep(1, 100)), # 1 for Left, 0 for Right
                        "sex_demean" = rep(0, 200),
                        "scanner_zscored" = rep(0, 200))
  pred_data = predict(model, newdata = new_data, type = "lpmatrix")
  index_col_left = grepl("hemi\\)1", colnames(pred_data))
  index_col_right = grepl("hemi\\)0", colnames(pred_data))
  index_row_left = with(new_data, hemi == 1)
  index_row_right = with(new_data, hemi == 0)
  # Compute Asymmetry...........................................................
  x_asym = pred_data[index_row_left, ] - pred_data[index_row_right, ]
  x_asym[, !(index_col_left | index_col_right)] = 0
  x_asym[, !grepl("s\\(", colnames(pred_data))] = 0
  pred_asym = x_asym %*% coef(model)
  se_asym = sqrt(rowSums((x_asym %*% vcov(model, unconditional = T)) * x_asym))
  fit_traj[, colnames(fit_traj) == paste0(current_roi, " - Asym - pred")] = pred_asym
  fit_traj[, colnames(fit_traj) == paste0(current_roi, " - Asym - SE")] = se_asym
  # Left Trajectory.............................................................
  x_left = pred_data[index_row_left, ]
  x_left[, !index_col_left] = 0
  x_left[, !grepl("s\\(", colnames(pred_data))] = 0
  pred_left = x_left %*% coef(model)
  se_left = sqrt(rowSums((x_left %*% vcov(model, unconditional = T)) * x_left))
  fit_traj[, colnames(fit_traj) == paste0(current_roi, " - L - pred")] = pred_left
  fit_traj[, colnames(fit_traj) == paste0(current_roi, " - L - SE")] = se_left
  # Right Trajectory............................................................
  x_right = pred_data[index_row_right, ]
  x_right[, !index_col_right] = 0
  x_right[, !grepl("s\\(", colnames(pred_data))] = 0
  pred_right = x_right %*% coef(model)
  se_right = sqrt(rowSums((x_right %*% vcov(model, unconditional = T)) * x_right))
  fit_traj[, colnames(fit_traj) == paste0(current_roi, " - R - pred")] = pred_right
  fit_traj[, colnames(fit_traj) == paste0(current_roi, " - R - SE")] = se_right
}

# Visualization of Trajectories.................................................
#...............................................................................
for (i in 1:dim(fit_traj)[2]){
  fit_traj[, i] = as.numeric(fit_traj[, i])
}
fit_traj_long = gather(fit_traj, key = "measurement", value = "value", -age)
split_measure = as.data.frame(do.call(rbind, 
                                      str_split(fit_traj_long$measurement,
                                                " - ")))
fit_traj_long$region = split_measure[, 1]
fit_traj_long$side = split_measure[, 2]
fit_traj_long$measure = split_measure[, 3]
fit_traj_long$measurement = NULL
# Wrangling Data for Visualization..............................................
pred_data = fit_traj_long %>% filter(measure == "pred") %>%
  mutate(side = factor(side, levels = c("L", "R", "Asym")))
se_data = fit_traj_long %>% filter(measure == "SE") %>%
  mutate(side = factor(side, levels = c("L", "R", "Asym")))
merged_data = merge(pred_data, se_data, 
                    by = c("age", "region", "side"),
                    suffixes = c("_pred", "_se"))
merged_data$measure_pred = merged_data$measure_se = NULL
# Visualization.................................................................
ggplot(merged_data, aes(x = age, y = value_pred, color = side)) +
  geom_line(linewidth = 0.7) +
  geom_ribbon(aes(ymin = value_pred - value_se, 
                  ymax = value_pred + value_se, 
                  fill = side), 
              alpha = 0.2,
              colour = NA) +
  geom_hline(aes(yintercept = 0), linetype='dotted') +
  scale_color_manual(values = c("L" = "#F94144",
                                "R" = "#90BE6D",
                                "Asym" = "#277DA1")) +
  scale_fill_manual(values = c("L" = "#F94144",
                               "R" = "#90BE6D",
                               "Asym" = "#277DA1")) +
  facet_wrap(~region, scales = "free", ncol = 7) +
  labs(title = "", x = "", y = "") +
  scale_x_continuous(limits = c(min(merged_data$age), max(merged_data$age)),
                     breaks = scales::pretty_breaks(n = 5),
                     expand = c(1e-2, 1e-2)) +
  scale_y_continuous(limits = c(min(merged_data$value_pred - merged_data$value_se),
                                max(merged_data$value_pred + merged_data$value_se)),
    breaks = scales::pretty_breaks(n = 5),
    expand = c(0, 0)) +
  theme_classic() +
  theme(strip.background = element_blank())

# Regions with a Significant Age by Hemisphere Effect...........................
#...............................................................................
beta_HemiByAge$sig_reg = p.adjust(beta_HemiByAge$pValue_ageBYhemi,
                                  method = "BY")
beta_HemiByAge$sig_reg_bin = beta_HemiByAge$sig_reg < 0.05

# Save Results..................................................................
#...............................................................................
write.csv(data.frame(ROI = rownames(beta_HemiByAge),
                     sig_roi = beta_HemiByAge$sig_reg_bin),
          file = here("Data", "effect_HemiByAge.csv"))
write.csv(fit_traj, 
          file = here("Data", "trajectories_hrois.csv"))
write.csv(beta_Hemi, 
          file = here("Data", "effect_Hemi.csv"))

#...............................................................................
# Further Reading...............................................................
# See Gavin Simpson's Smooth Term Comparison Procedure..........................
# https://fromthebottomoftheheap.net/2017/10/11/difference-splines-i/
# https://fromthebottomoftheheap.net/2017/12/14/difference-splines-ii/




# Dice Index....................................................................
all_db = read.csv(here("Data", "effect_HemiByAge.csv"))
camcan_db = read.csv(here("Data", "effect_HemiByAge_CamCAN.csv"))
intersection_db = sum(all_db$sig_reg_bin & camcan_db$sig_reg_bin)
sum_all_db = sum(all_db$sig_reg_bin)
sum_camcan_db = sum(camcan_db$sig_reg_bin)
dice_index = (2 * intersection_db) / (sum_all_db + sum_camcan_db)
# A Dice index of 0.92 means that there is a 92% overlap between
# the two sets, suggesting that the vectors are very similar. 
# This high value implies that most elements that are TRUE in one vector are 
# also TRUE in the other vector, and there are very few discrepancies 
# between them.
# Also curves are highly similar, minus the extrema, where the 2 others data 
# bases helped equilibriate the sampling. 