library(magrittr)

Fval = pp.spl = NULL

load("2_gamm_results.Rda")

Fval %<>% cbind(.,RR$Fval)
pp.spl %<>% cbind(.,RR$pp.spl)

#FDR correction
#Age x Hemi effects
FDRcorr = p.adjust(pp.spl,method = "BY")
table(FDRcorr < 0.05)

res_p = data.frame(roi = RR$roi_order,
                   p_correctd_ageXhemi = FDRcorr, 
                   p_correctd_ageXhemi_bool = ifelse(FDRcorr < 0.05, 1, 0),
                   F_stat_ageXhemi = t(Fval))

sig_ROI = res_p[res_p$p_correctd_ageXhemi_bool==1,]
rownames(sig_ROI) = c(1:dim(sig_ROI)[1])
sig_ROI