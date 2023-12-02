library(magrittr)
base="path"

setwd(base)

RRs = list.files(pattern = "gamm.results")

i=1
j = RRs[i]

if (i == 1) {
  #set outputs
  edf = Fval = p.spl = pp.spl = NULL 
  fit_valL = fit_valR = fit_valdiff = NULL
  hT = hCoef = hP = hPlog = NULL 
}
load(j)
edf %<>% cbind(.,RR$edf)
Fval %<>% cbind(.,RR$Fval)
p.spl %<>% cbind(.,RR$p.spl)
pp.spl %<>% cbind(.,RR$pp.spl)
fit_valL %<>% cbind(.,RR$fit_valL)
fit_valR %<>% cbind(.,RR$fit_valR)
fit_valdiff %<>% cbind(.,RR$fit_valdiff)
hT %<>% cbind(.,RR$hT)
hCoef %<>% cbind(.,RR$hCoef)
hP %<>% cbind(.,RR$hP)
hPlog %<>% cbind(.,RR$hPlog)


#FDR correction
#Age x Hemi effects
FDRcorr = p.adjust(pp.spl,method = "BY")
table(FDRcorr < 0.05)

#FDR correction
#Hemi effects
tmpHp = p.adjust(hP,method = "BY")
table(tmpHp < 0.05)

res_p = data.frame(roi = RR$roi_order,
                   p_correctd_ageXhemi = FDRcorr, 
                   p_correctd_ageXhemi_bool = ifelse(FDRcorr < 0.05, 1, 0),
                   F_stat_ageXhemi = t(Fval), 
                   p_correctd_hemi = tmpHp,
                   p_correctd_hemi_bool = ifelse(tmpHp < 0.05, 1, 0))

sig_ROI = res_p[res_p$p_correctd_ageXhemi_bool==1,]
rownames(sig_ROI) = c(1:dim(sig_ROI)[1])
sig_ROI