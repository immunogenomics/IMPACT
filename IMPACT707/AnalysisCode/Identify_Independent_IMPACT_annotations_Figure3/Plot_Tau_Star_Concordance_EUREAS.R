length(traits_EUR) 
length(traits_EAS) 
pops <- c("EUR","EAS")

#compute from other scripts
mytau_star_p_fdr <- cbind(mytau_star_p_EUR_IMPACT,mytau_star_p_EAS_IMPACT)
mytau_star <- cbind(mytau_star_EUR_IMPACT,mytau_star_EAS_IMPACT)
a <- mytau_star_p_fdr
a[mytau_star < 0] <- log(x = mytau_star_p_fdr[mytau_star < 0], base = 10) #x (-1)
a[mytau_star > 0] <- -log(x = mytau_star_p_fdr[mytau_star > 0], base = 10)
a[mytau_star_p_fdr > 0.05] <- 0
a[a < 0] <- 0 #5,141 positive associations. 
a_EUR <- a[,1:length(traits_EUR)]
a_EAS <- a[,(length(traits_EUR)+1):ncol(a)]

dim(a_EUR) 
dim(a_EAS) 
ncol(a_EUR) == length(traits_EUR)
ncol(a_EAS) == length(traits_EAS)

commontraits_EAS <- c("Baso","BMI","DBP","Eosino","Hb","Height","Ht","Lym","MCHC","MCH","MCV","Mono","Neutro","Plt","RBC","SBP","WBC","Arrhythmia","Asthma","BrCa","CAD","Cataract","CHF","Glaucoma","IS","Osteoporosis","PrCa","RA","T2D")
commontraits_EUR <- c("Basophil","BMI","DBP","Eosinophil","Hb","Height","Ht","Lymphocyte","MCHC","MCH","MCV","Monocyte","Neutrophil","Platelet","UKB2_145K.blood_RED_COUNT","SBP","WBC","AF","Asthma.ukbb","BrCa","CAD","Cataract","CHF","Glaucoma","IS","Osteoporosis","PrCa","PASS_Rheumatoid_Arthritis","T2D")
commontraits <- cbind(commontraits_EAS,commontraits_EUR)
commontraits <- commontraits[order(commontraits[,1],decreasing = F),] #alphabetical order
commontraits_EAS <- commontraits[,1]
commontraits_EUR <- commontraits[,2]
m_eas <- match(commontraits_EAS, traits_EAS)
m_eur <- match(commontraits_EUR, traits_EUR)

mytaustar_IMPACT <- cbind(mytau_star_EUR_IMPACT[,m_eur],mytau_star_EAS_IMPACT[,m_eas])#[-removethese,]
mytaustarse_IMPACT <- cbind(mytau_star_se_EUR_IMPACT[,m_eur],mytau_star_se_EAS_IMPACT[,m_eas])#[-removethese,]
a_EUR_m <- a_EUR[,m_eur]
a_EAS_m <- a_EAS[,m_eas]
mytau_star_logp_na <- cbind(a_EUR_m,a_EAS_m)

pdf("EURvsEAS_IndAnnots_PerTrait.pdf", height = 12, width = 12*1.75)
par(mfrow = c(5,7))
commontraits_nicenames <- commontraits_EAS 
eas_trait <- c(eas_trait, "T2D")
eas_trait <- sort(eas_trait, decreasing = F)
how_many_eur_ind_annots <- c()
how_many_eureas_ind_annots <- c()
for (i in 1:length(commontraits_EAS)){
  Correlation_Clumping <- get(paste0("Correlation_Clumping_",commontraits_EAS[i]))
  set_ind_annots <- Correlation_Clumping$Lead #independent annotations, leads here are based on RA 
  
  plot(mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)],mytaustar_IMPACT[set_ind_annots,i], pch = 19, col = "white", ylim = c(min(mytaustar_IMPACT[,i] - 1.96*mytaustarse_IMPACT[,i]),max(mytaustar_IMPACT[,i] + 1.96*mytaustarse_IMPACT[,i])), xlim = c(min(mytaustar_IMPACT[,i+length(commontraits_EAS)] - 1.96*mytaustarse_IMPACT[,i+length(commontraits_EAS)]),max(mytaustar_IMPACT[,i+length(commontraits_EAS)] + 1.96*mytaustarse_IMPACT[,i+length(commontraits_EAS)])), ylab = "EUR tau*",xlab ="EAS tau*", main = commontraits_nicenames[i], cex.lab = 1.2, cex.axis = 1.2)
  segments(x0 = mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)], x1 = mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)], y0 = mytaustar_IMPACT[set_ind_annots,i] - 1.96*mytaustarse_IMPACT[set_ind_annots,i], y1 = mytaustar_IMPACT[set_ind_annots,i] + 1.96*mytaustarse_IMPACT[set_ind_annots,i], col = "snow2")
  segments(x0 = mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)] - 1.96*mytaustarse_IMPACT[set_ind_annots,i+length(commontraits_EAS)], x1 = mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)]+ 1.96*mytaustarse_IMPACT[set_ind_annots,i+length(commontraits_EAS)], y0 = mytaustar_IMPACT[set_ind_annots,i], y1 = mytaustar_IMPACT[set_ind_annots,i], col = "snow2")
  weur <- which(mytau_star_logp_na[set_ind_annots,i] > 0)
  weas <- which(mytau_star_logp_na[set_ind_annots,i+length(commontraits_EAS)] > 0)
  mycols <- rep("snow3",nrow(mytau_star_logp_na[set_ind_annots,]))
  mycols[weas] <- rgb(red = 0.6, green = 0.2, blue = 0.2, alpha = 1)
  mycols[weur] <- rgb(red = 0.2, green = 0.6, blue = 0.2, alpha = 1)
  mycols[intersect(weas,weur)] <- rgb(red = 0.2, green = 0.2, blue = 0.6, alpha = 1)
  how_many_eur_ind_annots <- c(how_many_eur_ind_annots, length(weur))
  how_many_eureas_ind_annots <- c(how_many_eureas_ind_annots, length(intersect(weas,weur)))
  mybgcols <- rep("snow3",nrow(mytau_star_logp_na[set_ind_annots,]))
  mybgcols[weas] <- rgb(red = 0.6, green = 0.2, blue = 0.2, alpha = 0.5)
  mybgcols[weur] <- rgb(red = 0.2, green = 0.6, blue = 0.2, alpha = 0.5)
  mybgcols[intersect(weas,weur)] <- rgb(red = 0.2, green = 0.2, blue = 0.6, alpha = 0.5)
  
  plotfirst <- which(mycols == "snow3")
  plotsecond <- which(mycols != "snow3")
  
  points(mytaustar_IMPACT[set_ind_annots[plotfirst],i+length(commontraits_EAS)],mytaustar_IMPACT[set_ind_annots[plotfirst],i], bg = "snow2",pch = 21, col = mycols[plotfirst]) #col specifies border color
  segments(x0 = -100,y0 = -100, x1 = 100, y1 = 100, lty = 4, col = "black")
  segments(x0 = -100,y0 = 0, x1 = 100, y1 = 0, lty = 4, col = "gray")
  segments(x0 = 0,y0 = -100, x1 = 0, y1 = 100, lty = 4, col = "gray")
  segments(x0 = mytaustar_IMPACT[set_ind_annots[plotsecond],i+length(commontraits_EAS)], x1 = mytaustar_IMPACT[set_ind_annots[plotsecond],i+length(commontraits_EAS)], y0 = mytaustar_IMPACT[set_ind_annots[plotsecond],i] - 1.96*mytaustarse_IMPACT[set_ind_annots[plotsecond],i], y1 = mytaustar_IMPACT[set_ind_annots[plotsecond],i] + 1.96*mytaustarse_IMPACT[set_ind_annots[plotsecond],i], col = mycols[plotsecond])
  segments(x0 = mytaustar_IMPACT[set_ind_annots[plotsecond],i+length(commontraits_EAS)] - 1.96*mytaustarse_IMPACT[set_ind_annots[plotsecond],i+length(commontraits_EAS)], x1 = mytaustar_IMPACT[set_ind_annots[plotsecond],i+length(commontraits_EAS)]+ 1.96*mytaustarse_IMPACT[set_ind_annots[plotsecond],i+length(commontraits_EAS)], y0 = mytaustar_IMPACT[set_ind_annots[plotsecond],i], y1 = mytaustar_IMPACT[set_ind_annots[plotsecond],i], col = mycols[plotsecond])
  
  points(mytaustar_IMPACT[set_ind_annots[plotsecond],i+length(commontraits_EAS)],mytaustar_IMPACT[set_ind_annots[plotsecond],i], pch = 21, bg = "gray", col = mycols[plotsecond]) #mybgcols[plotsecond]
  
  ct <- cor.test(mytaustar_IMPACT[set_ind_annots,i],mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)])
  ct$p.val
  ct$estimate
  
  cor(mytaustar_IMPACT[set_ind_annots[which(mycols == "#333399FF")],i],mytaustar_IMPACT[set_ind_annots[which(mycols == "#333399FF")],i+length(commontraits_EAS)])
  
  fit <- deming::deming(formula = mytaustar_IMPACT[set_ind_annots,i] ~ mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)], xstd = mytaustarse_IMPACT[set_ind_annots,i+length(commontraits_EAS)], ystd = mytaustarse_IMPACT[set_ind_annots,i])
  beta <- fit$coefficients[2]
  beta_se <- (beta - fit$ci[2,1])/1.96
}
dev.off()

#plot all together
p <- c()
betas <- c()
eas_points <- c()
eur_points <- c()
cols_points <- c()
eur_sds <- c()
eas_sds <- c()
eur_points_traits <- c()
set_ind_annots_list <- c()
for (i in 1:length(commontraits_EAS)){
  print(i)
  Correlation_Clumping <- get(paste0("Correlation_Clumping_",commontraits_EAS[i]))
  set_ind_annots <- Correlation_Clumping$Lead #independent annotations, leads here are based on RA 
  
  eur_points <- c(eur_points, mytaustar_IMPACT[set_ind_annots,i])
  eas_points <- c(eas_points, mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)])
  eur_points_traits <- c(eur_points_traits, rep(commontraits_EAS[i],length(set_ind_annots)))
  eur_sds <- c(eur_sds, mytaustarse_IMPACT[set_ind_annots,i])
  eas_sds <- c(eas_sds, mytaustarse_IMPACT[set_ind_annots,i+length(commontraits_EAS)])
  set_ind_annots_list <- c(set_ind_annots_list, set_ind_annots)
  
  weur <- which(mytau_star_logp_na[set_ind_annots,i] > 0)
  weas <- which(mytau_star_logp_na[set_ind_annots,i+length(commontraits_EAS)] > 0)
  mycols <- rep("snow3",nrow(mytau_star_logp_na[set_ind_annots,]))
  mycols[weas] <- rgb(red = 0.6, green = 0.2, blue = 0.2, alpha = 1)
  mycols[weur] <- rgb(red = 0.2, green = 0.6, blue = 0.2, alpha = 1)
  mycols[intersect(weas,weur)] <- rgb(red = 0.2, green = 0.2, blue = 0.6, alpha = 1)
  cols_points <- c(cols_points,mycols)
  
  fit <- deming::deming(formula = mytaustar_IMPACT[set_ind_annots,i] ~ mytaustar_IMPACT[set_ind_annots,i+length(commontraits_EAS)], xstd = mytaustarse_IMPACT[set_ind_annots,i+length(commontraits_EAS)], ystd = mytaustarse_IMPACT[set_ind_annots,i])
  beta <- fit$coefficients[2]
  beta_se <- (beta - fit$ci[2,1])/1.96
  if(beta < 1){
    p[i] <- pnorm(q = 1, mean = beta, sd = beta_se, lower.tail = F)
  }else{
    p[i] <- pnorm(q = 1, mean = beta, sd = beta_se, lower.tail = T)
  }
  betas[i] <- beta
  
}

pdf("EURvsEAS_taustar_Concordance_allTraits.pdf", height = 7, width = 7)
par(mar = c(5.1, 5, 4.1, 2.1))
plot(0,0, xlim = c(-10,10), ylim = c(-10,10), col = "white", main = expression(tau*"* concordance across 29 EUR/EAS traits"), ylab = expression("EUR "*tau*"*"), xlab = expression("EAS "*tau*"*"), cex.axis = 1.3, cex.lab = 1.3, cex.main = 1.3)
segments(x0 = -100,y0 = 0, x1 = 100, y1 = 0, lty = 4, col = "gray")
segments(x0 = 0,y0 = -100, x1 = 0, y1 = 100, lty = 4, col = "gray")
w <- which(cols_points != "snow3")
points(eas_points[-w],eur_points[-w], col = alpha(cols_points,0.6)[-w], pch = 19)
points(eas_points[w],eur_points[w], col = alpha(cols_points,0.6)[w], pch = 19)
segments(x0 = -100,y0 = -100, x1 = 100, y1 = 100, lty = 4, col = "black")
legend("topleft", legend = c("Sig in Both","Sig in EUR","Sig in EAS","Sig in Neither"), title = expression("FDR 5% adjusted "*tau*"* p"), col = c("#333399FF", "#339933FF", "#993333FF","snow3"), bty = "n", pch = 19, adj = 0)
w <- which(cols_points == "snow3")
fit <- deming::deming(formula = eur_points[w] ~ eas_points[w], xstd = eas_sds[w], ystd = eur_sds[w])
w <- which(cols_points != "snow3")
fit_anysig <- deming::deming(formula = eur_points[w] ~ eas_points[w], xstd = eas_sds[w], ystd = eur_sds[w])
w <- which(cols_points == "#333399FF")
fit_bothsig <- deming::deming(formula = eur_points[w] ~ eas_points[w], xstd = eas_sds[w], ystd = eur_sds[w])
beta_sig <- fit_bothsig$coefficients[2]
int_sig <- fit_bothsig$coefficients[1]
beta_sig_se <- (beta_sig - fit_bothsig$ci[2,1])/1.96
beta_nonsig <- fit$coefficients[2]
int_nonsig <- fit$coefficients[1]
beta_nonsig_se <- (beta_nonsig - fit$ci[2,1])/1.96
pnorm(q = 1, mean = beta_sig, sd = beta_sig_se, lower.tail = F) #0.1826005 (any sig), 0.2773807 (both sig)
pnorm(q = 1, mean = beta_nonsig, sd = beta_nonsig_se, lower.tail = F) #1.493082e-09
abline(a =  int_sig, b = beta_sig, lwd = 2, col = "purple")
ci_x <- -10:10
ci_y_upper <- (beta_sig+1.96*beta_sig_se)*ci_x+int_sig
ci_y_lower <- (beta_sig-1.96*beta_sig_se)*rev(ci_x)+int_sig
polygon(c(ci_x,rev(ci_x)),c(ci_y_upper,ci_y_lower),border = NA, col=alpha("purple",0.4))

#annotations that are sig, higher slope than those that are not? 
mu_new <- beta_sig - beta_nonsig
sd_new <- sqrt(beta_sig_se^2 + beta_nonsig_se^2) 
pnorm(q = 0, mean = mu_new, sd = sd_new, lower.tail = T) #2.12e-5, significantly different from one another
legend("bottomright", legend = c("Identity","Slope (Sig in Both)"), col = c("black", "purple"), bty = "n", lty = c(4,1), lwd = 2, adj = 0)
text(-8.7,-9.3, col = "purple", expression(beta*" = "*"0.98 (sd = 0.04)"), adj = 0)
w <- which(cols_points != "snow3")
ct_nonsig <- cor.test(eur_points[-w],eas_points[-w]) #0.41
ct_sig <- cor.test(eur_points[w],eas_points[w]) #0.74

ct_sig_est <- ct_sig$estimate
ct_sig_sd <- (ct_sig$conf.int[2] - ct_sig_est) / 1.96

ct_nonsig_est <- ct_nonsig$estimate
ct_nonsig_sd <- (ct_nonsig$conf.int[2] - ct_nonsig_est) / 1.96

mu_new <- ct_sig_est - ct_nonsig_est
sd_new <- sqrt(ct_sig_sd^2 + ct_nonsig_sd^2) / sqrt(2)
pnorm(q = 0, mean = mu_new, sd = sd_new) #5.5e-62
#conclude: significant annotations are significantly more correlated in tau* than nonsignificant annotations. 
#conclude: not only are they more significantly correlated, but the regression slope is not significantly different from 1, suggesting signficant annotations capture the same amount of h2 in both populations. 

dev.off()

mu_new <- eur_points - eas_points
sd_new <- sqrt(eur_sds^2 + eas_sds^2)
hettest <- sapply(1:length(mu_new), function(x) ifelse(mu_new[x] > 0, pnorm(q = 0, mean = mu_new[x], sd = sd_new[x], lower.tail = T), pnorm(q = 0, mean = mu_new[x], sd = sd_new[x], lower.tail = F)))
length(which(hettest < 0.05))
hettest_fdr <- p.adjust(hettest, method = "fdr", n = length(hettest))
length(which(hettest_fdr > 0.05)) / length(hettest_fdr) #no evidence of heterogeneity 
