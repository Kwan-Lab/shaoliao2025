## dendrite plasticity: 
## lingxiao's 2p calcium imaging of dendrite branches and spines, acute pre/post exposure, it and pt neurons
# 
# notes:
# 
# 
# 20230107 $nks kwanlab @ cornell bme

## remove workspace variables
rm(list=ls())

## load some packages
library(car)
library(afex)
library(ez)
library(reshape2)
library(lme4)
library(multcomp)
library(emmeans)
library(lmtest)
library(lmerTest)
library(ggplot2)

## directories
proj_dir <- "/Users/neilsavalia/Documents/yale/kwanlab/project/dendrites_plasticity"
data_dir <- paste(proj_dir, "data", "psilocybin_pt_it_2p_calcium", sep = "/")
script_dir <- paste(proj_dir, "scripts", sep = "/")
output_dir <- paste(proj_dir, "output", "psilocybin_pt_it_2p_calcium", sep = "/")
figure_dir <- paste(proj_dir, "figure", sep = "/")

## load data table(s)
branch    <- read.csv(paste(data_dir, "dendrites_plasticity_pt_it_data_branches_20230107.csv", sep = "/"), header = TRUE)
spine     <- read.csv(paste(data_dir, "dendrites_plasticity_pt_it_data_spines_20230107.csv", sep = "/"), header = TRUE)
spine_reg <- read.csv(paste(data_dir, "dendrites_plasticity_pt_it_data_spines_reg_20230107.csv", sep = "/"), header = TRUE)

## calculate delta values ( [post - pre] / pre )
branch$eventrate_delta <- (branch$eventrate_post - branch$eventrate_pre) / branch$eventrate_pre
branch$amplitude_delta <- (branch$amplitude_post - branch$amplitude_pre) / branch$amplitude_pre
branch$frequency_delta <- (branch$frequency_post - branch$frequency_pre) / branch$frequency_pre
spine$eventrate_delta <- (spine$eventrate_post - spine$eventrate_pre) / spine$eventrate_pre
spine$amplitude_delta <- (spine$amplitude_post - spine$amplitude_pre) / spine$amplitude_pre
spine$frequency_delta <- (spine$frequency_post - spine$frequency_pre) / spine$frequency_pre
spine_reg$eventrate_delta <- (spine_reg$eventrate_post - spine_reg$eventrate_pre) / spine_reg$eventrate_pre
spine_reg$amplitude_delta <- (spine_reg$amplitude_post - spine_reg$amplitude_pre) / spine_reg$amplitude_pre
spine_reg$frequency_delta <- (spine_reg$frequency_post - spine_reg$frequency_pre) / spine_reg$frequency_pre

## LME: branch
for (analyze_branch in 1) 
{
  ## data setup
  dat <- branch
  dat.name <- "eventrate"
  dat$dendrite_id <- factor(paste(dat$subject, dat$treatment, dat$fov, dat$branch, sep = "_"))
  dat$subject <- factor(dat$subject)
  dat$treatment <- factor(dat$treatment)
  dat$cell_type <- factor(dat$cell_type)
  dat$order <- factor(dat$order)
  dat$fov <- factor(dat$fov)
  dat[is.na(dat) | dat=="Inf"] = NA
  
  # remove eventrate outlier from branch data
  datfilt <- dat[dat[, "eventrate_delta"] < 8, ]
  
  ## LME: branch
  mod3 <- lmer(eventrate_delta ~ treatment * cell_type + order + (1 | subject / fov), data = datfilt)
  anova(mod3)
  # summary(mod3)
  # plot(mod3) # visual check on residual linearity & variance
  # qqnorm(residuals(mod3)) # visual check on normality of residuals
  # fixef(mod3)
  # ranef(mod3)
  # coef(mod3)
  
  ## write out anova table
  write.csv(x = anova(mod3), paste(output_dir, paste("psilocybin_pt_it_gcamp6f_2p_R_lme_model_output_branch_", dat.name, "_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))
  
  # ## dot plot branch data
  # dodge <- position_dodge(width = 0.75)
  # ggplot(dat, aes(x = factor(cell_type), y = 100*eventrate_delta, color = factor(treatment))) +
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  #   geom_violin(trim = FALSE, width = .75, position = dodge) + 
  #   geom_boxplot(outlier.shape = NA, width = 0.2, position = dodge) +
  #   # geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 1, alpha = 0.5) +
  #   scale_color_manual(values = c("Red", "Black")) +
  #   ylim(-75, 250) +
  #   theme_bw() +
  #   theme(legend.position="none") +
  #   theme(panel.border = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5, size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #   labs(title = "Change in Branch Event Rate by Cell Type", x = "Cell Type", y = expression(Delta * "Event Rate (%) [(Post - Pre) / Pre]"))
  
  # grouped sample sizes
  n.psi = c(n.all = sum(datfilt$treatment == "psilocybin"))
  n.psi.pt = c(n.all = sum(datfilt$treatment == "psilocybin" & datfilt$cell_type == "PT"))
  n.psi.it = c(n.all = sum(datfilt$treatment == "psilocybin" & datfilt$cell_type == "IT"))
  n.sal = c(n.all = sum(datfilt$treatment == "saline"))
  n.sal.pt = c(n.all = sum(datfilt$treatment == "saline" & datfilt$cell_type == "PT"))
  n.sal.it = c(n.all = sum(datfilt$treatment == "saline" & datfilt$cell_type == "IT"))
  
  # 2-sample t-test of means (psilocybin minus saline)
  meas <- datfilt$eventrate_delta
  all_x <- meas[datfilt$treatment == "psilocybin"]
  all_y <- meas[datfilt$treatment == "saline"]
  it_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "IT"]
  it_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "IT"]
  pt_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "PT"]
  pt_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "PT"]
  testmeans.all <- t.test(x = all_x, y = all_y, alternative = "two.sided")
  testmeans.it <- t.test(x = it_x, y = it_y, alternative = "two.sided")
  testmeans.pt <- t.test(x = pt_x, y = pt_y, alternative = "two.sided")
  all.tab <- c(testmeans.all$estimate[1], testmeans.all$estimate[2], diff = testmeans.all$estimate[1] - testmeans.all$estimate[2], stderr = testmeans.all$stderr, n1 = n.psi[1], n2 = n.sal[1], testmeans.all$statistic, testmeans.all$parameter, testmeans.all$p.value, testmeans.all$conf.int[1], testmeans.all$conf.int[2])
  names(all.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
  it.tab <- c(testmeans.it$estimate[1], testmeans.it$estimate[2], diff = testmeans.it$estimate[1] - testmeans.it$estimate[2], stderr = testmeans.it$stderr, n1 = n.psi.it[1], n2 = n.sal.it[1], testmeans.it$statistic, testmeans.it$parameter, testmeans.it$p.value, testmeans.it$conf.int[1], testmeans.it$conf.int[2])
  names(all.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
  pt.tab <- c(testmeans.pt$estimate[1], testmeans.pt$estimate[2], diff = testmeans.pt$estimate[1] - testmeans.pt$estimate[2], stderr = testmeans.pt$stderr, n1 = n.psi.pt[1], n2 = n.sal.pt[1], testmeans.pt$statistic, testmeans.pt$parameter, testmeans.pt$p.value, testmeans.pt$conf.int[1], testmeans.pt$conf.int[2])
  names(all.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
  
  # accumulate & write out post hoc 2-sample t-test tables
  conc.all.tab <- rbind(all.tab, it.tab, pt.tab)
  # print(conc.all.tab)
  write.csv(x = conc.all.tab, paste(output_dir, paste("psilocybin_pt_it_gcamp6f_2p_R_lme_post_hoc_ttests_branch_", dat.name, "_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))
  
  # # 2-sample wilcoxon ranksum tests of medians (psilocybin minus saline)
  # meas <- datfilt$eventrate_delta
  # all_x <- meas[datfilt$treatment == "psilocybin"]
  # all_y <- meas[datfilt$treatment == "saline"]
  # it_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "IT"]
  # it_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "IT"]
  # pt_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "PT"]
  # pt_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "PT"]
  # testmeds.all <- wilcox.test(x = all_x, y = all_y)
  # testmeds.it <- wilcox.test(x = it_x, y = it_y)
  # testmeds.pt <- wilcox.test(x = pt_x, y = pt_y)
  
  # # Various 2-sample variance tests (psilocybin vs. saline)
  # meas <- datfilt$eventrate_delta
  # all_x <- meas[datfilt$treatment == "psilocybin"]
  # all_y <- meas[datfilt$treatment == "saline"]
  # it_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "IT"]
  # it_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "IT"]
  # pt_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "PT"]
  # pt_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "PT"]
  # # var.test(x = it_x, y = it_y, alternative = "two.sided")
  # # var.test(x = pt_x, y = pt_y, alternative = "two.sided")
  # # bartlett.test(x = c(it_x, it_y), g = c(rep(1, length(it_x)), rep(2, length(it_y))), alternative = "two.sided")
  # # bartlett.test(x = c(pt_x, pt_y), g = c(rep(1, length(pt_x)), rep(2, length(pt_y))), alternative = "two.sided")
  # testvars.all <- leveneTest(y = c(all_x, all_y), group = factor(c(rep(1, length(all_x)), rep(2, length(all_y)))))
  # testvars.it <- leveneTest(y = c(it_x, it_y), group = factor(c(rep(1, length(it_x)), rep(2, length(it_y)))))
  # testvars.pt <- leveneTest(y = c(pt_x, pt_y), group = factor(c(rep(1, length(pt_x)), rep(2, length(pt_y)))))
  
}

## LME: spine (branch-regressed)
for (analyze_spine_reg in 1)
{
  
  ## data setup
  dat <- spine_reg
  dat.name <- "eventrate"
  dat$dendrite_id <- factor(paste(dat$subject, dat$treatment, dat$fov, dat$branch, sep = "_"))
  dat$subject <- factor(dat$subject)
  dat$treatment <- factor(dat$treatment)
  dat$cell_type <- factor(dat$cell_type)
  dat$order <- factor(dat$order)
  dat$fov <- factor(dat$fov)
  dat[is.na(dat) | dat=="Inf"] = NA
  
  # remove eventrate outliers from spine_reg data
  datfilt <- dat[dat[, "eventrate_delta"] < 50, ]
  
  ## LME: spine (branch-regressed)
  mod3 <- lmer(eventrate_delta ~ treatment * cell_type + order + (1 | subject / fov / dendrite_id), data = datfilt)
  anova(mod3)
  # summary(mod3)
  # plot(mod3) # visual check on residual linearity & variance
  # qqnorm(residuals(mod3)) # visual check on normality of residuals
  # fixef(mod3)
  # ranef(mod3)
  # coef(mod3)
  
  ## write out anova table
  write.csv(x = anova(mod3), paste(output_dir, paste("psilocybin_pt_it_gcamp6f_2p_R_lme_model_output_spinereg_", dat.name, "_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))
  
  # # ## dot plot spine_reg data
  # dodge <- position_dodge(width = 0.75)
  # ggplot(dat, aes(x = factor(cell_type), y = 100*eventrate_delta, color = factor(treatment))) +
  #   geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  #   geom_violin(trim = FALSE, width = .75, position = dodge) +
  #   geom_boxplot(outlier.shape = NA, width = .2, position = dodge) +
  #   # geom_point(position = position_jitterdodge(jitter.width = 0.4), size = .1, alpha = 0.25) +
  #   scale_color_manual(values = c("Red", "Black")) +
  #   ylim(-100, 1000) +
  #   theme_bw() +
  #   theme(legend.position="none") +
  #   theme(panel.border = element_blank(), plot.title = element_text(hjust = 0.5, size = 16), axis.text = element_text(size = 16), axis.title = element_text(size = 16)) +
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  #   labs(title = "Change in Spine Event Rate by Cell Type", x = "Cell Type", y = expression(Delta * "Event Rate (%) [(Post - Pre) / Pre]"))
  
  # grouped sample sizes
  n.psi = c(n.all = sum(datfilt$treatment == "psilocybin", na.rm = TRUE))
  n.psi.pt = c(n.all = sum(datfilt$treatment == "psilocybin" & datfilt$cell_type == "PT", na.rm = TRUE))
  n.psi.it = c(n.all = sum(datfilt$treatment == "psilocybin" & datfilt$cell_type == "IT", na.rm = TRUE))
  n.sal = c(n.all = sum(datfilt$treatment == "saline", na.rm = TRUE))
  n.sal.pt = c(n.all = sum(datfilt$treatment == "saline" & datfilt$cell_type == "PT", na.rm = TRUE))
  n.sal.it = c(n.all = sum(datfilt$treatment == "saline" & datfilt$cell_type == "IT", na.rm = TRUE))
  
  # 2-sample t-test of means (psilocybin minus saline)
  meas <- datfilt$eventrate_delta
  all_x <- meas[datfilt$treatment == "psilocybin"]
  all_y <- meas[datfilt$treatment == "saline"]
  it_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "IT"]
  it_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "IT"]
  pt_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "PT"]
  pt_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "PT"]
  testmeans.all <- t.test(x = all_x, y = all_y, alternative = "two.sided")
  testmeans.it <- t.test(x = it_x, y = it_y, alternative = "two.sided")
  testmeans.pt <- t.test(x = pt_x, y = pt_y, alternative = "two.sided")
  all.tab <- c(testmeans.all$estimate[1], testmeans.all$estimate[2], diff = testmeans.all$estimate[1] - testmeans.all$estimate[2], stderr = testmeans.all$stderr, n1 = n.psi[1], n2 = n.sal[1], testmeans.all$statistic, testmeans.all$parameter, testmeans.all$p.value, testmeans.all$conf.int[1], testmeans.all$conf.int[2])
  names(all.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
  it.tab <- c(testmeans.it$estimate[1], testmeans.it$estimate[2], diff = testmeans.it$estimate[1] - testmeans.it$estimate[2], stderr = testmeans.it$stderr, n1 = n.psi.it[1], n2 = n.sal.it[1], testmeans.it$statistic, testmeans.it$parameter, testmeans.it$p.value, testmeans.it$conf.int[1], testmeans.it$conf.int[2])
  names(all.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
  pt.tab <- c(testmeans.pt$estimate[1], testmeans.pt$estimate[2], diff = testmeans.pt$estimate[1] - testmeans.pt$estimate[2], stderr = testmeans.pt$stderr, n1 = n.psi.pt[1], n2 = n.sal.pt[1], testmeans.pt$statistic, testmeans.pt$parameter, testmeans.pt$p.value, testmeans.pt$conf.int[1], testmeans.pt$conf.int[2])
  names(all.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
  
  # accumulate & write out post hoc 2-sample t-test tables
  conc.all.tab <- rbind(all.tab, it.tab, pt.tab)
  # print(conc.all.tab)
  write.csv(x = conc.all.tab, paste(output_dir, paste("psilocybin_pt_it_gcamp6f_2p_R_lme_post_hoc_ttests_spinereg_", dat.name, "_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))
  
  # # 2-sample wilcoxon ranksum tests of medians (psilocybin minus saline)
  # meas <- datfilt$eventrate_delta
  # all_x <- meas[datfilt$treatment == "psilocybin"]
  # all_y <- meas[datfilt$treatment == "saline"]
  # it_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "IT"]
  # it_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "IT"]
  # pt_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "PT"]
  # pt_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "PT"]
  # testmeds.all <- wilcox.test(x = all_x, y = all_y)
  # testmeds.it <- wilcox.test(x = it_x, y = it_y)
  # testmeds.pt <- wilcox.test(x = pt_x, y = pt_y)
  
  # # Various 2-sample variance tests (psilocybin vs. saline)
  # meas <- datfilt$eventrate_delta
  # all_x <- meas[datfilt$treatment == "psilocybin"]
  # all_y <- meas[datfilt$treatment == "saline"]
  # it_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "IT"]
  # it_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "IT"]
  # pt_x <- meas[datfilt$treatment == "psilocybin" & datfilt$cell_type == "PT"]
  # pt_y <- meas[datfilt$treatment == "saline" & datfilt$cell_type == "PT"]
  # # var.test(x = it_x, y = it_y, alternative = "two.sided")
  # # var.test(x = pt_x, y = pt_y, alternative = "two.sided")
  # # bartlett.test(x = c(it_x, it_y), g = c(rep(1, length(it_x)), rep(2, length(it_y))), alternative = "two.sided")
  # # bartlett.test(x = c(pt_x, pt_y), g = c(rep(1, length(pt_x)), rep(2, length(pt_y))), alternative = "two.sided")
  # testvars.all <- leveneTest(y = c(all_x, all_y), group = factor(c(rep(1, length(all_x)), rep(2, length(all_y)))))
  # testvars.it <- leveneTest(y = c(it_x, it_y), group = factor(c(rep(1, length(it_x)), rep(2, length(it_y)))))
  # testvars.pt <- leveneTest(y = c(pt_x, pt_y), group = factor(c(rep(1, length(pt_x)), rep(2, length(pt_y)))))
  
}





