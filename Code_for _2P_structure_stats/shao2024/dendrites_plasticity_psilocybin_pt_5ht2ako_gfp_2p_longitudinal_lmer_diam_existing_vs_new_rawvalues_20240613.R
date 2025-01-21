## dendrites plasticity: 
## lingxiao longitudinal 2P imaging data of cg1/m2 pt neuron structure, psilocybin vs. saline, pt neurons in 5ht2a f/f mice (5ht2ako's)
# 
# measures: spine head width for pre-existing vs newly-formed spines
# 
# notes:
# run linear mixed effects models of spine head width (raw value) contrasting pre-existing vs newly-formed spines
# 
# 20240610 $nks kwanlab @ cornell bme

## remove workspace variables
rm(list=ls())

## load some packages
library(reshape2)
library(lme4)
library(multcomp)
library(emmeans)
library(lmtest)
library(lmerTest)
library(ggplot2)

## directories
proj_dir <- "/Users/neilsavalia/Documents/yale/kwanlab/project/dendrites_plasticity"
data_dir <- paste(proj_dir, "data", "psilocybin_pt_it_gfp_2p_structure", "20240608_pt_c57_vs_5ht2ako", sep = "/")
script_dir <- paste(proj_dir, "scripts", sep = "/")
output_dir <- paste(proj_dir, "output", "psilocybin_pt_it_gfp_2p_structure", "pt_5ht2ako_rawvalues", sep = "/")
figure_dir <- paste(proj_dir, "figure", sep = "/")

## load data table(s)
diam  <- read.csv(paste(data_dir, "gfp_2p_psilocybin_pt_5ht2ako_spine_headwidth_rawvalues_existingvsnew_long_20240608.csv", sep = "/"), header = TRUE)

## select data and reformat
#   - table must be in long format (time columnized)
dat <- diam
dat.name <- "headwidth"
dat$mouse_id <- factor(sub(" ", "_", dat$mouse_id))
dat$genotype <- factor(dat$genotype)
dat$cell_id <- factor(paste(dat$mouse_id, dat$cell_id, sep = "_"))
dat$cell_id[grepl("NA", dat$cell_id)] = NA
dat$branch_id <- factor(paste(dat$mouse_id, dat$branch_id, sep = "_"))
dat$branch_id[grepl("NA", dat$branch_id)] = NA
dat$multibranch_id <- factor(paste(dat$mouse_id, dat$multibranch_id, sep = "_"))
dat$multibranch_id[grepl("NA", dat$multibranch_id)] = NA
dat$sex <- factor(dat$sex)
dat$new <- factor(dat$new)
dat$treatment <- factor(dat$treatment)
dat$soma_depth[grepl("\\+", dat$soma_depth)] = NA
dat$soma_depth <- as.numeric(dat$soma_depth)
dat$fov <- factor(dat$fov)
datLong <- melt(data = dat,
                id.vars = c("mouse_id", "sex", "treatment", "genotype", "fov", "soma_depth", "cell_id", "branch_id", "multibranch_id", "new"),
                measure.vars = c("day_1", "day_3"),
                variable.name = "time",
                value.name = "value",
                variable.factor = TRUE)
# tapply(datLong$value, datLong$sex, mean, na.rm = TRUE)
# tapply(datLong$value, datLong$treatment, mean, na.rm = TRUE)
# tapply(datLong$value, datLong$time, mean, na.rm = TRUE)
# tapply(datLong$value, datLong$genotype, mean, na.rm = TRUE)
datLong_ko <- subset(datLong, genotype == "ko")
datLong_c57 <- subset(datLong, genotype == "c57")
datLong.tcont <- datLong
levels(datLong.tcont$time)[levels(datLong.tcont$time) == "day_1"] <- as.numeric(1)
levels(datLong.tcont$time)[levels(datLong.tcont$time) == "day_3"] <- as.numeric(3)
datLong.tcont$time <- as.numeric(paste(datLong.tcont$time))

## LME - value ~ treatment * time * genotype + controlling for nested random effects of mouse_id / cell_id / multibranch_id
# mod3 <- lmer(value ~ treatment * time * genotype + (1 | mouse_id / cell_id / multibranch_id), data = datLong.tcont)  # multibranch_id drops info about known cells that we only measure a single branch, but these should be treated separately from NAs with no cell provenance info.! 06/11/24
mod3 <- lmer(value ~ treatment * time * genotype * new + (1 | mouse_id / cell_id / branch_id), data = datLong.tcont)
anova(mod3)
# summary(mod3)
# plot(mod3) # visual check on residual linearity & variance
# qqnorm(residuals(mod3)) # visual check on normality of residuals
# fixef(mod3)
# ranef(mod3)
# coef(mod3)
write.csv(x = anova(mod3), paste(output_dir, paste("gfp_2P_psilocybin_pt_5ht2ako_R_lme_model_output_spine_headwidth_existing_vs_new_foldchange_bigmodel_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))


# ## FULL post hoc t-tests comparison table, shorthand
# mm <- emmeans(mod3, specs = c("treatment", "time", "genotype", "sex"))
# cc <- summary(contrast(mm, method = "pairwise", adjust = "bonferroni"))


## LME - value ~ treatment * time * genotype + controlling for nested random effects of mouse_id / cell_id / multibranch_id
# mod3 <- lmer(value ~ treatment * time * genotype + (1 | mouse_id / cell_id / multibranch_id), data = datLong.tcont)  # multibranch_id drops info about known cells that we only measure a single branch, but these should be treated separately from NAs with no cell provenance info.! 06/11/24
mod3 <- lmer(value ~ treatment * genotype * new + (1 | mouse_id / cell_id / branch_id), data = datLong.tcont)
anova(mod3)
# summary(mod3)
# plot(mod3) # visual check on residual linearity & variance
# qqnorm(residuals(mod3)) # visual check on normality of residuals
# fixef(mod3)
# ranef(mod3)
# coef(mod3)
write.csv(x = anova(mod3), paste(output_dir, paste("gfp_2P_psilocybin_pt_5ht2ako_R_lme_model_output_spine_headwidth_existing_vs_new_foldchange_smallmodel_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ################################# Selected post hoc t-tests ################################# #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Get sample sizes for various groupings and sub-groupings
#####
n.psi = c(n.all = sum(dat$treatment == "psilocybin"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin"])))
n.psi.old = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 0), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 0])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 0])))
n.psi.new = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 1), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 1])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 1])))
n.psi.old.c57 = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "c57"])))
n.psi.new.c57 = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "c57"])))
n.psi.old.ko = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "ko"])))
n.psi.new.ko = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "ko"])))
n.sal = c(n.all = sum(dat$treatment == "saline"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline"])))
n.sal.old = c(n.all = sum(dat$treatment == "saline" & dat$new == 0), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 0])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 0])))
n.sal.new = c(n.all = sum(dat$treatment == "saline" & dat$new == 1), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 1])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 1])))
n.sal.old.c57 = c(n.all = sum(dat$treatment == "saline" & dat$new == 0 & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 0 & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 0 & dat$genotype == "c57"])))
n.sal.new.c57 = c(n.all = sum(dat$treatment == "saline" & dat$new == 1 & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 1 & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 1 & dat$genotype == "c57"])))
n.sal.old.ko = c(n.all = sum(dat$treatment == "saline" & dat$new == 0 & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 0 & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 0 & dat$genotype == "ko"])))
n.sal.new.ko = c(n.all = sum(dat$treatment == "saline" & dat$new == 1 & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 1 & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 1 & dat$genotype == "ko"])))




n.old = c(n.all = sum(dat$new == 0), n1 = sum(!is.na(dat$day_1[dat$new == 0])), n3 = sum(!is.na(dat$day_3[dat$new == 0])))
n.old.psi = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 0), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 0])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 0])))
n.old.sal = c(n.all = sum(dat$treatment == "saline" & dat$new == 0), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 1])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 0])))
n.old.psi.c57 = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "c57"])))
n.old.sal.c57 = c(n.all = sum(dat$treatment == "saline" & dat$new == 0 & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 0 & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 0 & dat$genotype == "c57"])))
n.old.psi.ko = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 0 & dat$genotype == "ko"])))
n.old.sal.ko = c(n.all = sum(dat$treatment == "saline" & dat$new == 0 & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 0 & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 0 & dat$genotype == "ko"])))
n.new = c(n.all = sum(dat$new == 1), n1 = sum(!is.na(dat$day_1[dat$new == 1])), n3 = sum(!is.na(dat$day_3[dat$new == 1])))
n.new.psi = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 1), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 1])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 1])))
n.new.sal = c(n.all = sum(dat$treatment == "saline" & dat$new == 1), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 1])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 1])))
n.new.psi.c57 = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "c57"])))
n.new.sal.c57 = c(n.all = sum(dat$treatment == "saline" & dat$new == 1 & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 1 & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 1 & dat$genotype == "c57"])))
n.new.psi.ko = c(n.all = sum(dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$new == 1 & dat$genotype == "ko"])))
n.new.sal.ko = c(n.all = sum(dat$treatment == "saline" & dat$new == 1 & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$new == 1 & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$new == 1 & dat$genotype == "ko"])))



#####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# TWO-SAMPLE T-TESTS: treatment effect over time, altogether and split up by genotype and pre-existing vs new
#####

## TWO-SAMPLE T-TESTS: treatment effect over time, agnostic to genotype, pre-existing spines
condt.old.1 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$new == 0], y = dat$day_1[dat$treatment == "saline" & dat$new == 0], alternative = "two.sided")
condt.old.3 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$new == 0], y = dat$day_3[dat$treatment == "saline" & dat$new == 0], alternative = "two.sided")
condt.old.1.tab <- c(condt.old.1$estimate[1], condt.old.1$estimate[2], diff = condt.old.1$estimate[1] - condt.old.1$estimate[2], stderr = condt.old.1$stderr, n1 = n.psi.old[2], n2 = n.sal.old[2], condt.old.1$statistic, condt.old.1$parameter, condt.old.1$p.value, condt.old.1$conf.int[1], condt.old.1$conf.int[2])
names(condt.old.1.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.old.3.tab <- c(condt.old.3$estimate[1], condt.old.3$estimate[2], diff = condt.old.3$estimate[1] - condt.old.3$estimate[2], stderr = condt.old.3$stderr, n1 = n.psi.old[3], n2 = n.sal.old[3], condt.old.3$statistic, condt.old.3$parameter, condt.old.3$p.value, condt.old.3$conf.int[1], condt.old.3$conf.int[2])
names(condt.old.3.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effect over time, agnostic to genotype, newly-formed spines
condt.new.1 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$new == 1], y = dat$day_1[dat$treatment == "saline" & dat$new == 1], alternative = "two.sided")
condt.new.3 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$new == 1], y = dat$day_3[dat$treatment == "saline" & dat$new == 1], alternative = "two.sided")
condt.new.1.tab <- c(condt.new.1$estimate[1], condt.new.1$estimate[2], diff = condt.new.1$estimate[1] - condt.new.1$estimate[2], stderr = condt.new.1$stderr, n1 = n.psi.new[2], n2 = n.sal.new[2], condt.new.1$statistic, condt.new.1$parameter, condt.new.1$p.value, condt.new.1$conf.int[1], condt.new.1$conf.int[2])
names(condt.new.1.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.new.3.tab <- c(condt.new.3$estimate[1], condt.new.3$estimate[2], diff = condt.new.3$estimate[1] - condt.new.3$estimate[2], stderr = condt.new.3$stderr, n1 = n.psi.new[3], n2 = n.sal.new[3], condt.new.3$statistic, condt.new.3$parameter, condt.new.3$p.value, condt.new.3$conf.int[1], condt.new.3$conf.int[2])
names(condt.new.3.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effects over time for c57 pre-existing spines
condt.old.1.c57 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "c57" & dat$new == 0], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "c57" & dat$new == 0], alternative = "two.sided")
condt.old.3.c57 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "c57" & dat$new == 0], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "c57" & dat$new == 0], alternative = "two.sided")
condt.old.1.c57.tab <- c(condt.old.1.c57$estimate[1], condt.old.1.c57$estimate[2], diff = condt.old.1.c57$estimate[1] - condt.old.1.c57$estimate[2], stderr = condt.old.1.c57$stderr, n1 = (n.psi.old.c57[2]), n2 = (n.sal.old.c57[2]), condt.old.1.c57$statistic, condt.old.1.c57$parameter, condt.old.1.c57$p.value, condt.old.1.c57$conf.int[1], condt.old.1.c57$conf.int[2])
names(condt.old.1.c57.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.old.3.c57.tab <- c(condt.old.3.c57$estimate[1], condt.old.3.c57$estimate[2], diff = condt.old.3.c57$estimate[1] - condt.old.3.c57$estimate[2], stderr = condt.old.3.c57$stderr, n1 = (n.psi.old.c57[3]), n2 = (n.sal.old.c57[3]), condt.old.3.c57$statistic, condt.old.3.c57$parameter, condt.old.3.c57$p.value, condt.old.3.c57$conf.int[1], condt.old.3.c57$conf.int[2])
names(condt.old.3.c57.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effects over time for c57 newly-formed spines
condt.new.1.c57 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "c57" & dat$new == 1], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "c57" & dat$new == 1], alternative = "two.sided")
condt.new.3.c57 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "c57" & dat$new == 1], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "c57" & dat$new == 1], alternative = "two.sided")
condt.new.1.c57.tab <- c(condt.new.1.c57$estimate[1], condt.new.1.c57$estimate[2], diff = condt.new.1.c57$estimate[1] - condt.new.1.c57$estimate[2], stderr = condt.new.1.c57$stderr, n1 = (n.psi.new.c57[2]), n2 = (n.sal.new.c57[2]), condt.new.1.c57$statistic, condt.new.1.c57$parameter, condt.new.1.c57$p.value, condt.new.1.c57$conf.int[1], condt.new.1.c57$conf.int[2])
names(condt.new.1.c57.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.new.3.c57.tab <- c(condt.new.3.c57$estimate[1], condt.new.3.c57$estimate[2], diff = condt.new.3.c57$estimate[1] - condt.new.3.c57$estimate[2], stderr = condt.new.3.c57$stderr, n1 = (n.psi.new.c57[3]), n2 = (n.sal.new.c57[3]), condt.new.3.c57$statistic, condt.new.3.c57$parameter, condt.new.3.c57$p.value, condt.new.3.c57$conf.int[1], condt.new.3.c57$conf.int[2])
names(condt.new.3.c57.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effects over time for KO pre-existing spines
condt.old.1.ko <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "ko" & dat$new == 0], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "ko" & dat$new == 0], alternative = "two.sided")
condt.old.3.ko <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "ko" & dat$new == 0], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "ko" & dat$new == 0], alternative = "two.sided")
condt.old.1.ko.tab <- c(condt.old.1.ko$estimate[1], condt.old.1.ko$estimate[2], diff = condt.old.1.ko$estimate[1] - condt.old.1.ko$estimate[2], stderr = condt.old.1.ko$stderr, n1 = (n.psi.old.ko[2]), n2 = (n.sal.old.ko[2]), condt.old.1.ko$statistic, condt.old.1.ko$parameter, condt.old.1.ko$p.value, condt.old.1.ko$conf.int[1], condt.old.1.ko$conf.int[2])
names(condt.old.1.ko.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.old.3.ko.tab <- c(condt.old.3.ko$estimate[1], condt.old.3.ko$estimate[2], diff = condt.old.3.ko$estimate[1] - condt.old.3.ko$estimate[2], stderr = condt.old.3.ko$stderr, n1 = (n.psi.old.ko[3]), n2 = (n.sal.old.ko[3]), condt.old.3.ko$statistic, condt.old.3.ko$parameter, condt.old.3.ko$p.value, condt.old.3.ko$conf.int[1], condt.old.3.ko$conf.int[2])
names(condt.old.3.ko.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effects over time for KO newly-formed spines
condt.new.1.ko <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "ko" & dat$new == 1], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "ko" & dat$new == 1], alternative = "two.sided")
condt.new.3.ko <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "ko" & dat$new == 1], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "ko" & dat$new == 1], alternative = "two.sided")
condt.new.1.ko.tab <- c(condt.new.1.ko$estimate[1], condt.new.1.ko$estimate[2], diff = condt.new.1.ko$estimate[1] - condt.new.1.ko$estimate[2], stderr = condt.new.1.ko$stderr, n1 = (n.psi.new.ko[2]), n2 = (n.sal.new.ko[2]), condt.new.1.ko$statistic, condt.new.1.ko$parameter, condt.new.1.ko$p.value, condt.new.1.ko$conf.int[1], condt.new.1.ko$conf.int[2])
names(condt.new.1.ko.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.new.3.ko.tab <- c(condt.new.3.ko$estimate[1], condt.new.3.ko$estimate[2], diff = condt.new.3.ko$estimate[1] - condt.new.3.ko$estimate[2], stderr = condt.new.3.ko$stderr, n1 = (n.psi.new.ko[3]), n2 = (n.sal.new.ko[3]), condt.new.3.ko$statistic, condt.new.3.ko$parameter, condt.new.3.ko$p.value, condt.new.3.ko$conf.int[1], condt.new.3.ko$conf.int[2])
names(condt.new.3.ko.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")


# ACCUMULATE
cond.all.tab <- rbind(
  condt.old.1.tab, condt.old.3.tab, 
  condt.new.1.tab, condt.new.3.tab, 
  condt.old.1.c57.tab, condt.old.3.c57.tab, 
  condt.new.1.c57.tab, condt.new.3.c57.tab, 
  condt.old.1.ko.tab, condt.old.3.ko.tab, 
  condt.new.1.ko.tab, condt.new.3.ko.tab)
# print(cond.all.tab)

#####

write.csv(x = cond.all.tab, paste(output_dir, paste("gfp_2P_psilocybin_pt_5ht2ako_R_lme_post_hoc_ttests_spine_headwidth_treatment_effect_rawvalues_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# TWO-SAMPLE T-TESTS: pre-existing vs newly-formed effect over time, split up by treatment and genotype
#####

## TWO-SAMPLE T-TESTS: pre-existing vs. newly-formed spines, compare within-treatment (psilocybin) over time, both genotypes
condt.psi.1 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$new == 0], y = dat$day_1[dat$treatment == "psilocybin" & dat$new == 1], alternative = "two.sided")
condt.psi.3 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$new == 0], y = dat$day_3[dat$treatment == "psilocybin" & dat$new == 1], alternative = "two.sided")
condt.psi.1.tab <- c(condt.psi.1$estimate[1], condt.psi.1$estimate[2], diff = condt.psi.1$estimate[1] - condt.psi.1$estimate[2], stderr = condt.psi.1$stderr, n1 = n.old.psi[2], n2 = n.new.psi[2], condt.psi.1$statistic, condt.psi.1$parameter, condt.psi.1$p.value, condt.psi.1$conf.int[1], condt.psi.1$conf.int[2])
names(condt.psi.1.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.psi.3.tab <- c(condt.psi.3$estimate[1], condt.psi.3$estimate[2], diff = condt.psi.3$estimate[1] - condt.psi.3$estimate[2], stderr = condt.psi.3$stderr, n1 = n.old.psi[3], n2 = n.new.psi[3], condt.psi.3$statistic, condt.psi.3$parameter, condt.psi.3$p.value, condt.psi.3$conf.int[1], condt.psi.3$conf.int[2])
names(condt.psi.3.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: pre-existing vs. newly-formed spines, compare within-treatment (saline) over time, both genotypes
condt.sal.1 <- t.test(x = dat$day_1[dat$treatment == "saline" & dat$new == 0], y = dat$day_1[dat$treatment == "saline" & dat$new == 1], alternative = "two.sided")
condt.sal.3 <- t.test(x = dat$day_3[dat$treatment == "saline" & dat$new == 0], y = dat$day_3[dat$treatment == "saline" & dat$new == 1], alternative = "two.sided")
condt.sal.1.tab <- c(condt.sal.1$estimate[1], condt.sal.1$estimate[2], diff = condt.sal.1$estimate[1] - condt.sal.1$estimate[2], stderr = condt.sal.1$stderr, n1 = n.old.sal[2], n2 = n.new.sal[2], condt.sal.1$statistic, condt.sal.1$parameter, condt.sal.1$p.value, condt.sal.1$conf.int[1], condt.sal.1$conf.int[2])
names(condt.sal.1.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.sal.3.tab <- c(condt.sal.3$estimate[1], condt.sal.3$estimate[2], diff = condt.sal.3$estimate[1] - condt.sal.3$estimate[2], stderr = condt.sal.3$stderr, n1 = n.old.sal[3], n2 = n.new.sal[3], condt.sal.3$statistic, condt.sal.3$parameter, condt.sal.3$p.value, condt.sal.3$conf.int[1], condt.sal.3$conf.int[2])
names(condt.sal.3.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: pre-existing vs. newly-formed spines, compare within-treatment (psilocybin) over time, wt
condt.psi.1.c57 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "c57" & dat$new == 0], y = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "c57" & dat$new == 1], alternative = "two.sided")
condt.psi.3.c57 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "c57" & dat$new == 0], y = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "c57" & dat$new == 1], alternative = "two.sided")
condt.psi.1.c57.tab <- c(condt.psi.1.c57$estimate[1], condt.psi.1.c57$estimate[2], diff = condt.psi.1.c57$estimate[1] - condt.psi.1.c57$estimate[2], stderr = condt.psi.1.c57$stderr, n1 = (n.old.psi.c57[2]), n2 = (n.new.psi.c57[2]), condt.psi.1.c57$statistic, condt.psi.1.c57$parameter, condt.psi.1.c57$p.value, condt.psi.1.c57$conf.int[1], condt.psi.1.c57$conf.int[2])
names(condt.psi.1.c57.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.psi.3.c57.tab <- c(condt.psi.3.c57$estimate[1], condt.psi.3.c57$estimate[2], diff = condt.psi.3.c57$estimate[1] - condt.psi.3.c57$estimate[2], stderr = condt.psi.3.c57$stderr, n1 = (n.old.psi.c57[3]), n2 = (n.new.psi.c57[3]), condt.psi.3.c57$statistic, condt.psi.3.c57$parameter, condt.psi.3.c57$p.value, condt.psi.3.c57$conf.int[1], condt.psi.3.c57$conf.int[2])
names(condt.psi.3.c57.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: pre-existing vs. newly-formed spines, compare within-treatment (saline) over time, wt
condt.sal.1.c57 <- t.test(x = dat$day_1[dat$treatment == "saline" & dat$genotype == "c57" & dat$new == 0], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "c57" & dat$new == 1], alternative = "two.sided")
condt.sal.3.c57 <- t.test(x = dat$day_3[dat$treatment == "saline" & dat$genotype == "c57" & dat$new == 0], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "c57" & dat$new == 1], alternative = "two.sided")
condt.sal.1.c57.tab <- c(condt.sal.1.c57$estimate[1], condt.sal.1.c57$estimate[2], diff = condt.sal.1.c57$estimate[1] - condt.sal.1.c57$estimate[2], stderr = condt.sal.1.c57$stderr, n1 = (n.old.sal.c57[2]), n2 = (n.new.sal.c57[2]), condt.sal.1.c57$statistic, condt.sal.1.c57$parameter, condt.sal.1.c57$p.value, condt.sal.1.c57$conf.int[1], condt.sal.1.c57$conf.int[2])
names(condt.sal.1.c57.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.sal.3.c57.tab <- c(condt.sal.3.c57$estimate[1], condt.sal.3.c57$estimate[2], diff = condt.sal.3.c57$estimate[1] - condt.sal.3.c57$estimate[2], stderr = condt.sal.3.c57$stderr, n1 = (n.old.sal.c57[3]), n2 = (n.new.sal.c57[3]), condt.sal.3.c57$statistic, condt.sal.3.c57$parameter, condt.sal.3.c57$p.value, condt.sal.3.c57$conf.int[1], condt.sal.3.c57$conf.int[2])
names(condt.sal.3.c57.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: pre-existing vs. newly-formed spines, compare within-treatment (psilocybin) over time, ko
condt.psi.1.ko <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "ko" & dat$new == 0], y = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "ko" & dat$new == 1], alternative = "two.sided")
condt.psi.3.ko <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "ko" & dat$new == 0], y = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "ko" & dat$new == 1], alternative = "two.sided")
condt.psi.1.ko.tab <- c(condt.psi.1.ko$estimate[1], condt.psi.1.ko$estimate[2], diff = condt.psi.1.ko$estimate[1] - condt.psi.1.ko$estimate[2], stderr = condt.psi.1.ko$stderr, n1 = (n.old.psi.ko[2]), n2 = (n.new.psi.ko[2]), condt.psi.1.ko$statistic, condt.psi.1.ko$parameter, condt.psi.1.ko$p.value, condt.psi.1.ko$conf.int[1], condt.psi.1.ko$conf.int[2])
names(condt.psi.1.ko.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.psi.3.ko.tab <- c(condt.psi.3.ko$estimate[1], condt.psi.3.ko$estimate[2], diff = condt.psi.3.ko$estimate[1] - condt.psi.3.ko$estimate[2], stderr = condt.psi.3.ko$stderr, n1 = (n.old.psi.ko[3]), n2 = (n.new.psi.ko[3]), condt.psi.3.ko$statistic, condt.psi.3.ko$parameter, condt.psi.3.ko$p.value, condt.psi.3.ko$conf.int[1], condt.psi.3.ko$conf.int[2])
names(condt.psi.3.ko.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: pre-existing vs. newly-formed spines, compare within-treatment (saline) over time, ko
condt.sal.1.ko <- t.test(x = dat$day_1[dat$treatment == "saline" & dat$genotype == "ko" & dat$new == 0], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "ko" & dat$new == 1], alternative = "two.sided")
condt.sal.3.ko <- t.test(x = dat$day_3[dat$treatment == "saline" & dat$genotype == "ko" & dat$new == 0], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "ko" & dat$new == 1], alternative = "two.sided")
condt.sal.1.ko.tab <- c(condt.sal.1.ko$estimate[1], condt.sal.1.ko$estimate[2], diff = condt.sal.1.ko$estimate[1] - condt.sal.1.ko$estimate[2], stderr = condt.sal.1.ko$stderr, n1 = (n.old.sal.ko[2]), n2 = (n.new.sal.ko[2]), condt.sal.1.ko$statistic, condt.sal.1.ko$parameter, condt.sal.1.ko$p.value, condt.sal.1.ko$conf.int[1], condt.sal.1.ko$conf.int[2])
names(condt.sal.1.ko.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.sal.3.ko.tab <- c(condt.sal.3.ko$estimate[1], condt.sal.3.ko$estimate[2], diff = condt.sal.3.ko$estimate[1] - condt.sal.3.ko$estimate[2], stderr = condt.sal.3.ko$stderr, n1 = (n.old.sal.ko[3]), n2 = (n.new.sal.ko[3]), condt.sal.3.ko$statistic, condt.sal.3.ko$parameter, condt.sal.3.ko$p.value, condt.sal.3.ko$conf.int[1], condt.sal.3.ko$conf.int[2])
names(condt.sal.3.ko.tab) = c("pre-existing", "newly-formed", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")


# ACCUMULATE
cond.all.tab <- rbind(
  condt.psi.1.tab, condt.psi.3.tab, 
  condt.sal.1.tab, condt.sal.3.tab, 
  condt.psi.1.c57.tab, condt.psi.3.c57.tab, 
  condt.sal.1.c57.tab, condt.sal.3.c57.tab, 
  condt.psi.1.ko.tab, condt.psi.3.ko.tab, 
  condt.sal.1.ko.tab, condt.sal.3.ko.tab)
# print(cond.all.tab)

#####

write.csv(x = cond.all.tab, paste(output_dir, paste("gfp_2P_psilocybin_pt_5ht2ako_R_lme_post_hoc_ttests_spine_headwidth_existing_vs_new_rawvalues_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))





