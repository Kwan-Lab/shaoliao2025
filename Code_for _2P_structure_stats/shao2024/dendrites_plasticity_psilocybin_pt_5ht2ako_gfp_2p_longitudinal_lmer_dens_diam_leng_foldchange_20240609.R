## dendrites plasticity: 
## lingxiao longitudinal 2P imaging data of cg1/m2 pt neuron structure, psilocybin vs. saline, pt neurons in 5ht2a f/f mice (5ht2ako's)
# 
# measures: density, head width, head length
# 
# notes:
# run linear mixed effects models of spine density, diameter, and length (fold change)
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
output_dir <- paste(proj_dir, "output", "psilocybin_pt_it_gfp_2p_structure", "pt_5ht2ako_foldchange", sep = "/")
figure_dir <- paste(proj_dir, "figure", sep = "/")

## load data table(s)
dens  <- read.csv(paste(data_dir, "gfp_2p_psilocybin_pt_5ht2ako_spine_density_foldchange_20240608.csv", sep = "/"), header = TRUE)
diam  <- read.csv(paste(data_dir, "gfp_2p_psilocybin_pt_5ht2ako_spine_headwidth_foldchange_20240608.csv", sep = "/"), header = TRUE)
leng  <- read.csv(paste(data_dir, "gfp_2p_psilocybin_pt_5ht2ako_spine_headlength_foldchange_20240608.csv", sep = "/"), header = TRUE)

## select data and reformat
#   - table must be in long format (time columnized)
dat <- dens
dat.name <- "density"
dat$mouse_id <- factor(sub(" ", "_", dat$mouse_id))
dat$genotype <- factor(dat$genotype)
dat$cell_id <- factor(paste(dat$mouse_id, dat$cell_id, sep = "_"))
dat$cell_id[grepl("NA", dat$cell_id)] = NA
dat$branch_id <- factor(paste(dat$mouse_id, dat$branch_id, sep = "_"))
dat$branch_id[grepl("NA", dat$branch_id)] = NA
dat$multibranch_id <- factor(paste(dat$mouse_id, dat$multibranch_id, sep = "_"))
dat$multibranch_id[grepl("NA", dat$multibranch_id)] = NA
dat$sex <- factor(dat$sex)
dat$treatment <- factor(dat$treatment)
dat$soma_depth[grepl("\\+", dat$soma_depth)] = NA
dat$soma_depth <- as.numeric(dat$soma_depth)
dat$fov <- factor(dat$fov)
datLong <- melt(data = dat,
                id.vars = c("mouse_id", "sex", "treatment", "genotype", "fov", "soma_depth", "cell_id", "branch_id", "multibranch_id"),
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
mod3 <- lmer(value ~ treatment * time * genotype + (1 | mouse_id / cell_id / branch_id), data = datLong.tcont)
anova(mod3)
# summary(mod3)
# plot(mod3) # visual check on residual linearity & variance
# qqnorm(residuals(mod3)) # visual check on normality of residuals
# fixef(mod3)
# ranef(mod3)
# coef(mod3)
write.csv(x = anova(mod3), paste(output_dir, paste("gfp_2P_psilocybin_pt_5ht2ako_R_lme_model_output_spine_", dat.name, "_foldchange_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))


# ## FULL post hoc t-tests comparison table, shorthand
# mm <- emmeans(mod3, specs = c("treatment", "time", "genotype", "sex"))
# cc <- summary(contrast(mm, method = "pairwise", adjust = "bonferroni"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ################################# Selected post hoc t-tests ################################# #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Get sample sizes for various groupings and sub-groupings
#####
n.psi = c(n.all = sum(dat$treatment == "psilocybin"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin"])))
n.psi.c57 = c(n.all = sum(dat$treatment == "psilocybin" & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "c57"])))
n.psi.ko = c(n.all = sum(dat$treatment == "psilocybin" & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "ko"])))
n.sal = c(n.all = sum(dat$treatment == "saline"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline"])))
n.sal.c57 = c(n.all = sum(dat$treatment == "saline" & dat$genotype == "c57"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$genotype == "c57"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$genotype == "c57"])))
n.sal.ko = c(n.all = sum(dat$treatment == "saline" & dat$genotype == "ko"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$genotype == "ko"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$genotype == "ko"])))

#####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# TWO-SAMPLE T-TESTS: treatment effect over time, altogether and split up by genotype
#####

## TWO-SAMPLE T-TESTS: treatment effect over time, agnostic to cell type, both sexes
condt.1 <- t.test(x = dat$day_1[dat$treatment == "psilocybin"], y = dat$day_1[dat$treatment == "saline"], alternative = "two.sided")
condt.3 <- t.test(x = dat$day_3[dat$treatment == "psilocybin"], y = dat$day_3[dat$treatment == "saline"], alternative = "two.sided")
condt.1.tab <- c(condt.1$estimate[1], condt.1$estimate[2], diff = condt.1$estimate[1] - condt.1$estimate[2], stderr = condt.1$stderr, n1 = n.psi[2], n2 = n.sal[2], condt.1$statistic, condt.1$parameter, condt.1$p.value, condt.1$conf.int[1], condt.1$conf.int[2])
names(condt.1.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.tab <- c(condt.3$estimate[1], condt.3$estimate[2], diff = condt.3$estimate[1] - condt.3$estimate[2], stderr = condt.3$stderr, n1 = n.psi[3], n2 = n.sal[3], condt.3$statistic, condt.3$parameter, condt.3$p.value, condt.3$conf.int[1], condt.3$conf.int[2])
names(condt.3.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effects over time for PT neurons, both sexes
condt.1.c57 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "c57"], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "c57"], alternative = "two.sided")
condt.3.c57 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "c57"], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "c57"], alternative = "two.sided")
condt.1.c57.tab <- c(condt.1.c57$estimate[1], condt.1.c57$estimate[2], diff = condt.1.c57$estimate[1] - condt.1.c57$estimate[2], stderr = condt.1.c57$stderr, n1 = (n.psi.c57[2]), n2 = (n.sal.c57[2]), condt.1.c57$statistic, condt.1.c57$parameter, condt.1.c57$p.value, condt.1.c57$conf.int[1], condt.1.c57$conf.int[2])
names(condt.1.c57.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.c57.tab <- c(condt.3.c57$estimate[1], condt.3.c57$estimate[2], diff = condt.3.c57$estimate[1] - condt.3.c57$estimate[2], stderr = condt.3.c57$stderr, n1 = (n.psi.c57[3]), n2 = (n.sal.c57[3]), condt.3.c57$statistic, condt.3.c57$parameter, condt.3.c57$p.value, condt.3.c57$conf.int[1], condt.3.c57$conf.int[2])
names(condt.3.c57.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effects over time for it neurons, both sexes
condt.1.ko <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "ko"], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "ko"], alternative = "two.sided")
condt.3.ko <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "ko"], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "ko"], alternative = "two.sided")
condt.1.ko.tab <- c(condt.1.ko$estimate[1], condt.1.ko$estimate[2], diff = condt.1.ko$estimate[1] - condt.1.ko$estimate[2], stderr = condt.1.ko$stderr, n1 = (n.psi.ko[2]), n2 = (n.sal.ko[2]), condt.1.ko$statistic, condt.1.ko$parameter, condt.1.ko$p.value, condt.1.ko$conf.int[1], condt.1.ko$conf.int[2])
names(condt.1.ko.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.ko.tab <- c(condt.3.ko$estimate[1], condt.3.ko$estimate[2], diff = condt.3.ko$estimate[1] - condt.3.ko$estimate[2], stderr = condt.3.ko$stderr, n1 = (n.psi.ko[3]), n2 = (n.sal.ko[3]), condt.3.ko$statistic, condt.3.ko$parameter, condt.3.ko$p.value, condt.3.ko$conf.int[1], condt.3.ko$conf.int[2])
names(condt.3.ko.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

# ACCUMULATE
cond.all.tab <- rbind(
  condt.1.tab, condt.3.tab, 
  condt.1.c57.tab, condt.3.c57.tab, 
  condt.1.ko.tab, condt.3.ko.tab)
# print(cond.all.tab)

#####

write.csv(x = cond.all.tab, paste(output_dir, paste("gfp_2P_psilocybin_pt_5ht2ako_R_lme_post_hoc_ttests_spine_", dat.name, "_foldchange_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# TWO-SAMPLE T-TESTS: genotype-based difference for same treatment, both sexes / females / males
#####

# TWO-SAMPLE T-TESTS: cell type-based difference for same treatment, both sexes
condt.1.c57_vs_ko.psi <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "c57"], y = dat$day_1[dat$treatment == "psilocybin" & dat$genotype == "ko"], alternative = "two.sided")
condt.3.c57_vs_ko.psi <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "c57"], y = dat$day_3[dat$treatment == "psilocybin" & dat$genotype == "ko"], alternative = "two.sided")
condt.1.c57_vs_ko.sal <- t.test(x = dat$day_1[dat$treatment == "saline" & dat$genotype == "c57"], y = dat$day_1[dat$treatment == "saline" & dat$genotype == "ko"], alternative = "two.sided")
condt.3.c57_vs_ko.sal <- t.test(x = dat$day_3[dat$treatment == "saline" & dat$genotype == "c57"], y = dat$day_3[dat$treatment == "saline" & dat$genotype == "ko"], alternative = "two.sided")
condt.1.c57_vs_ko.psi.tab <- c(condt.1.c57_vs_ko.psi$estimate[1], condt.1.c57_vs_ko.psi$estimate[2], diff = condt.1.c57_vs_ko.psi$estimate[1] - condt.1.c57_vs_ko.psi$estimate[2], stderr = condt.1.c57_vs_ko.psi$stderr, n1 = n.psi.c57[2], n2 = n.psi.ko[2], condt.1.c57_vs_ko.psi$statistic, condt.1.c57_vs_ko.psi$parameter, condt.1.c57_vs_ko.psi$p.value, condt.1.c57_vs_ko.psi$conf.int[1], condt.1.c57_vs_ko.psi$conf.int[2])
names(condt.1.c57_vs_ko.psi.tab) = c("c57", "ko", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.c57_vs_ko.psi.tab <- c(condt.3.c57_vs_ko.psi$estimate[1], condt.3.c57_vs_ko.psi$estimate[2], diff = condt.3.c57_vs_ko.psi$estimate[1] - condt.3.c57_vs_ko.psi$estimate[2], stderr = condt.3.c57_vs_ko.psi$stderr, n1 = n.psi.c57[3], n2 = n.psi.ko[3], condt.3.c57_vs_ko.psi$statistic, condt.3.c57_vs_ko.psi$parameter, condt.3.c57_vs_ko.psi$p.value, condt.3.c57_vs_ko.psi$conf.int[1], condt.3.c57_vs_ko.psi$conf.int[2])
names(condt.3.c57_vs_ko.psi.tab) = c("c57", "ko", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.1.c57_vs_ko.sal.tab <- c(condt.1.c57_vs_ko.sal$estimate[1], condt.1.c57_vs_ko.sal$estimate[2], diff = condt.1.c57_vs_ko.sal$estimate[1] - condt.1.c57_vs_ko.sal$estimate[2], stderr = condt.1.c57_vs_ko.sal$stderr, n1 = n.sal.c57[2], n2 = n.sal.ko[2], condt.1.c57_vs_ko.sal$statistic, condt.1.c57_vs_ko.sal$parameter, condt.1.c57_vs_ko.sal$p.value, condt.1.c57_vs_ko.sal$conf.int[1], condt.1.c57_vs_ko.sal$conf.int[2])
names(condt.1.c57_vs_ko.sal.tab) = c("c57", "ko", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.c57_vs_ko.sal.tab <- c(condt.3.c57_vs_ko.sal$estimate[1], condt.3.c57_vs_ko.sal$estimate[2], diff = condt.3.c57_vs_ko.sal$estimate[1] - condt.3.c57_vs_ko.sal$estimate[2], stderr = condt.3.c57_vs_ko.sal$stderr, n1 = n.sal.c57[3], n2 = n.sal.ko[3], condt.3.c57_vs_ko.sal$statistic, condt.3.c57_vs_ko.sal$parameter, condt.3.c57_vs_ko.sal$p.value, condt.3.c57_vs_ko.sal$conf.int[1], condt.3.c57_vs_ko.sal$conf.int[2])
names(condt.3.c57_vs_ko.sal.tab) = c("c57", "ko", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

cond.c57_vs_ko.tab <- rbind(
  condt.1.c57_vs_ko.psi.tab, condt.3.c57_vs_ko.psi.tab, 
  condt.1.c57_vs_ko.sal.tab, condt.3.c57_vs_ko.sal.tab)
# print(cond.c57_vs_ko.tab)

#####

write.csv(x = cond.c57_vs_ko.tab, paste(output_dir, paste("gfp_2P_psilocybin_pt_5ht2ako_R_lme_post_hoc_ttests_spine_", dat.name, "_foldchange_c57_vs_ko_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #







