## dendrites plasticity: 
## lingxiao longitudinal 2P imaging data of cg1/m2 pt & it neuron structure, psilocybin vs. saline, pt vs. it neurons
# 
# measures: density, head width, head length RAW VALUES
# 
# notes:
# run linear mixed effects models of spine density, diameter, and length
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
data_dir <- paste(proj_dir, "data", "psilocybin_pt_it_gfp_2p_structure", "20240610_pt_vs_it", sep = "/")
script_dir <- paste(proj_dir, "scripts", sep = "/")
output_dir <- paste(proj_dir, "output", "psilocybin_pt_it_gfp_2p_structure", "pt_vs_it_rawvalues", sep = "/")
figure_dir <- paste(proj_dir, "figure", sep = "/")

## load data table(s)
dens  <- read.csv(paste(data_dir, "gfp_2p_psilocybin_pt_it_spine_density_rawvalues_20240610.csv", sep = "/"), header = TRUE)
diam  <- read.csv(paste(data_dir, "gfp_2p_psilocybin_pt_it_spine_headwidth_rawvalues_20240610.csv", sep = "/"), header = TRUE)
leng   <- read.csv(paste(data_dir, "gfp_2p_psilocybin_pt_it_spine_headlength_rawvalues_20240610.csv", sep = "/"), header = TRUE)

## select data and reformat
#   - table must be in long format (time columnized)
dat <- dens
dat.name <- "density"
dat$cell_type <- factor(dat$cell_type)
dat$cell_id <- factor(paste(dat$mouse_id, dat$cell_id, sep = "_"))
dat$cell_id[grepl("NA", dat$cell_id)] = NA
dat$branch_id <- factor(paste(dat$mouse_id, dat$branch_id, sep = "_"))
dat$branch_id[grepl("NA", dat$branch_id)] = NA
dat$multibranch_id <- factor(paste(dat$mouse_id, dat$multibranch_id, sep = "_"))
dat$multibranch_id[grepl("NA", dat$multibranch_id)] = NA
dat$mouse_id <- factor(dat$mouse_id)
dat$sex <- factor(dat$sex)
dat$treatment <- factor(dat$treatment)
dat$soma_depth[grepl("\\+", dat$soma_depth)] = NA
dat$soma_depth <- as.numeric(dat$soma_depth)
dat$fov <- factor(dat$fov)
datLong <- melt(data = dat,
                id.vars = c("mouse_id", "sex", "treatment", "cell_type", "fov", "soma_depth", "cell_id", "branch_id", "multibranch_id"),
                measure.vars = c("day_1", "day_3", "day_5", "day_7", "day_35", "day_65"),
                variable.name = "time",
                value.name = "value",
                variable.factor = TRUE)
# tapply(datLong$value, datLong$sex, mean, na.rm = TRUE)
# tapply(datLong$value, datLong$treatment, mean, na.rm = TRUE)
# tapply(datLong$value, datLong$time, mean, na.rm = TRUE)
# tapply(datLong$value, datLong$cell_type, mean, na.rm = TRUE)
datLong_pt <- subset(datLong, cell_type == "pt")
datLong_it <- subset(datLong, cell_type == "it")
datLong.tcont <- datLong
levels(datLong.tcont$time)[levels(datLong.tcont$time) == "day_1"] <- as.numeric(1)
levels(datLong.tcont$time)[levels(datLong.tcont$time) == "day_3"] <- as.numeric(3)
levels(datLong.tcont$time)[levels(datLong.tcont$time) == "day_5"] <- as.numeric(5)
levels(datLong.tcont$time)[levels(datLong.tcont$time) == "day_7"] <- as.numeric(7)
levels(datLong.tcont$time)[levels(datLong.tcont$time) == "day_35"] <- as.numeric(35)
levels(datLong.tcont$time)[levels(datLong.tcont$time) == "day_65"] <- as.numeric(65)
datLong.tcont$time <- as.numeric(paste(datLong.tcont$time))

## LME - value ~ treatment * time * sex * cell_type + controlling for nested random effects of mouse_id / cell_id / multibranch_id
# mod3 <- lmer(value ~ treatment * time * sex * cell_type + (1 | mouse_id / cell_id / multibranch_id), data = datLong.tcont)  # multibranch_id drops info about known cells that we only measure a single branch, but these should be treated separately from NAs with no cell provenance info.! 06/11/24
mod3 <- lmer(value ~ treatment * time * sex * cell_type + (1 | mouse_id / cell_id / branch_id), data = datLong.tcont)
anova(mod3)
# summary(mod3)
# plot(mod3) # visual check on residual linearity & variance
# qqnorm(residuals(mod3)) # visual check on normality of residuals
# fixef(mod3)
# ranef(mod3)
# coef(mod3)
write.csv(x = anova(mod3), paste(output_dir, paste("gfp_2P_psilocybin_R_lme_model_output_spine_", dat.name, "_rawvalues_all_neurons_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))

## LME - IT neuron values ~ treatment * time * sex * soma_depth * + controlling for nested random effects of mouse_id / cell_id / multibranch_id
mod3_it <- lmer(value ~ treatment * time * soma_depth * sex + (1 | mouse_id / cell_id / branch_id), data = datLong_it)
anova(mod3_it)
# summary(mod3_it)
# plot(mod3_it) # visual check on residual linearity & variance
# qqnorm(residuals(mod3_it)) # visual check on normality of residuals
# fixef(mod3_it)
# ranef(mod3_it)
# coef(mod3_it)
write.csv(x = anova(mod3_it), paste(output_dir, paste("gfp_2P_psilocybin_it_only_soma_depth_R_lme_model_output_spine_", dat.name, "_foldchange_all_neurons_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))
# cor.test(x=datLong_it$soma_depth, y=datLong_it$value)
# plot(x=datLong_it$soma_depth, y=datLong_it$value)
# mean.by.depth <- aggregate(datLong_it$value, list(datLong_it$soma_depth), FUN=mean)
# plot(mean.by.depth)
# ggplot(datLong_it[datLong_it$time=="day_65",], aes(x = soma_depth, y = value, color = factor(treatment))) +
#   geom_point(size = 1, alpha = 0.5) +
#   scale_color_manual(values = c("Red", "Black")) +
#   ylim(-.25, .5) +
#   theme_bw() +
#   theme(legend.position="bottom") +
#   geom_smooth(method = lm, fullrange = TRUE) +
#   theme(panel.border = element_blank(), axis.line = element_line()) +
#   theme(panel.border = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5, size = 16), axis.text = element_text(size = 12), axis.title = element_text(size = 16)) +
#   labs(title = paste("IT neuron spine ", dat.name, " by depth, day 65", sep =""), x = "soma depth (um)", y = dat.name)
# ggplot(datLong_it[], aes(x = factor(soma_depth > 400), y = value, color = factor(treatment))) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitterdodge(jitter.width = 0.5), size = 1, alpha = 0.5) +
#   scale_color_manual(values = c("Red", "Black")) +
#   ylim(-.35, .5) +
#   theme_bw() +
#   theme(legend.position="bottom") +
#   theme(panel.border = element_blank(), axis.line = element_line()) +
#   theme(panel.border = element_blank(), axis.line = element_line(), plot.title = element_text(hjust = 0.5, size = 16), axis.text = element_text(size = 12), axis.title = element_text(size = 16)) +
#   labs(title = paste(dat.name, " by Depth", sep =""), x = "Depth > 400um", y = dat.name)

## LME - PT neuron values ~ treatment * time * sex + controlling for nested random effects of mouse_id / cell_id / multibranch_id
# mod3_pt <- lmer(value ~ treatment * time * sex + (1 | mouse_id / cell_id / multibranch_id), data = datLong_pt)
# anova(mod3_pt)
# summary(mod3_pt)
# plot(mod3_pt) # visual check on residual linearity & variance
# qqnorm(residuals(mod3_pt)) # visual check on normality of residuals
# fixef(mod3_pt)
# ranef(mod3_pt)
# coef(mod3_pt)

# ## FULL post hoc t-tests comparison table, shorthand
# mm <- emmeans(mod3, specs = c("treatment", "time", "cell_type", "sex"))
# cc <- summary(contrast(mm, method = "pairwise", adjust = "bonferroni"))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ################################# Selected post hoc t-tests ################################# #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Get sample sizes for various groupings and sub-groupings
#####
n.psi = c(n.all = sum(dat$treatment == "psilocybin"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin"])))
n.psi.f = c(n.all = sum(dat$treatment == "psilocybin" & dat$sex == "f"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$sex == "f"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$sex == "f"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin" & dat$sex == "f"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin" & dat$sex == "f"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin" & dat$sex == "f"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin" & dat$sex == "f"])))
n.psi.m = c(n.all = sum(dat$treatment == "psilocybin" & dat$sex == "m"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$sex == "m"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$sex == "m"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin" & dat$sex == "m"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin" & dat$sex == "m"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin" & dat$sex == "m"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin" & dat$sex == "m"])))
n.psi.pt = c(n.all = sum(dat$treatment == "psilocybin" & dat$cell_type == "pt"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "pt"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "pt"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "pt"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "pt"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "pt"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "pt"])))
n.psi.it = c(n.all = sum(dat$treatment == "psilocybin" & dat$cell_type == "it"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "it"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "it"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "it"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "it"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "it"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "it"])))
n.psi.f.pt = c(n.all = sum(dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"])))
n.psi.f.it = c(n.all = sum(dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"])))
n.psi.m.pt = c(n.all = sum(dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"])))
n.psi.m.it = c(n.all = sum(dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"), n1 = sum(!is.na(dat$day_1[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"])))
n.sal = c(n.all = sum(dat$treatment == "saline"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline"])))
n.sal.f = c(n.all = sum(dat$treatment == "saline" & dat$sex == "f"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$sex == "f"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$sex == "f"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline" & dat$sex == "f"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline" & dat$sex == "f"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline" & dat$sex == "f"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline" & dat$sex == "f"])))
n.sal.m = c(n.all = sum(dat$treatment == "saline" & dat$sex == "m"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$sex == "m"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$sex == "m"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline" & dat$sex == "m"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline" & dat$sex == "m"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline" & dat$sex == "m"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline" & dat$sex == "m"])))
n.sal.pt = c(n.all = sum(dat$treatment == "saline" & dat$cell_type == "pt"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$cell_type == "pt"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$cell_type == "pt"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline" & dat$cell_type == "pt"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline" & dat$cell_type == "pt"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline" & dat$cell_type == "pt"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline" & dat$cell_type == "pt"])))
n.sal.it = c(n.all = sum(dat$treatment == "saline" & dat$cell_type == "it"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$cell_type == "it"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$cell_type == "it"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline" & dat$cell_type == "it"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline" & dat$cell_type == "it"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline" & dat$cell_type == "it"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline" & dat$cell_type == "it"])))
n.sal.f.pt = c(n.all = sum(dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"])))
n.sal.f.it = c(n.all = sum(dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"])))
n.sal.m.pt = c(n.all = sum(dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"])))
n.sal.m.it = c(n.all = sum(dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"), n1 = sum(!is.na(dat$day_1[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"])), n3 = sum(!is.na(dat$day_3[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"])), n5 = sum(!is.na(dat$day_5[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"])), n7 = sum(!is.na(dat$day_7[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"])), n35 = sum(!is.na(dat$day_35[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"])), n65 = sum(!is.na(dat$day_65[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"])))

#####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# TWO-SAMPLE T-TESTS: treatment effect over time, agnostic to cell type, both sexes / females / males
#####

## TWO-SAMPLE T-TESTS: treatment effect over time, agnostic to cell type, both sexes
condt.1 <- t.test(x = dat$day_1[dat$treatment == "psilocybin"], y = dat$day_1[dat$treatment == "saline"], alternative = "two.sided")
condt.3 <- t.test(x = dat$day_3[dat$treatment == "psilocybin"], y = dat$day_3[dat$treatment == "saline"], alternative = "two.sided")
condt.5 <- t.test(x = dat$day_5[dat$treatment == "psilocybin"], y = dat$day_5[dat$treatment == "saline"], alternative = "two.sided")
condt.7 <- t.test(x = dat$day_7[dat$treatment == "psilocybin"], y = dat$day_7[dat$treatment == "saline"], alternative = "two.sided")
condt.35 <- t.test(x = dat$day_35[dat$treatment == "psilocybin"], y = dat$day_35[dat$treatment == "saline"], alternative = "two.sided")
condt.65 <- t.test(x = dat$day_65[dat$treatment == "psilocybin"], y = dat$day_65[dat$treatment == "saline"], alternative = "two.sided")
condt.1.tab <- c(condt.1$estimate[1], condt.1$estimate[2], diff = condt.1$estimate[1] - condt.1$estimate[2], stderr = condt.1$stderr, n1 = n.psi[2], n2 = n.sal[2], condt.1$statistic, condt.1$parameter, condt.1$p.value, condt.1$conf.int[1], condt.1$conf.int[2])
names(condt.1.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.tab <- c(condt.3$estimate[1], condt.3$estimate[2], diff = condt.3$estimate[1] - condt.3$estimate[2], stderr = condt.3$stderr, n1 = n.psi[3], n2 = n.sal[3], condt.3$statistic, condt.3$parameter, condt.3$p.value, condt.3$conf.int[1], condt.3$conf.int[2])
names(condt.3.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.5.tab <- c(condt.5$estimate[1], condt.5$estimate[2], diff = condt.5$estimate[1] - condt.5$estimate[2], stderr = condt.5$stderr, n1 = n.psi[4], n2 = n.sal[4], condt.5$statistic, condt.5$parameter, condt.5$p.value, condt.5$conf.int[1], condt.5$conf.int[2])
names(condt.5.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.7.tab <- c(condt.7$estimate[1], condt.7$estimate[2], diff = condt.7$estimate[1] - condt.7$estimate[2], stderr = condt.7$stderr, n1 = n.psi[5], n2 = n.sal[5], condt.7$statistic, condt.7$parameter, condt.7$p.value, condt.7$conf.int[1], condt.7$conf.int[2])
names(condt.7.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.35.tab <- c(condt.35$estimate[1], condt.35$estimate[2], diff = condt.35$estimate[1] - condt.35$estimate[2], stderr = condt.35$stderr, n1 = n.psi[6], n2 = n.sal[6], condt.35$statistic, condt.35$parameter, condt.35$p.value, condt.35$conf.int[1], condt.35$conf.int[2])
names(condt.35.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.65.tab <- c(condt.65$estimate[1], condt.65$estimate[2], diff = condt.65$estimate[1] - condt.65$estimate[2], stderr = condt.65$stderr, n1 = n.psi[7], n2 = n.sal[7], condt.65$statistic, condt.65$parameter, condt.65$p.value, condt.65$conf.int[1], condt.65$conf.int[2])
names(condt.65.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effect over time, agnostic to cell type, FEMALES
Fcondt.1 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$sex == "f"], y = dat$day_1[dat$treatment == "saline" & dat$sex == "f"], alternative = "two.sided")
Fcondt.3 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$sex == "f"], y = dat$day_3[dat$treatment == "saline" & dat$sex == "f"], alternative = "two.sided")
Fcondt.5 <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$sex == "f"], y = dat$day_5[dat$treatment == "saline" & dat$sex == "f"], alternative = "two.sided")
Fcondt.7 <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$sex == "f"], y = dat$day_7[dat$treatment == "saline" & dat$sex == "f"], alternative = "two.sided")
Fcondt.35 <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$sex == "f"], y = dat$day_35[dat$treatment == "saline" & dat$sex == "f"], alternative = "two.sided")
Fcondt.65 <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$sex == "f"], y = dat$day_65[dat$treatment == "saline" & dat$sex == "f"], alternative = "two.sided")
Fcondt.1.tab <- c(Fcondt.1$estimate[1], Fcondt.1$estimate[2], diff = Fcondt.1$estimate[1] - Fcondt.1$estimate[2], stderr = Fcondt.1$stderr, n1 = n.psi.f[2], n2 = n.sal.f[2], Fcondt.1$statistic, Fcondt.1$parameter, Fcondt.1$p.value, Fcondt.1$conf.int[1], Fcondt.1$conf.int[2])
names(Fcondt.1.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.3.tab <- c(Fcondt.3$estimate[1], Fcondt.3$estimate[2], diff = Fcondt.3$estimate[1] - Fcondt.3$estimate[2], stderr = Fcondt.3$stderr, n1 = n.psi.f[3], n2 = n.sal.f[3], Fcondt.3$statistic, Fcondt.3$parameter, Fcondt.3$p.value, Fcondt.3$conf.int[1], Fcondt.3$conf.int[2])
names(Fcondt.3.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.5.tab <- c(Fcondt.5$estimate[1], Fcondt.5$estimate[2], diff = Fcondt.5$estimate[1] - Fcondt.5$estimate[2], stderr = Fcondt.5$stderr, n1 = n.psi.f[4], n2 = n.sal.f[4], Fcondt.5$statistic, Fcondt.5$parameter, Fcondt.5$p.value, Fcondt.5$conf.int[1], Fcondt.5$conf.int[2])
names(Fcondt.5.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.7.tab <- c(Fcondt.7$estimate[1], Fcondt.7$estimate[2], diff = Fcondt.7$estimate[1] - Fcondt.7$estimate[2], stderr = Fcondt.7$stderr, n1 = n.psi.f[5], n2 = n.sal.f[5], Fcondt.7$statistic, Fcondt.7$parameter, Fcondt.7$p.value, Fcondt.7$conf.int[1], Fcondt.7$conf.int[2])
names(Fcondt.7.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.35.tab <- c(Fcondt.35$estimate[1], Fcondt.35$estimate[2], diff = Fcondt.35$estimate[1] - Fcondt.35$estimate[2], stderr = Fcondt.35$stderr, n1 = n.psi.f[6], n2 = n.sal.f[6], Fcondt.35$statistic, Fcondt.35$parameter, Fcondt.35$p.value, Fcondt.35$conf.int[1], Fcondt.35$conf.int[2])
names(Fcondt.35.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.65.tab <- c(Fcondt.65$estimate[1], Fcondt.65$estimate[2], diff = Fcondt.65$estimate[1] - Fcondt.65$estimate[2], stderr = Fcondt.65$stderr, n1 = n.psi.f[7], n2 = n.sal.f[7], Fcondt.65$statistic, Fcondt.65$parameter, Fcondt.65$p.value, Fcondt.65$conf.int[1], Fcondt.65$conf.int[2])
names(Fcondt.65.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## TWO-SAMPLE T-TESTS: treatment effect over time, agnostic to cell type, MALES
Mcondt.1 <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$sex == "m"], y = dat$day_1[dat$treatment == "saline" & dat$sex == "m"], alternative = "two.sided")
Mcondt.3 <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$sex == "m"], y = dat$day_3[dat$treatment == "saline" & dat$sex == "m"], alternative = "two.sided")
Mcondt.5 <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$sex == "m"], y = dat$day_5[dat$treatment == "saline" & dat$sex == "m"], alternative = "two.sided")
Mcondt.7 <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$sex == "m"], y = dat$day_7[dat$treatment == "saline" & dat$sex == "m"], alternative = "two.sided")
Mcondt.35 <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$sex == "m"], y = dat$day_35[dat$treatment == "saline" & dat$sex == "m"], alternative = "two.sided")
Mcondt.65 <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$sex == "m"], y = dat$day_65[dat$treatment == "saline" & dat$sex == "m"], alternative = "two.sided")
Mcondt.1.tab <- c(Mcondt.1$estimate[1], Mcondt.1$estimate[2], diff = Mcondt.1$estimate[1] - Mcondt.1$estimate[2], stderr = Mcondt.1$stderr, n1 = n.psi.m[2], n2 = n.sal.m[2], Mcondt.1$statistic, Mcondt.1$parameter, Mcondt.1$p.value, Mcondt.1$conf.int[1], Mcondt.1$conf.int[2])
names(Mcondt.1.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.3.tab <- c(Mcondt.3$estimate[1], Mcondt.3$estimate[2], diff = Mcondt.3$estimate[1] - Mcondt.3$estimate[2], stderr = Mcondt.3$stderr, n1 = n.psi.m[3], n2 = n.sal.m[3], Mcondt.3$statistic, Mcondt.3$parameter, Mcondt.3$p.value, Mcondt.3$conf.int[1], Mcondt.3$conf.int[2])
names(Mcondt.3.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.5.tab <- c(Mcondt.5$estimate[1], Mcondt.5$estimate[2], diff = Mcondt.5$estimate[1] - Mcondt.5$estimate[2], stderr = Mcondt.5$stderr, n1 = n.psi.m[4], n2 = n.sal.m[4], Mcondt.5$statistic, Mcondt.5$parameter, Mcondt.5$p.value, Mcondt.5$conf.int[1], Mcondt.5$conf.int[2])
names(Mcondt.5.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.7.tab <- c(Mcondt.7$estimate[1], Mcondt.7$estimate[2], diff = Mcondt.7$estimate[1] - Mcondt.7$estimate[2], stderr = Mcondt.7$stderr, n1 = n.psi.m[5], n2 = n.sal.m[5], Mcondt.7$statistic, Mcondt.7$parameter, Mcondt.7$p.value, Mcondt.7$conf.int[1], Mcondt.7$conf.int[2])
names(Mcondt.7.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.35.tab <- c(Mcondt.35$estimate[1], Mcondt.35$estimate[2], diff = Mcondt.35$estimate[1] - Mcondt.35$estimate[2], stderr = Mcondt.35$stderr, n1 = n.psi.m[6], n2 = n.sal.m[6], Mcondt.35$statistic, Mcondt.35$parameter, Mcondt.35$p.value, Mcondt.35$conf.int[1], Mcondt.35$conf.int[2])
names(Mcondt.35.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.65.tab <- c(Mcondt.65$estimate[1], Mcondt.65$estimate[2], diff = Mcondt.65$estimate[1] - Mcondt.65$estimate[2], stderr = Mcondt.65$stderr, n1 = n.psi.m[7], n2 = n.sal.m[7], Mcondt.65$statistic, Mcondt.65$parameter, Mcondt.65$p.value, Mcondt.65$conf.int[1], Mcondt.65$conf.int[2])
names(Mcondt.65.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

# ACCUMULATE
cond.all.tab <- rbind(
  condt.1.tab, condt.3.tab, condt.5.tab, condt.7.tab, condt.35.tab, condt.65.tab,
  Fcondt.1.tab, Fcondt.3.tab, Fcondt.5.tab, Fcondt.7.tab, Fcondt.35.tab, Fcondt.65.tab,
  Mcondt.1.tab, Mcondt.3.tab, Mcondt.5.tab, Mcondt.7.tab, Mcondt.35.tab, Mcondt.65.tab)
# print(cond.all.tab)

#####

write.csv(x = cond.all.tab, paste(output_dir, paste("gfp_2P_psilocybin_R_lme_post_hoc_ttests_spine_", dat.name, "_rawvalues_all_neurons_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# TWO-SAMPLE T-TESTS: treatment effect over time, PT neurons only, both sexes / females / males
#####

## TWO-SAMPLE T-TESTS: treatment effects over time for PT neurons, both sexes
condt.1.pt <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_1[dat$treatment == "saline" & dat$cell_type == "pt"], alternative = "two.sided")
condt.3.pt <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_3[dat$treatment == "saline" & dat$cell_type == "pt"], alternative = "two.sided")
condt.5.pt <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_5[dat$treatment == "saline" & dat$cell_type == "pt"], alternative = "two.sided")
condt.7.pt <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_7[dat$treatment == "saline" & dat$cell_type == "pt"], alternative = "two.sided")
condt.35.pt <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_35[dat$treatment == "saline" & dat$cell_type == "pt"], alternative = "two.sided")
condt.65.pt <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_65[dat$treatment == "saline" & dat$cell_type == "pt"], alternative = "two.sided")
condt.1.pt.tab <- c(condt.1.pt$estimate[1], condt.1.pt$estimate[2], diff = condt.1.pt$estimate[1] - condt.1.pt$estimate[2], stderr = condt.1.pt$stderr, n1 = (n.psi.m.pt[2] + n.psi.f.pt[2]), n2 = (n.sal.m.pt[2] + n.sal.f.pt[2]), condt.1.pt$statistic, condt.1.pt$parameter, condt.1.pt$p.value, condt.1.pt$conf.int[1], condt.1.pt$conf.int[2])
names(condt.1.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.pt.tab <- c(condt.3.pt$estimate[1], condt.3.pt$estimate[2], diff = condt.3.pt$estimate[1] - condt.3.pt$estimate[2], stderr = condt.3.pt$stderr, n1 = (n.psi.m.pt[3] + n.psi.f.pt[3]), n2 = (n.sal.m.pt[3] + n.sal.f.pt[3]), condt.3.pt$statistic, condt.3.pt$parameter, condt.3.pt$p.value, condt.3.pt$conf.int[1], condt.3.pt$conf.int[2])
names(condt.3.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.5.pt.tab <- c(condt.5.pt$estimate[1], condt.5.pt$estimate[2], diff = condt.5.pt$estimate[1] - condt.5.pt$estimate[2], stderr = condt.5.pt$stderr, n1 = (n.psi.m.pt[4] + n.psi.f.pt[4]), n2 = (n.sal.m.pt[4] + n.sal.f.pt[4]), condt.5.pt$statistic, condt.5.pt$parameter, condt.5.pt$p.value, condt.5.pt$conf.int[1], condt.5.pt$conf.int[2])
names(condt.5.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.7.pt.tab <- c(condt.7.pt$estimate[1], condt.7.pt$estimate[2], diff = condt.7.pt$estimate[1] - condt.7.pt$estimate[2], stderr = condt.7.pt$stderr, n1 = (n.psi.m.pt[5] + n.psi.f.pt[5]), n2 = (n.sal.m.pt[5] + n.sal.f.pt[5]), condt.7.pt$statistic, condt.7.pt$parameter, condt.7.pt$p.value, condt.7.pt$conf.int[1], condt.7.pt$conf.int[2])
names(condt.7.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.35.pt.tab <- c(condt.35.pt$estimate[1], condt.35.pt$estimate[2], diff = condt.35.pt$estimate[1] - condt.35.pt$estimate[2], stderr = condt.35.pt$stderr, n1 = (n.psi.m.pt[6] + n.psi.f.pt[6]), n2 = (n.sal.m.pt[6] + n.sal.f.pt[6]), condt.35.pt$statistic, condt.35.pt$parameter, condt.35.pt$p.value, condt.35.pt$conf.int[1], condt.35.pt$conf.int[2])
names(condt.35.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.65.pt.tab <- c(condt.65.pt$estimate[1], condt.65.pt$estimate[2], diff = condt.65.pt$estimate[1] - condt.65.pt$estimate[2], stderr = condt.65.pt$stderr, n1 = (n.psi.m.pt[7] + n.psi.f.pt[7]), n2 = (n.sal.m.pt[7] + n.sal.f.pt[7]), condt.65.pt$statistic, condt.65.pt$parameter, condt.65.pt$p.value, condt.65.pt$conf.int[1], condt.65.pt$conf.int[2])
names(condt.65.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

# TWO-SAMPLE T-TESTS: treatment effects over time for PT neurons, FEMALES
Fcondt.1.pt <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"], y = dat$day_1[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"], alternative = "two.sided")
Fcondt.3.pt <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"], y = dat$day_3[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"], alternative = "two.sided")
Fcondt.5.pt <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"], y = dat$day_5[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"], alternative = "two.sided")
Fcondt.7.pt <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"], y = dat$day_7[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"], alternative = "two.sided")
Fcondt.35.pt <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"], y = dat$day_35[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"], alternative = "two.sided")
Fcondt.65.pt <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "pt"], y = dat$day_65[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "pt"], alternative = "two.sided")
Fcondt.1.pt.tab <- c(Fcondt.1.pt$estimate[1], Fcondt.1.pt$estimate[2], diff = Fcondt.1.pt$estimate[1] - Fcondt.1.pt$estimate[2], stderr = Fcondt.1.pt$stderr, n1 = n.psi.f.pt[2], n2 = n.sal.f.pt[2], Fcondt.1.pt$statistic, Fcondt.1.pt$parameter, Fcondt.1.pt$p.value, Fcondt.1.pt$conf.int[1], Fcondt.1.pt$conf.int[2])
names(Fcondt.1.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.3.pt.tab <- c(Fcondt.3.pt$estimate[1], Fcondt.3.pt$estimate[2], diff = Fcondt.3.pt$estimate[1] - Fcondt.3.pt$estimate[2], stderr = Fcondt.3.pt$stderr, n1 = n.psi.f.pt[3], n2 = n.sal.f.pt[3], Fcondt.3.pt$statistic, Fcondt.3.pt$parameter, Fcondt.3.pt$p.value, Fcondt.3.pt$conf.int[1], Fcondt.3.pt$conf.int[2])
names(Fcondt.3.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.5.pt.tab <- c(Fcondt.5.pt$estimate[1], Fcondt.5.pt$estimate[2], diff = Fcondt.5.pt$estimate[1] - Fcondt.5.pt$estimate[2], stderr = Fcondt.5.pt$stderr, n1 = n.psi.f.pt[4], n2 = n.sal.f.pt[4], Fcondt.5.pt$statistic, Fcondt.5.pt$parameter, Fcondt.5.pt$p.value, Fcondt.5.pt$conf.int[1], Fcondt.5.pt$conf.int[2])
names(Fcondt.5.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.7.pt.tab <- c(Fcondt.7.pt$estimate[1], Fcondt.7.pt$estimate[2], diff = Fcondt.7.pt$estimate[1] - Fcondt.7.pt$estimate[2], stderr = Fcondt.7.pt$stderr, n1 = n.psi.f.pt[5], n2 = n.sal.f.pt[5], Fcondt.7.pt$statistic, Fcondt.7.pt$parameter, Fcondt.7.pt$p.value, Fcondt.7.pt$conf.int[1], Fcondt.7.pt$conf.int[2])
names(Fcondt.7.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.35.pt.tab <- c(Fcondt.35.pt$estimate[1], Fcondt.35.pt$estimate[2], diff = Fcondt.35.pt$estimate[1] - Fcondt.35.pt$estimate[2], stderr = Fcondt.35.pt$stderr, n1 = n.psi.f.pt[6], n2 = n.sal.f.pt[6], Fcondt.35.pt$statistic, Fcondt.35.pt$parameter, Fcondt.35.pt$p.value, Fcondt.35.pt$conf.int[1], Fcondt.35.pt$conf.int[2])
names(Fcondt.35.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.65.pt.tab <- c(Fcondt.65.pt$estimate[1], Fcondt.65.pt$estimate[2], diff = Fcondt.65.pt$estimate[1] - Fcondt.65.pt$estimate[2], stderr = Fcondt.65.pt$stderr, n1 = n.psi.f.pt[7], n2 = n.sal.f.pt[7], Fcondt.65.pt$statistic, Fcondt.65.pt$parameter, Fcondt.65.pt$p.value, Fcondt.65.pt$conf.int[1], Fcondt.65.pt$conf.int[2])
names(Fcondt.65.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

# TWO-SAMPLE T-TESTS: treatment effects over time for PT neurons, MALES
Mcondt.1.pt <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"], y = dat$day_1[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"], alternative = "two.sided")
Mcondt.3.pt <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"], y = dat$day_3[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"], alternative = "two.sided")
Mcondt.5.pt <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"], y = dat$day_5[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"], alternative = "two.sided")
Mcondt.7.pt <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"], y = dat$day_7[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"], alternative = "two.sided")
Mcondt.35.pt <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"], y = dat$day_35[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"], alternative = "two.sided")
Mcondt.65.pt <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "pt"], y = dat$day_65[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "pt"], alternative = "two.sided")
Mcondt.1.pt.tab <- c(Mcondt.1.pt$estimate[1], Mcondt.1.pt$estimate[2], diff = Mcondt.1.pt$estimate[1] - Mcondt.1.pt$estimate[2], stderr = Mcondt.1.pt$stderr, n1 = n.psi.m.pt[2], n2 = n.sal.m.pt[2], Mcondt.1.pt$statistic, Mcondt.1.pt$parameter, Mcondt.1.pt$p.value, Mcondt.1.pt$conf.int[1], Mcondt.1.pt$conf.int[2])
names(Mcondt.1.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.3.pt.tab <- c(Mcondt.3.pt$estimate[1], Mcondt.3.pt$estimate[2], diff = Mcondt.3.pt$estimate[1] - Mcondt.3.pt$estimate[2], stderr = Mcondt.3.pt$stderr, n1 = n.psi.m.pt[3], n2 = n.sal.m.pt[3], Mcondt.3.pt$statistic, Mcondt.3.pt$parameter, Mcondt.3.pt$p.value, Mcondt.3.pt$conf.int[1], Mcondt.3.pt$conf.int[2])
names(Mcondt.3.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.5.pt.tab <- c(Mcondt.5.pt$estimate[1], Mcondt.5.pt$estimate[2], diff = Mcondt.5.pt$estimate[1] - Mcondt.5.pt$estimate[2], stderr = Mcondt.5.pt$stderr, n1 = n.psi.m.pt[4], n2 = n.sal.m.pt[4], Mcondt.5.pt$statistic, Mcondt.5.pt$parameter, Mcondt.5.pt$p.value, Mcondt.5.pt$conf.int[1], Mcondt.5.pt$conf.int[2])
names(Mcondt.5.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.7.pt.tab <- c(Mcondt.7.pt$estimate[1], Mcondt.7.pt$estimate[2], diff = Mcondt.7.pt$estimate[1] - Mcondt.7.pt$estimate[2], stderr = Mcondt.7.pt$stderr, n1 = n.psi.m.pt[5], n2 = n.sal.m.pt[5], Mcondt.7.pt$statistic, Mcondt.7.pt$parameter, Mcondt.7.pt$p.value, Mcondt.7.pt$conf.int[1], Mcondt.7.pt$conf.int[2])
names(Mcondt.7.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.35.pt.tab <- c(Mcondt.35.pt$estimate[1], Mcondt.35.pt$estimate[2], diff = Mcondt.35.pt$estimate[1] - Mcondt.35.pt$estimate[2], stderr = Mcondt.35.pt$stderr, n1 = n.psi.m.pt[6], n2 = n.sal.m.pt[6], Mcondt.35.pt$statistic, Mcondt.35.pt$parameter, Mcondt.35.pt$p.value, Mcondt.35.pt$conf.int[1], Mcondt.35.pt$conf.int[2])
names(Mcondt.35.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.65.pt.tab <- c(Mcondt.65.pt$estimate[1], Mcondt.65.pt$estimate[2], diff = Mcondt.65.pt$estimate[1] - Mcondt.65.pt$estimate[2], stderr = Mcondt.65.pt$stderr, n1 = n.psi.m.pt[7], n2 = n.sal.m.pt[7], Mcondt.65.pt$statistic, Mcondt.65.pt$parameter, Mcondt.65.pt$p.value, Mcondt.65.pt$conf.int[1], Mcondt.65.pt$conf.int[2])
names(Mcondt.65.pt.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

# ACCUMULATE
cond.pt.tab <- rbind(
  condt.1.pt.tab, condt.3.pt.tab, condt.5.pt.tab, condt.7.pt.tab, condt.35.pt.tab, condt.65.pt.tab,
  Fcondt.1.pt.tab, Fcondt.3.pt.tab, Fcondt.5.pt.tab, Fcondt.7.pt.tab, Fcondt.35.pt.tab, Fcondt.65.pt.tab,
  Mcondt.1.pt.tab, Mcondt.3.pt.tab, Mcondt.5.pt.tab, Mcondt.7.pt.tab, Mcondt.35.pt.tab, Mcondt.65.pt.tab)
# print(cond.pt.tab)

#####

write.csv(x = cond.pt.tab, paste(output_dir, paste("gfp_2P_psilocybin_R_lme_post_hoc_ttests_spine_", dat.name, "_rawvalues_pt_neurons_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# TWO-SAMPLE T-TESTS: treatment effects over time,  IT neurons only, both sexes / females / males
#####

## TWO-SAMPLE T-TESTS: treatment effects over time for it neurons, both sexes
condt.1.it <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "it"], y = dat$day_1[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.3.it <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "it"], y = dat$day_3[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.5.it <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "it"], y = dat$day_5[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.7.it <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "it"], y = dat$day_7[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.35.it <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "it"], y = dat$day_35[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.65.it <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "it"], y = dat$day_65[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.1.it.tab <- c(condt.1.it$estimate[1], condt.1.it$estimate[2], diff = condt.1.it$estimate[1] - condt.1.it$estimate[2], stderr = condt.1.it$stderr, n1 = (n.psi.m.it[2] + n.psi.f.it[2]), n2 = (n.sal.m.it[2] + n.sal.f.it[2]), condt.1.it$statistic, condt.1.it$parameter, condt.1.it$p.value, condt.1.it$conf.int[1], condt.1.it$conf.int[2])
names(condt.1.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.it.tab <- c(condt.3.it$estimate[1], condt.3.it$estimate[2], diff = condt.3.it$estimate[1] - condt.3.it$estimate[2], stderr = condt.3.it$stderr, n1 = (n.psi.m.it[3] + n.psi.f.it[3]), n2 = (n.sal.m.it[3] + n.sal.f.it[3]), condt.3.it$statistic, condt.3.it$parameter, condt.3.it$p.value, condt.3.it$conf.int[1], condt.3.it$conf.int[2])
names(condt.3.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.5.it.tab <- c(condt.5.it$estimate[1], condt.5.it$estimate[2], diff = condt.5.it$estimate[1] - condt.5.it$estimate[2], stderr = condt.5.it$stderr, n1 = (n.psi.m.it[4] + n.psi.f.it[4]), n2 = (n.sal.m.it[4] + n.sal.f.it[4]), condt.5.it$statistic, condt.5.it$parameter, condt.5.it$p.value, condt.5.it$conf.int[1], condt.5.it$conf.int[2])
names(condt.5.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.7.it.tab <- c(condt.7.it$estimate[1], condt.7.it$estimate[2], diff = condt.7.it$estimate[1] - condt.7.it$estimate[2], stderr = condt.7.it$stderr, n1 = (n.psi.m.it[5] + n.psi.f.it[5]), n2 = (n.sal.m.it[5] + n.sal.f.it[5]), condt.7.it$statistic, condt.7.it$parameter, condt.7.it$p.value, condt.7.it$conf.int[1], condt.7.it$conf.int[2])
names(condt.7.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.35.it.tab <- c(condt.35.it$estimate[1], condt.35.it$estimate[2], diff = condt.35.it$estimate[1] - condt.35.it$estimate[2], stderr = condt.35.it$stderr, n1 = (n.psi.m.it[6] + n.psi.f.it[6]), n2 = (n.sal.m.it[6] + n.sal.f.it[6]), condt.35.it$statistic, condt.35.it$parameter, condt.35.it$p.value, condt.35.it$conf.int[1], condt.35.it$conf.int[2])
names(condt.35.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.65.it.tab <- c(condt.65.it$estimate[1], condt.65.it$estimate[2], diff = condt.65.it$estimate[1] - condt.65.it$estimate[2], stderr = condt.65.it$stderr, n1 = (n.psi.m.it[7] + n.psi.f.it[7]), n2 = (n.sal.m.it[7] + n.sal.f.it[7]), condt.65.it$statistic, condt.65.it$parameter, condt.65.it$p.value, condt.65.it$conf.int[1], condt.65.it$conf.int[2])
names(condt.65.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

# TWO-SAMPLE T-TESTS: treatment effects over time for it neurons, FEMALES
Fcondt.1.it <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"], y = dat$day_1[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"], alternative = "two.sided")
Fcondt.3.it <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"], y = dat$day_3[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"], alternative = "two.sided")
Fcondt.5.it <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"], y = dat$day_5[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"], alternative = "two.sided")
Fcondt.7.it <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"], y = dat$day_7[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"], alternative = "two.sided")
Fcondt.35.it <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"], y = dat$day_35[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"], alternative = "two.sided")
Fcondt.65.it <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$sex == "f" & dat$cell_type == "it"], y = dat$day_65[dat$treatment == "saline" & dat$sex == "f" & dat$cell_type == "it"], alternative = "two.sided")
Fcondt.1.it.tab <- c(Fcondt.1.it$estimate[1], Fcondt.1.it$estimate[2], diff = Fcondt.1.it$estimate[1] - Fcondt.1.it$estimate[2], stderr = Fcondt.1.it$stderr, n1 = n.psi.f.it[2], n2 = n.sal.f.it[2], Fcondt.1.it$statistic, Fcondt.1.it$parameter, Fcondt.1.it$p.value, Fcondt.1.it$conf.int[1], Fcondt.1.it$conf.int[2])
names(Fcondt.1.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.3.it.tab <- c(Fcondt.3.it$estimate[1], Fcondt.3.it$estimate[2], diff = Fcondt.3.it$estimate[1] - Fcondt.3.it$estimate[2], stderr = Fcondt.3.it$stderr, n1 = n.psi.f.it[3], n2 = n.sal.f.it[3], Fcondt.3.it$statistic, Fcondt.3.it$parameter, Fcondt.3.it$p.value, Fcondt.3.it$conf.int[1], Fcondt.3.it$conf.int[2])
names(Fcondt.3.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.5.it.tab <- c(Fcondt.5.it$estimate[1], Fcondt.5.it$estimate[2], diff = Fcondt.5.it$estimate[1] - Fcondt.5.it$estimate[2], stderr = Fcondt.5.it$stderr, n1 = n.psi.f.it[4], n2 = n.sal.f.it[4], Fcondt.5.it$statistic, Fcondt.5.it$parameter, Fcondt.5.it$p.value, Fcondt.5.it$conf.int[1], Fcondt.5.it$conf.int[2])
names(Fcondt.5.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.7.it.tab <- c(Fcondt.7.it$estimate[1], Fcondt.7.it$estimate[2], diff = Fcondt.7.it$estimate[1] - Fcondt.7.it$estimate[2], stderr = Fcondt.7.it$stderr, n1 = n.psi.f.it[5], n2 = n.sal.f.it[5], Fcondt.7.it$statistic, Fcondt.7.it$parameter, Fcondt.7.it$p.value, Fcondt.7.it$conf.int[1], Fcondt.7.it$conf.int[2])
names(Fcondt.7.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.35.it.tab <- c(Fcondt.35.it$estimate[1], Fcondt.35.it$estimate[2], diff = Fcondt.35.it$estimate[1] - Fcondt.35.it$estimate[2], stderr = Fcondt.35.it$stderr, n1 = n.psi.f.it[6], n2 = n.sal.f.it[6], Fcondt.35.it$statistic, Fcondt.35.it$parameter, Fcondt.35.it$p.value, Fcondt.35.it$conf.int[1], Fcondt.35.it$conf.int[2])
names(Fcondt.35.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.65.it.tab <- c(Fcondt.65.it$estimate[1], Fcondt.65.it$estimate[2], diff = Fcondt.65.it$estimate[1] - Fcondt.65.it$estimate[2], stderr = Fcondt.65.it$stderr, n1 = n.psi.f.it[7], n2 = n.sal.f.it[7], Fcondt.65.it$statistic, Fcondt.65.it$parameter, Fcondt.65.it$p.value, Fcondt.65.it$conf.int[1], Fcondt.65.it$conf.int[2])
names(Fcondt.65.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

# TWO-SAMPLE T-TESTS: treatment effects over time for it neurons, MALES
Mcondt.1.it <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"], y = dat$day_1[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"], alternative = "two.sided")
Mcondt.3.it <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"], y = dat$day_3[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"], alternative = "two.sided")
Mcondt.5.it <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"], y = dat$day_5[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"], alternative = "two.sided")
Mcondt.7.it <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"], y = dat$day_7[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"], alternative = "two.sided")
Mcondt.35.it <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"], y = dat$day_35[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"], alternative = "two.sided")
Mcondt.65.it <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$sex == "m" & dat$cell_type == "it"], y = dat$day_65[dat$treatment == "saline" & dat$sex == "m" & dat$cell_type == "it"], alternative = "two.sided")
Mcondt.1.it.tab <- c(Mcondt.1.it$estimate[1], Mcondt.1.it$estimate[2], diff = Mcondt.1.it$estimate[1] - Mcondt.1.it$estimate[2], stderr = Mcondt.1.it$stderr, n1 = n.psi.m.it[2], n2 = n.sal.m.it[2], Mcondt.1.it$statistic, Mcondt.1.it$parameter, Mcondt.1.it$p.value, Mcondt.1.it$conf.int[1], Mcondt.1.it$conf.int[2])
names(Mcondt.1.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.3.it.tab <- c(Mcondt.3.it$estimate[1], Mcondt.3.it$estimate[2], diff = Mcondt.3.it$estimate[1] - Mcondt.3.it$estimate[2], stderr = Mcondt.3.it$stderr, n1 = n.psi.m.it[3], n2 = n.sal.m.it[3], Mcondt.3.it$statistic, Mcondt.3.it$parameter, Mcondt.3.it$p.value, Mcondt.3.it$conf.int[1], Mcondt.3.it$conf.int[2])
names(Mcondt.3.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.5.it.tab <- c(Mcondt.5.it$estimate[1], Mcondt.5.it$estimate[2], diff = Mcondt.5.it$estimate[1] - Mcondt.5.it$estimate[2], stderr = Mcondt.5.it$stderr, n1 = n.psi.m.it[4], n2 = n.sal.m.it[4], Mcondt.5.it$statistic, Mcondt.5.it$parameter, Mcondt.5.it$p.value, Mcondt.5.it$conf.int[1], Mcondt.5.it$conf.int[2])
names(Mcondt.5.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.7.it.tab <- c(Mcondt.7.it$estimate[1], Mcondt.7.it$estimate[2], diff = Mcondt.7.it$estimate[1] - Mcondt.7.it$estimate[2], stderr = Mcondt.7.it$stderr, n1 = n.psi.m.it[5], n2 = n.sal.m.it[5], Mcondt.7.it$statistic, Mcondt.7.it$parameter, Mcondt.7.it$p.value, Mcondt.7.it$conf.int[1], Mcondt.7.it$conf.int[2])
names(Mcondt.7.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.35.it.tab <- c(Mcondt.35.it$estimate[1], Mcondt.35.it$estimate[2], diff = Mcondt.35.it$estimate[1] - Mcondt.35.it$estimate[2], stderr = Mcondt.35.it$stderr, n1 = n.psi.m.it[6], n2 = n.sal.m.it[6], Mcondt.35.it$statistic, Mcondt.35.it$parameter, Mcondt.35.it$p.value, Mcondt.35.it$conf.int[1], Mcondt.35.it$conf.int[2])
names(Mcondt.35.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.65.it.tab <- c(Mcondt.65.it$estimate[1], Mcondt.65.it$estimate[2], diff = Mcondt.65.it$estimate[1] - Mcondt.65.it$estimate[2], stderr = Mcondt.65.it$stderr, n1 = n.psi.m.it[7], n2 = n.sal.m.it[7], Mcondt.65.it$statistic, Mcondt.65.it$parameter, Mcondt.65.it$p.value, Mcondt.65.it$conf.int[1], Mcondt.65.it$conf.int[2])
names(Mcondt.65.it.tab) = c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

## ACCUMULATE
cond.it.tab <- rbind(
  condt.1.it.tab, condt.3.it.tab, condt.5.it.tab, condt.7.it.tab, condt.35.it.tab, condt.65.it.tab,
  Fcondt.1.it.tab, Fcondt.3.it.tab, Fcondt.5.it.tab, Fcondt.7.it.tab, Fcondt.35.it.tab, Fcondt.65.it.tab,
  Mcondt.1.it.tab, Mcondt.3.it.tab, Mcondt.5.it.tab, Mcondt.7.it.tab, Mcondt.35.it.tab, Mcondt.65.it.tab)
# print(cond.it.tab)

#####

write.csv(x = cond.it.tab, paste(output_dir, paste("gfp_2P_psilocybin_R_lme_post_hoc_ttests_spine_", dat.name, "_rawvalues_it_neurons_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# TWO-SAMPLE T-TESTS: cell type-based difference for same treatment, both sexes / females / males
#####

# TWO-SAMPLE T-TESTS: cell type-based difference for same treatment, both sexes
condt.1.pt_vs_it.psi <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "it"], alternative = "two.sided")
condt.3.pt_vs_it.psi <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "it"], alternative = "two.sided")
condt.5.pt_vs_it.psi <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "it"], alternative = "two.sided")
condt.7.pt_vs_it.psi <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "it"], alternative = "two.sided")
condt.35.pt_vs_it.psi <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "it"], alternative = "two.sided")
condt.65.pt_vs_it.psi <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "pt"], y = dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "it"], alternative = "two.sided")
condt.1.pt_vs_it.sal <- t.test(x = dat$day_1[dat$treatment == "saline" & dat$cell_type == "pt"], y = dat$day_1[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.3.pt_vs_it.sal <- t.test(x = dat$day_3[dat$treatment == "saline" & dat$cell_type == "pt"], y = dat$day_3[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.5.pt_vs_it.sal <- t.test(x = dat$day_5[dat$treatment == "saline" & dat$cell_type == "pt"], y = dat$day_5[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.7.pt_vs_it.sal <- t.test(x = dat$day_7[dat$treatment == "saline" & dat$cell_type == "pt"], y = dat$day_7[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.35.pt_vs_it.sal <- t.test(x = dat$day_35[dat$treatment == "saline" & dat$cell_type == "pt"], y = dat$day_35[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.65.pt_vs_it.sal <- t.test(x = dat$day_65[dat$treatment == "saline" & dat$cell_type == "pt"], y = dat$day_65[dat$treatment == "saline" & dat$cell_type == "it"], alternative = "two.sided")
condt.1.pt_vs_it.psi.tab <- c(condt.1.pt_vs_it.psi$estimate[1], condt.1.pt_vs_it.psi$estimate[2], diff = condt.1.pt_vs_it.psi$estimate[1] - condt.1.pt_vs_it.psi$estimate[2], stderr = condt.1.pt_vs_it.psi$stderr, n1 = n.psi.pt[2], n2 = n.psi.it[2], condt.1.pt_vs_it.psi$statistic, condt.1.pt_vs_it.psi$parameter, condt.1.pt_vs_it.psi$p.value, condt.1.pt_vs_it.psi$conf.int[1], condt.1.pt_vs_it.psi$conf.int[2])
names(condt.1.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.pt_vs_it.psi.tab <- c(condt.3.pt_vs_it.psi$estimate[1], condt.3.pt_vs_it.psi$estimate[2], diff = condt.3.pt_vs_it.psi$estimate[1] - condt.3.pt_vs_it.psi$estimate[2], stderr = condt.3.pt_vs_it.psi$stderr, n1 = n.psi.pt[3], n2 = n.psi.it[3], condt.3.pt_vs_it.psi$statistic, condt.3.pt_vs_it.psi$parameter, condt.3.pt_vs_it.psi$p.value, condt.3.pt_vs_it.psi$conf.int[1], condt.3.pt_vs_it.psi$conf.int[2])
names(condt.3.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.5.pt_vs_it.psi.tab <- c(condt.5.pt_vs_it.psi$estimate[1], condt.5.pt_vs_it.psi$estimate[2], diff = condt.5.pt_vs_it.psi$estimate[1] - condt.5.pt_vs_it.psi$estimate[2], stderr = condt.5.pt_vs_it.psi$stderr, n1 = n.psi.pt[4], n2 = n.psi.it[4], condt.5.pt_vs_it.psi$statistic, condt.5.pt_vs_it.psi$parameter, condt.5.pt_vs_it.psi$p.value, condt.5.pt_vs_it.psi$conf.int[1], condt.5.pt_vs_it.psi$conf.int[2])
names(condt.5.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.7.pt_vs_it.psi.tab <- c(condt.7.pt_vs_it.psi$estimate[1], condt.7.pt_vs_it.psi$estimate[2], diff = condt.7.pt_vs_it.psi$estimate[1] - condt.7.pt_vs_it.psi$estimate[2], stderr = condt.7.pt_vs_it.psi$stderr, n1 = n.psi.pt[5], n2 = n.psi.it[5], condt.7.pt_vs_it.psi$statistic, condt.7.pt_vs_it.psi$parameter, condt.7.pt_vs_it.psi$p.value, condt.7.pt_vs_it.psi$conf.int[1], condt.7.pt_vs_it.psi$conf.int[2])
names(condt.7.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.35.pt_vs_it.psi.tab <- c(condt.35.pt_vs_it.psi$estimate[1], condt.35.pt_vs_it.psi$estimate[2], diff = condt.35.pt_vs_it.psi$estimate[1] - condt.35.pt_vs_it.psi$estimate[2], stderr = condt.35.pt_vs_it.psi$stderr, n1 = n.psi.pt[6], n2 = n.psi.it[6], condt.35.pt_vs_it.psi$statistic, condt.35.pt_vs_it.psi$parameter, condt.35.pt_vs_it.psi$p.value, condt.35.pt_vs_it.psi$conf.int[1], condt.35.pt_vs_it.psi$conf.int[2])
names(condt.35.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.65.pt_vs_it.psi.tab <- c(condt.65.pt_vs_it.psi$estimate[1], condt.65.pt_vs_it.psi$estimate[2], diff = condt.65.pt_vs_it.psi$estimate[1] - condt.65.pt_vs_it.psi$estimate[2], stderr = condt.65.pt_vs_it.psi$stderr, n1 = n.psi.pt[7], n2 = n.psi.it[7], condt.65.pt_vs_it.psi$statistic, condt.65.pt_vs_it.psi$parameter, condt.65.pt_vs_it.psi$p.value, condt.65.pt_vs_it.psi$conf.int[1], condt.65.pt_vs_it.psi$conf.int[2])
names(condt.65.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.1.pt_vs_it.sal.tab <- c(condt.1.pt_vs_it.sal$estimate[1], condt.1.pt_vs_it.sal$estimate[2], diff = condt.1.pt_vs_it.sal$estimate[1] - condt.1.pt_vs_it.sal$estimate[2], stderr = condt.1.pt_vs_it.sal$stderr, n1 = n.sal.pt[2], n2 = n.sal.it[2], condt.1.pt_vs_it.sal$statistic, condt.1.pt_vs_it.sal$parameter, condt.1.pt_vs_it.sal$p.value, condt.1.pt_vs_it.sal$conf.int[1], condt.1.pt_vs_it.sal$conf.int[2])
names(condt.1.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.3.pt_vs_it.sal.tab <- c(condt.3.pt_vs_it.sal$estimate[1], condt.3.pt_vs_it.sal$estimate[2], diff = condt.3.pt_vs_it.sal$estimate[1] - condt.3.pt_vs_it.sal$estimate[2], stderr = condt.3.pt_vs_it.sal$stderr, n1 = n.sal.pt[3], n2 = n.sal.it[3], condt.3.pt_vs_it.sal$statistic, condt.3.pt_vs_it.sal$parameter, condt.3.pt_vs_it.sal$p.value, condt.3.pt_vs_it.sal$conf.int[1], condt.3.pt_vs_it.sal$conf.int[2])
names(condt.3.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.5.pt_vs_it.sal.tab <- c(condt.5.pt_vs_it.sal$estimate[1], condt.5.pt_vs_it.sal$estimate[2], diff = condt.5.pt_vs_it.sal$estimate[1] - condt.5.pt_vs_it.sal$estimate[2], stderr = condt.5.pt_vs_it.sal$stderr, n1 = n.sal.pt[4], n2 = n.sal.it[4], condt.5.pt_vs_it.sal$statistic, condt.5.pt_vs_it.sal$parameter, condt.5.pt_vs_it.sal$p.value, condt.5.pt_vs_it.sal$conf.int[1], condt.5.pt_vs_it.sal$conf.int[2])
names(condt.5.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.7.pt_vs_it.sal.tab <- c(condt.7.pt_vs_it.sal$estimate[1], condt.7.pt_vs_it.sal$estimate[2], diff = condt.7.pt_vs_it.sal$estimate[1] - condt.7.pt_vs_it.sal$estimate[2], stderr = condt.7.pt_vs_it.sal$stderr, n1 = n.sal.pt[5], n2 = n.sal.it[5], condt.7.pt_vs_it.sal$statistic, condt.7.pt_vs_it.sal$parameter, condt.7.pt_vs_it.sal$p.value, condt.7.pt_vs_it.sal$conf.int[1], condt.7.pt_vs_it.sal$conf.int[2])
names(condt.7.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.35.pt_vs_it.sal.tab <- c(condt.35.pt_vs_it.sal$estimate[1], condt.35.pt_vs_it.sal$estimate[2], diff = condt.35.pt_vs_it.sal$estimate[1] - condt.35.pt_vs_it.sal$estimate[2], stderr = condt.35.pt_vs_it.sal$stderr, n1 = n.sal.pt[6], n2 = n.sal.it[6], condt.35.pt_vs_it.sal$statistic, condt.35.pt_vs_it.sal$parameter, condt.35.pt_vs_it.sal$p.value, condt.35.pt_vs_it.sal$conf.int[1], condt.35.pt_vs_it.sal$conf.int[2])
names(condt.35.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
condt.65.pt_vs_it.sal.tab <- c(condt.65.pt_vs_it.sal$estimate[1], condt.65.pt_vs_it.sal$estimate[2], diff = condt.65.pt_vs_it.sal$estimate[1] - condt.65.pt_vs_it.sal$estimate[2], stderr = condt.65.pt_vs_it.sal$stderr, n1 = n.sal.pt[7], n2 = n.sal.it[7], condt.65.pt_vs_it.sal$statistic, condt.65.pt_vs_it.sal$parameter, condt.65.pt_vs_it.sal$p.value, condt.65.pt_vs_it.sal$conf.int[1], condt.65.pt_vs_it.sal$conf.int[2])
names(condt.65.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

# TWO-SAMPLE T-TESTS: cell type-based difference for same treatment, FEMALES
Fcondt.1.pt_vs_it.psi <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.3.pt_vs_it.psi <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.5.pt_vs_it.psi <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.7.pt_vs_it.psi <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.35.pt_vs_it.psi <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.65.pt_vs_it.psi <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.1.pt_vs_it.sal <- t.test(x = dat$day_1[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_1[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.3.pt_vs_it.sal <- t.test(x = dat$day_3[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_3[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.5.pt_vs_it.sal <- t.test(x = dat$day_5[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_5[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.7.pt_vs_it.sal <- t.test(x = dat$day_7[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_7[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.35.pt_vs_it.sal <- t.test(x = dat$day_35[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_35[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.65.pt_vs_it.sal <- t.test(x = dat$day_65[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "f"], y = dat$day_65[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "f"], alternative = "two.sided")
Fcondt.1.pt_vs_it.psi.tab <- c(Fcondt.1.pt_vs_it.psi$estimate[1], Fcondt.1.pt_vs_it.psi$estimate[2], diff = Fcondt.1.pt_vs_it.psi$estimate[1] - Fcondt.1.pt_vs_it.psi$estimate[2], stderr = Fcondt.1.pt_vs_it.psi$stderr, n1 = n.psi.f.pt[2], n2 = n.psi.f.it[2], Fcondt.1.pt_vs_it.psi$statistic, Fcondt.1.pt_vs_it.psi$parameter, Fcondt.1.pt_vs_it.psi$p.value, Fcondt.1.pt_vs_it.psi$conf.int[1], Fcondt.1.pt_vs_it.psi$conf.int[2])
names(Fcondt.1.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.3.pt_vs_it.psi.tab <- c(Fcondt.3.pt_vs_it.psi$estimate[1], Fcondt.3.pt_vs_it.psi$estimate[2], diff = Fcondt.3.pt_vs_it.psi$estimate[1] - Fcondt.3.pt_vs_it.psi$estimate[2], stderr = Fcondt.3.pt_vs_it.psi$stderr, n1 = n.psi.f.pt[3], n2 = n.psi.f.it[3], Fcondt.3.pt_vs_it.psi$statistic, Fcondt.3.pt_vs_it.psi$parameter, Fcondt.3.pt_vs_it.psi$p.value, Fcondt.3.pt_vs_it.psi$conf.int[1], Fcondt.3.pt_vs_it.psi$conf.int[2])
names(Fcondt.3.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.5.pt_vs_it.psi.tab <- c(Fcondt.5.pt_vs_it.psi$estimate[1], Fcondt.5.pt_vs_it.psi$estimate[2], diff = Fcondt.5.pt_vs_it.psi$estimate[1] - Fcondt.5.pt_vs_it.psi$estimate[2], stderr = Fcondt.5.pt_vs_it.psi$stderr, n1 = n.psi.f.pt[4], n2 = n.psi.f.it[4], Fcondt.5.pt_vs_it.psi$statistic, Fcondt.5.pt_vs_it.psi$parameter, Fcondt.5.pt_vs_it.psi$p.value, Fcondt.5.pt_vs_it.psi$conf.int[1], Fcondt.5.pt_vs_it.psi$conf.int[2])
names(Fcondt.5.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.7.pt_vs_it.psi.tab <- c(Fcondt.7.pt_vs_it.psi$estimate[1], Fcondt.7.pt_vs_it.psi$estimate[2], diff = Fcondt.7.pt_vs_it.psi$estimate[1] - Fcondt.7.pt_vs_it.psi$estimate[2], stderr = Fcondt.7.pt_vs_it.psi$stderr, n1 = n.psi.f.pt[5], n2 = n.psi.f.it[5], Fcondt.7.pt_vs_it.psi$statistic, Fcondt.7.pt_vs_it.psi$parameter, Fcondt.7.pt_vs_it.psi$p.value, Fcondt.7.pt_vs_it.psi$conf.int[1], Fcondt.7.pt_vs_it.psi$conf.int[2])
names(Fcondt.7.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.35.pt_vs_it.psi.tab <- c(Fcondt.35.pt_vs_it.psi$estimate[1], Fcondt.35.pt_vs_it.psi$estimate[2], diff = Fcondt.35.pt_vs_it.psi$estimate[1] - Fcondt.35.pt_vs_it.psi$estimate[2], stderr = Fcondt.35.pt_vs_it.psi$stderr, n1 = n.psi.f.pt[6], n2 = n.psi.f.it[6], Fcondt.35.pt_vs_it.psi$statistic, Fcondt.35.pt_vs_it.psi$parameter, Fcondt.35.pt_vs_it.psi$p.value, Fcondt.35.pt_vs_it.psi$conf.int[1], Fcondt.35.pt_vs_it.psi$conf.int[2])
names(Fcondt.35.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.65.pt_vs_it.psi.tab <- c(Fcondt.65.pt_vs_it.psi$estimate[1], Fcondt.65.pt_vs_it.psi$estimate[2], diff = Fcondt.65.pt_vs_it.psi$estimate[1] - Fcondt.65.pt_vs_it.psi$estimate[2], stderr = Fcondt.65.pt_vs_it.psi$stderr, n1 = n.psi.f.pt[7], n2 = n.psi.f.it[7], Fcondt.65.pt_vs_it.psi$statistic, Fcondt.65.pt_vs_it.psi$parameter, Fcondt.65.pt_vs_it.psi$p.value, Fcondt.65.pt_vs_it.psi$conf.int[1], Fcondt.65.pt_vs_it.psi$conf.int[2])
names(Fcondt.65.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.1.pt_vs_it.sal.tab <- c(Fcondt.1.pt_vs_it.sal$estimate[1], Fcondt.1.pt_vs_it.sal$estimate[2], diff = Fcondt.1.pt_vs_it.sal$estimate[1] - Fcondt.1.pt_vs_it.sal$estimate[2], stderr = Fcondt.1.pt_vs_it.sal$stderr, n1 = n.sal.f.pt[2], n2 = n.sal.f.it[2], Fcondt.1.pt_vs_it.sal$statistic, Fcondt.1.pt_vs_it.sal$parameter, Fcondt.1.pt_vs_it.sal$p.value, Fcondt.1.pt_vs_it.sal$conf.int[1], Fcondt.1.pt_vs_it.sal$conf.int[2])
names(Fcondt.1.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.3.pt_vs_it.sal.tab <- c(Fcondt.3.pt_vs_it.sal$estimate[1], Fcondt.3.pt_vs_it.sal$estimate[2], diff = Fcondt.3.pt_vs_it.sal$estimate[1] - Fcondt.3.pt_vs_it.sal$estimate[2], stderr = Fcondt.3.pt_vs_it.sal$stderr, n1 = n.sal.f.pt[3], n2 = n.sal.f.it[3], Fcondt.3.pt_vs_it.sal$statistic, Fcondt.3.pt_vs_it.sal$parameter, Fcondt.3.pt_vs_it.sal$p.value, Fcondt.3.pt_vs_it.sal$conf.int[1], Fcondt.3.pt_vs_it.sal$conf.int[2])
names(Fcondt.3.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.5.pt_vs_it.sal.tab <- c(Fcondt.5.pt_vs_it.sal$estimate[1], Fcondt.5.pt_vs_it.sal$estimate[2], diff = Fcondt.5.pt_vs_it.sal$estimate[1] - Fcondt.5.pt_vs_it.sal$estimate[2], stderr = Fcondt.5.pt_vs_it.sal$stderr, n1 = n.sal.f.pt[4], n2 = n.sal.f.it[4], Fcondt.5.pt_vs_it.sal$statistic, Fcondt.5.pt_vs_it.sal$parameter, Fcondt.5.pt_vs_it.sal$p.value, Fcondt.5.pt_vs_it.sal$conf.int[1], Fcondt.5.pt_vs_it.sal$conf.int[2])
names(Fcondt.5.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.7.pt_vs_it.sal.tab <- c(Fcondt.7.pt_vs_it.sal$estimate[1], Fcondt.7.pt_vs_it.sal$estimate[2], diff = Fcondt.7.pt_vs_it.sal$estimate[1] - Fcondt.7.pt_vs_it.sal$estimate[2], stderr = Fcondt.7.pt_vs_it.sal$stderr, n1 = n.sal.f.pt[5], n2 = n.sal.f.it[5], Fcondt.7.pt_vs_it.sal$statistic, Fcondt.7.pt_vs_it.sal$parameter, Fcondt.7.pt_vs_it.sal$p.value, Fcondt.7.pt_vs_it.sal$conf.int[1], Fcondt.7.pt_vs_it.sal$conf.int[2])
names(Fcondt.7.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.35.pt_vs_it.sal.tab <- c(Fcondt.35.pt_vs_it.sal$estimate[1], Fcondt.35.pt_vs_it.sal$estimate[2], diff = Fcondt.35.pt_vs_it.sal$estimate[1] - Fcondt.35.pt_vs_it.sal$estimate[2], stderr = Fcondt.35.pt_vs_it.sal$stderr, n1 = n.sal.f.pt[6], n2 = n.sal.f.it[6], Fcondt.35.pt_vs_it.sal$statistic, Fcondt.35.pt_vs_it.sal$parameter, Fcondt.35.pt_vs_it.sal$p.value, Fcondt.35.pt_vs_it.sal$conf.int[1], Fcondt.35.pt_vs_it.sal$conf.int[2])
names(Fcondt.35.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Fcondt.65.pt_vs_it.sal.tab <- c(Fcondt.65.pt_vs_it.sal$estimate[1], Fcondt.65.pt_vs_it.sal$estimate[2], diff = Fcondt.65.pt_vs_it.sal$estimate[1] - Fcondt.65.pt_vs_it.sal$estimate[2], stderr = Fcondt.65.pt_vs_it.sal$stderr, n1 = n.sal.f.pt[7], n2 = n.sal.f.it[7], Fcondt.65.pt_vs_it.sal$statistic, Fcondt.65.pt_vs_it.sal$parameter, Fcondt.65.pt_vs_it.sal$p.value, Fcondt.65.pt_vs_it.sal$conf.int[1], Fcondt.65.pt_vs_it.sal$conf.int[2])
names(Fcondt.65.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")


# TWO-SAMPLE T-TESTS: cell type-based difference for same treatment, MALES
Mcondt.1.pt_vs_it.psi <- t.test(x = dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_1[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.3.pt_vs_it.psi <- t.test(x = dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_3[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.5.pt_vs_it.psi <- t.test(x = dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_5[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.7.pt_vs_it.psi <- t.test(x = dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_7[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.35.pt_vs_it.psi <- t.test(x = dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_35[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.65.pt_vs_it.psi <- t.test(x = dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_65[dat$treatment == "psilocybin" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.1.pt_vs_it.sal <- t.test(x = dat$day_1[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_1[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.3.pt_vs_it.sal <- t.test(x = dat$day_3[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_3[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.5.pt_vs_it.sal <- t.test(x = dat$day_5[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_5[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.7.pt_vs_it.sal <- t.test(x = dat$day_7[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_7[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.35.pt_vs_it.sal <- t.test(x = dat$day_35[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_35[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.65.pt_vs_it.sal <- t.test(x = dat$day_65[dat$treatment == "saline" & dat$cell_type == "pt" & dat$sex == "m"], y = dat$day_65[dat$treatment == "saline" & dat$cell_type == "it" & dat$sex == "m"], alternative = "two.sided")
Mcondt.1.pt_vs_it.psi.tab <- c(Mcondt.1.pt_vs_it.psi$estimate[1], Mcondt.1.pt_vs_it.psi$estimate[2], diff = Mcondt.1.pt_vs_it.psi$estimate[1] - Mcondt.1.pt_vs_it.psi$estimate[2], stderr = Mcondt.1.pt_vs_it.psi$stderr, n1 = n.psi.m.pt[2], n2 = n.psi.m.it[2], Mcondt.1.pt_vs_it.psi$statistic, Mcondt.1.pt_vs_it.psi$parameter, Mcondt.1.pt_vs_it.psi$p.value, Mcondt.1.pt_vs_it.psi$conf.int[1], Mcondt.1.pt_vs_it.psi$conf.int[2])
names(Mcondt.1.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.3.pt_vs_it.psi.tab <- c(Mcondt.3.pt_vs_it.psi$estimate[1], Mcondt.3.pt_vs_it.psi$estimate[2], diff = Mcondt.3.pt_vs_it.psi$estimate[1] - Mcondt.3.pt_vs_it.psi$estimate[2], stderr = Mcondt.3.pt_vs_it.psi$stderr, n1 = n.psi.m.pt[3], n2 = n.psi.m.it[3], Mcondt.3.pt_vs_it.psi$statistic, Mcondt.3.pt_vs_it.psi$parameter, Mcondt.3.pt_vs_it.psi$p.value, Mcondt.3.pt_vs_it.psi$conf.int[1], Mcondt.3.pt_vs_it.psi$conf.int[2])
names(Mcondt.3.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.5.pt_vs_it.psi.tab <- c(Mcondt.5.pt_vs_it.psi$estimate[1], Mcondt.5.pt_vs_it.psi$estimate[2], diff = Mcondt.5.pt_vs_it.psi$estimate[1] - Mcondt.5.pt_vs_it.psi$estimate[2], stderr = Mcondt.5.pt_vs_it.psi$stderr, n1 = n.psi.m.pt[4], n2 = n.psi.m.it[4], Mcondt.5.pt_vs_it.psi$statistic, Mcondt.5.pt_vs_it.psi$parameter, Mcondt.5.pt_vs_it.psi$p.value, Mcondt.5.pt_vs_it.psi$conf.int[1], Mcondt.5.pt_vs_it.psi$conf.int[2])
names(Mcondt.5.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.7.pt_vs_it.psi.tab <- c(Mcondt.7.pt_vs_it.psi$estimate[1], Mcondt.7.pt_vs_it.psi$estimate[2], diff = Mcondt.7.pt_vs_it.psi$estimate[1] - Mcondt.7.pt_vs_it.psi$estimate[2], stderr = Mcondt.7.pt_vs_it.psi$stderr, n1 = n.psi.m.pt[5], n2 = n.psi.m.it[5], Mcondt.7.pt_vs_it.psi$statistic, Mcondt.7.pt_vs_it.psi$parameter, Mcondt.7.pt_vs_it.psi$p.value, Mcondt.7.pt_vs_it.psi$conf.int[1], Mcondt.7.pt_vs_it.psi$conf.int[2])
names(Mcondt.7.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.35.pt_vs_it.psi.tab <- c(Mcondt.35.pt_vs_it.psi$estimate[1], Mcondt.35.pt_vs_it.psi$estimate[2], diff = Mcondt.35.pt_vs_it.psi$estimate[1] - Mcondt.35.pt_vs_it.psi$estimate[2], stderr = Mcondt.35.pt_vs_it.psi$stderr, n1 = n.psi.m.pt[6], n2 = n.psi.m.it[6], Mcondt.35.pt_vs_it.psi$statistic, Mcondt.35.pt_vs_it.psi$parameter, Mcondt.35.pt_vs_it.psi$p.value, Mcondt.35.pt_vs_it.psi$conf.int[1], Mcondt.35.pt_vs_it.psi$conf.int[2])
names(Mcondt.35.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.65.pt_vs_it.psi.tab <- c(Mcondt.65.pt_vs_it.psi$estimate[1], Mcondt.65.pt_vs_it.psi$estimate[2], diff = Mcondt.65.pt_vs_it.psi$estimate[1] - Mcondt.65.pt_vs_it.psi$estimate[2], stderr = Mcondt.65.pt_vs_it.psi$stderr, n1 = n.psi.m.pt[7], n2 = n.psi.m.it[7], Mcondt.65.pt_vs_it.psi$statistic, Mcondt.65.pt_vs_it.psi$parameter, Mcondt.65.pt_vs_it.psi$p.value, Mcondt.65.pt_vs_it.psi$conf.int[1], Mcondt.65.pt_vs_it.psi$conf.int[2])
names(Mcondt.65.pt_vs_it.psi.tab) = c("psilocybin.pt", "psilocybin.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.1.pt_vs_it.sal.tab <- c(Mcondt.1.pt_vs_it.sal$estimate[1], Mcondt.1.pt_vs_it.sal$estimate[2], diff = Mcondt.1.pt_vs_it.sal$estimate[1] - Mcondt.1.pt_vs_it.sal$estimate[2], stderr = Mcondt.1.pt_vs_it.sal$stderr, n1 = n.sal.m.pt[2], n2 = n.sal.m.it[2], Mcondt.1.pt_vs_it.sal$statistic, Mcondt.1.pt_vs_it.sal$parameter, Mcondt.1.pt_vs_it.sal$p.value, Mcondt.1.pt_vs_it.sal$conf.int[1], Mcondt.1.pt_vs_it.sal$conf.int[2])
names(Mcondt.1.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.3.pt_vs_it.sal.tab <- c(Mcondt.3.pt_vs_it.sal$estimate[1], Mcondt.3.pt_vs_it.sal$estimate[2], diff = Mcondt.3.pt_vs_it.sal$estimate[1] - Mcondt.3.pt_vs_it.sal$estimate[2], stderr = Mcondt.3.pt_vs_it.sal$stderr, n1 = n.sal.m.pt[3], n2 = n.sal.m.it[3], Mcondt.3.pt_vs_it.sal$statistic, Mcondt.3.pt_vs_it.sal$parameter, Mcondt.3.pt_vs_it.sal$p.value, Mcondt.3.pt_vs_it.sal$conf.int[1], Mcondt.3.pt_vs_it.sal$conf.int[2])
names(Mcondt.3.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.5.pt_vs_it.sal.tab <- c(Mcondt.5.pt_vs_it.sal$estimate[1], Mcondt.5.pt_vs_it.sal$estimate[2], diff = Mcondt.5.pt_vs_it.sal$estimate[1] - Mcondt.5.pt_vs_it.sal$estimate[2], stderr = Mcondt.5.pt_vs_it.sal$stderr, n1 = n.sal.m.pt[4], n2 = n.sal.m.it[4], Mcondt.5.pt_vs_it.sal$statistic, Mcondt.5.pt_vs_it.sal$parameter, Mcondt.5.pt_vs_it.sal$p.value, Mcondt.5.pt_vs_it.sal$conf.int[1], Mcondt.5.pt_vs_it.sal$conf.int[2])
names(Mcondt.5.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.7.pt_vs_it.sal.tab <- c(Mcondt.7.pt_vs_it.sal$estimate[1], Mcondt.7.pt_vs_it.sal$estimate[2], diff = Mcondt.7.pt_vs_it.sal$estimate[1] - Mcondt.7.pt_vs_it.sal$estimate[2], stderr = Mcondt.7.pt_vs_it.sal$stderr, n1 = n.sal.m.pt[5], n2 = n.sal.m.it[5], Mcondt.7.pt_vs_it.sal$statistic, Mcondt.7.pt_vs_it.sal$parameter, Mcondt.7.pt_vs_it.sal$p.value, Mcondt.7.pt_vs_it.sal$conf.int[1], Mcondt.7.pt_vs_it.sal$conf.int[2])
names(Mcondt.7.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.35.pt_vs_it.sal.tab <- c(Mcondt.35.pt_vs_it.sal$estimate[1], Mcondt.35.pt_vs_it.sal$estimate[2], diff = Mcondt.35.pt_vs_it.sal$estimate[1] - Mcondt.35.pt_vs_it.sal$estimate[2], stderr = Mcondt.35.pt_vs_it.sal$stderr, n1 = n.sal.m.pt[6], n2 = n.sal.m.it[6], Mcondt.35.pt_vs_it.sal$statistic, Mcondt.35.pt_vs_it.sal$parameter, Mcondt.35.pt_vs_it.sal$p.value, Mcondt.35.pt_vs_it.sal$conf.int[1], Mcondt.35.pt_vs_it.sal$conf.int[2])
names(Mcondt.35.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")
Mcondt.65.pt_vs_it.sal.tab <- c(Mcondt.65.pt_vs_it.sal$estimate[1], Mcondt.65.pt_vs_it.sal$estimate[2], diff = Mcondt.65.pt_vs_it.sal$estimate[1] - Mcondt.65.pt_vs_it.sal$estimate[2], stderr = Mcondt.65.pt_vs_it.sal$stderr, n1 = n.sal.m.pt[7], n2 = n.sal.m.it[7], Mcondt.65.pt_vs_it.sal$statistic, Mcondt.65.pt_vs_it.sal$parameter, Mcondt.65.pt_vs_it.sal$p.value, Mcondt.65.pt_vs_it.sal$conf.int[1], Mcondt.65.pt_vs_it.sal$conf.int[2])
names(Mcondt.65.pt_vs_it.sal.tab) = c("saline.pt", "saline.it", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi")

cond.pt_vs_it.tab <- rbind(
  condt.1.pt_vs_it.psi.tab, condt.3.pt_vs_it.psi.tab, condt.5.pt_vs_it.psi.tab, condt.7.pt_vs_it.psi.tab, condt.35.pt_vs_it.psi.tab, condt.65.pt_vs_it.psi.tab,
  condt.1.pt_vs_it.sal.tab, condt.3.pt_vs_it.sal.tab, condt.5.pt_vs_it.sal.tab, condt.7.pt_vs_it.sal.tab, condt.35.pt_vs_it.sal.tab, condt.65.pt_vs_it.sal.tab, 
  Fcondt.1.pt_vs_it.psi.tab, Fcondt.3.pt_vs_it.psi.tab, Fcondt.5.pt_vs_it.psi.tab, Fcondt.7.pt_vs_it.psi.tab, Fcondt.35.pt_vs_it.psi.tab, Fcondt.65.pt_vs_it.psi.tab,
  Fcondt.1.pt_vs_it.sal.tab, Fcondt.3.pt_vs_it.sal.tab, Fcondt.5.pt_vs_it.sal.tab, Fcondt.7.pt_vs_it.sal.tab, Fcondt.35.pt_vs_it.sal.tab, Fcondt.65.pt_vs_it.sal.tab, 
  Mcondt.1.pt_vs_it.psi.tab, Mcondt.3.pt_vs_it.psi.tab, Mcondt.5.pt_vs_it.psi.tab, Mcondt.7.pt_vs_it.psi.tab, Mcondt.35.pt_vs_it.psi.tab, Mcondt.65.pt_vs_it.psi.tab, 
  Mcondt.1.pt_vs_it.sal.tab, Mcondt.3.pt_vs_it.sal.tab, Mcondt.5.pt_vs_it.sal.tab, Mcondt.7.pt_vs_it.sal.tab, Mcondt.35.pt_vs_it.sal.tab, Mcondt.65.pt_vs_it.sal.tab)
# print(cond.pt_vs_it.tab)

#####

write.csv(x = cond.pt_vs_it.tab, paste(output_dir, paste("gfp_2P_psilocybin_R_lme_post_hoc_ttests_spine_", dat.name, "_rawvalues_pt_vs_it_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #






#####
# TEMP TRYING IT NEURON DEPTH DEPENDENCE 2-SAMPLE T-TESTS:
# Get sample sizes for various groupings and sub-groupings
#####
depth.thresh <- 400
datLong_it.lo <- subset(datLong_it, soma_depth < depth.thresh)
datLong_it.hi <- subset(datLong_it, soma_depth >= depth.thresh)
n.it.lo.psi = c(n.all = sum(datLong_it.lo$treatment == "psilocybin"), n1 = sum(datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_1"), n3 = sum(datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_3"), n5 = sum(datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_5"), n7 = sum(datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_7"), n35 = sum(datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_35"), n65 = sum(datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_65"))
n.it.hi.psi = c(n.all = sum(datLong_it.hi$treatment == "psilocybin"), n1 = sum(datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_1"), n3 = sum(datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_3"), n5 = sum(datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_5"), n7 = sum(datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_7"), n35 = sum(datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_35"), n65 = sum(datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_65"))
n.it.lo.sal = c(n.all = sum(datLong_it.lo$treatment == "saline"), n1 = sum(datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_1"), n3 = sum(datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_3"), n5 = sum(datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_5"), n7 = sum(datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_7"), n35 = sum(datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_35"), n65 = sum(datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_65"))
n.it.hi.sal = c(n.all = sum(datLong_it.hi$treatment == "saline"), n1 = sum(datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_1"), n3 = sum(datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_3"), n5 = sum(datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_5"), n7 = sum(datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_7"), n35 = sum(datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_35"), n65 = sum(datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_65"))
condt.lo <- t.test(x = datLong_it.lo$value[datLong_it.lo$treatment == "psilocybin"], y = datLong_it.lo$value[datLong_it.lo$treatment == "saline"], alternative = "two.sided")
condt.hi <- t.test(x = datLong_it.hi$value[datLong_it.hi$treatment == "psilocybin"], y = datLong_it.hi$value[datLong_it.hi$treatment == "saline"], alternative = "two.sided")
condt.lo.1 <- t.test(x = datLong_it.lo$value[datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_1"], y = datLong_it.lo$value[datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_1"], alternative = "two.sided")
condt.hi.1 <- t.test(x = datLong_it.hi$value[datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_1"], y = datLong_it.hi$value[datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_1"], alternative = "two.sided")
condt.lo.3 <- t.test(x = datLong_it.lo$value[datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_3"], y = datLong_it.lo$value[datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_3"], alternative = "two.sided")
condt.hi.3 <- t.test(x = datLong_it.hi$value[datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_3"], y = datLong_it.hi$value[datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_3"], alternative = "two.sided")
condt.lo.5 <- t.test(x = datLong_it.lo$value[datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_5"], y = datLong_it.lo$value[datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_5"], alternative = "two.sided")
condt.hi.5 <- t.test(x = datLong_it.hi$value[datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_5"], y = datLong_it.hi$value[datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_5"], alternative = "two.sided")
condt.lo.7 <- t.test(x = datLong_it.lo$value[datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_7"], y = datLong_it.lo$value[datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_7"], alternative = "two.sided")
condt.hi.7 <- t.test(x = datLong_it.hi$value[datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_7"], y = datLong_it.hi$value[datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_7"], alternative = "two.sided")
condt.lo.35 <- t.test(x = datLong_it.lo$value[datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_35"], y = datLong_it.lo$value[datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_35"], alternative = "two.sided")
condt.hi.35 <- t.test(x = datLong_it.hi$value[datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_35"], y = datLong_it.hi$value[datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_35"], alternative = "two.sided")
condt.lo.65 <- t.test(x = datLong_it.lo$value[datLong_it.lo$treatment == "psilocybin" & datLong_it.lo$time == "day_65"], y = datLong_it.lo$value[datLong_it.lo$treatment == "saline" & datLong_it.lo$time == "day_65"], alternative = "two.sided")
condt.hi.65 <- t.test(x = datLong_it.hi$value[datLong_it.hi$treatment == "psilocybin" & datLong_it.hi$time == "day_65"], y = datLong_it.hi$value[datLong_it.hi$treatment == "saline" & datLong_it.hi$time == "day_65"], alternative = "two.sided")
condt.lo.tab <- c(condt.lo$estimate[1], condt.lo$estimate[2], diff = condt.lo$estimate[1] - condt.lo$estimate[2], stderr = condt.lo$stderr, n1 = n.it.lo.psi[1], n2 = n.it.lo.sal[1], condt.lo$statistic, condt.lo$parameter, condt.lo$p.value, condt.lo$conf.int[1], condt.lo$conf.int[2])
condt.hi.tab <- c(condt.hi$estimate[1], condt.hi$estimate[2], diff = condt.hi$estimate[1] - condt.hi$estimate[2], stderr = condt.hi$stderr, n1 = n.it.hi.psi[1], n2 = n.it.hi.sal[1], condt.hi$statistic, condt.hi$parameter, condt.hi$p.value, condt.hi$conf.int[1], condt.hi$conf.int[2])
condt.lo.1.tab <- c(condt.lo.1$estimate[1], condt.lo.1$estimate[2], diff = condt.lo.1$estimate[1] - condt.lo.1$estimate[2], stderr = condt.lo.1$stderr, n1 = n.it.lo.psi[2], n2 = n.it.lo.sal[2], condt.lo.1$statistic, condt.lo.1$parameter, condt.lo.1$p.value, condt.lo.1$conf.int[1], condt.lo.1$conf.int[2])
condt.hi.1.tab <- c(condt.hi.1$estimate[1], condt.hi.1$estimate[2], diff = condt.hi.1$estimate[1] - condt.hi.1$estimate[2], stderr = condt.hi.1$stderr, n1 = n.it.hi.psi[2], n2 = n.it.hi.sal[2], condt.hi.1$statistic, condt.hi.1$parameter, condt.hi.1$p.value, condt.hi.1$conf.int[1], condt.hi.1$conf.int[2])
condt.lo.3.tab <- c(condt.lo.3$estimate[1], condt.lo.3$estimate[2], diff = condt.lo.3$estimate[1] - condt.lo.3$estimate[2], stderr = condt.lo.3$stderr, n1 = n.it.lo.psi[3], n2 = n.it.lo.sal[3], condt.lo.3$statistic, condt.lo.3$parameter, condt.lo.3$p.value, condt.lo.3$conf.int[1], condt.lo.3$conf.int[2])
condt.hi.3.tab <- c(condt.hi.3$estimate[1], condt.hi.3$estimate[2], diff = condt.hi.3$estimate[1] - condt.hi.3$estimate[2], stderr = condt.hi.3$stderr, n1 = n.it.hi.psi[3], n2 = n.it.hi.sal[3], condt.hi.3$statistic, condt.hi.3$parameter, condt.hi.3$p.value, condt.hi.3$conf.int[1], condt.hi.3$conf.int[2])
condt.lo.5.tab <- c(condt.lo.5$estimate[1], condt.lo.5$estimate[2], diff = condt.lo.5$estimate[1] - condt.lo.5$estimate[2], stderr = condt.lo.5$stderr, n1 = n.it.lo.psi[4], n2 = n.it.lo.sal[4], condt.lo.5$statistic, condt.lo.5$parameter, condt.lo.5$p.value, condt.lo.5$conf.int[1], condt.lo.5$conf.int[2])
condt.hi.5.tab <- c(condt.hi.5$estimate[1], condt.hi.5$estimate[2], diff = condt.hi.5$estimate[1] - condt.hi.5$estimate[2], stderr = condt.hi.5$stderr, n1 = n.it.hi.psi[4], n2 = n.it.hi.sal[4], condt.hi.5$statistic, condt.hi.5$parameter, condt.hi.5$p.value, condt.hi.5$conf.int[1], condt.hi.5$conf.int[2])
condt.lo.7.tab <- c(condt.lo.7$estimate[1], condt.lo.7$estimate[2], diff = condt.lo.7$estimate[1] - condt.lo.7$estimate[2], stderr = condt.lo.7$stderr, n1 = n.it.lo.psi[5], n2 = n.it.lo.sal[5], condt.lo.7$statistic, condt.lo.7$parameter, condt.lo.7$p.value, condt.lo.7$conf.int[1], condt.lo.7$conf.int[2])
condt.hi.7.tab <- c(condt.hi.7$estimate[1], condt.hi.7$estimate[2], diff = condt.hi.7$estimate[1] - condt.hi.7$estimate[2], stderr = condt.hi.7$stderr, n1 = n.it.hi.psi[5], n2 = n.it.hi.sal[5], condt.hi.7$statistic, condt.hi.7$parameter, condt.hi.7$p.value, condt.hi.7$conf.int[1], condt.hi.7$conf.int[2])
condt.lo.35.tab <- c(condt.lo.35$estimate[1], condt.lo.35$estimate[2], diff = condt.lo.35$estimate[1] - condt.lo.35$estimate[2], stderr = condt.lo.35$stderr, n1 = n.it.lo.psi[6], n2 = n.it.lo.sal[6], condt.lo.35$statistic, condt.lo.35$parameter, condt.lo.35$p.value, condt.lo.35$conf.int[1], condt.lo.35$conf.int[2])
condt.hi.35.tab <- c(condt.hi.35$estimate[1], condt.hi.35$estimate[2], diff = condt.hi.35$estimate[1] - condt.hi.35$estimate[2], stderr = condt.hi.35$stderr, n1 = n.it.hi.psi[6], n2 = n.it.hi.sal[6], condt.hi.35$statistic, condt.hi.35$parameter, condt.hi.35$p.value, condt.hi.35$conf.int[1], condt.hi.35$conf.int[2])
condt.lo.65.tab <- c(condt.lo.65$estimate[1], condt.lo.65$estimate[2], diff = condt.lo.65$estimate[1] - condt.lo.65$estimate[2], stderr = condt.lo.65$stderr, n1 = n.it.lo.psi[7], n2 = n.it.lo.sal[7], condt.lo.65$statistic, condt.lo.65$parameter, condt.lo.65$p.value, condt.lo.65$conf.int[1], condt.lo.65$conf.int[2])
condt.hi.65.tab <- c(condt.hi.65$estimate[1], condt.hi.65$estimate[2], diff = condt.hi.65$estimate[1] - condt.hi.65$estimate[2], stderr = condt.hi.65$stderr, n1 = n.it.hi.psi[7], n2 = n.it.hi.sal[7], condt.hi.65$statistic, condt.hi.65$parameter, condt.hi.65$p.value, condt.hi.65$conf.int[1], condt.hi.65$conf.int[2])

cond.it.depth.tab <- rbind(c("psilocybin", "saline", "diff", "stderr", "n1", "n2", "t", "df", "pvalue", "diff.conf95.lo", "diff.conf95.hi"), 
  condt.lo.tab, condt.lo.1.tab, condt.lo.3.tab, condt.lo.5.tab, condt.lo.7.tab, condt.lo.35.tab, condt.lo.65.tab,
  condt.hi.tab, condt.hi.1.tab, condt.hi.3.tab, condt.hi.5.tab, condt.hi.7.tab, condt.hi.35.tab, condt.hi.65.tab)
# print(cond.it.depth.tab)

write.csv(x = cond.it.depth.tab, paste(output_dir, paste("gfp_2P_psilocybin_R_lme_post_hoc_ttests_spine_", dat.name, "_rawvalues_it_only_soma_depth_thresh_", depth.thresh, "_", format(Sys.time(), "%Y%m%d"), ".csv", sep = ""), sep = "/"))







