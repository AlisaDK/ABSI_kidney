library(openxlsx)
library(readr)
library(readxl)
library(vroom)
library(dplyr)
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(forestplot)
library(ggplot2)
library(patchwork)

`%notin%` <- Negate(`%in%`)

#################################################
#Genetic instruments: ABSI, WHI & HI 
#################################################
allometric <- read.csv("allometric_exposures/TABS9_snps_2021_04_09.csv") # https://www.nature.com/articles/s41598-021-89176-6#data-availability
distinct(as.data.frame(allometric$Index))

ABSI_w <- allometric %>%
  filter(Index=="ABSIUKB" & Sex=="Women" & gwasP<5e-8 & SexDiff==0) %>% 
  mutate(var=paste(chr, pos, sep=":"))
ABSI_m <- allometric %>%
  filter(Index=="ABSIUKB" & Sex=="Men" & gwasP<5e-8 & SexDiff==0) %>% 
  mutate(var=paste(chr, pos, sep=":"))

WHI_w <- allometric %>%
  filter(Index=="WHIUKB" & Sex=="Women" & gwasP<5e-8 & SexDiff==0) %>% 
  mutate(var=paste(chr, pos, sep=":"))
WHI_m <- allometric %>%
  filter(Index=="WHIUKB" & Sex=="Men" & gwasP<5e-8 & SexDiff==0) %>% 
  mutate(var=paste(chr, pos, sep=":"))

HI_w <- allometric %>%
  filter(Index=="HIUKB" & Sex=="Women" & gwasP<5e-8 & SexDiff==0) %>% 
  mutate(var=paste(chr, pos, sep=":"))
HI_m <- allometric %>%
  filter(Index=="HIUKB" & Sex=="Men" & gwasP<5e-8 & SexDiff==0) %>% 
  mutate(var=paste(chr, pos, sep=":"))


names(ABSI_w)
dat_exposure_w <- data.frame(SNP = ABSI_w$SNP, beta.exposure= ABSI_w$beta, se.exposure= ABSI_w$se, id.exposure=1, exposure=ABSI_w$Index, sex="Women",
                             units.exposure="SD", pval.exposure=ABSI_w$gwasP, samplesize.exposure=219872, var=ABSI_w$var,
                             effect_allele.exposure= ABSI_w$effect_allele, other_allele.exposure= ABSI_w$non_effect_allele, eaf.exposure= ABSI_w$MAF)
ABSI_w <- dat_exposure_w %>% 
  clump_data(., clump_r2 = 0.01) 

# men
dat_exposure_m <- data.frame(SNP = ABSI_m$SNP, beta.exposure= ABSI_m$beta, se.exposure= ABSI_m$se, id.exposure=2, exposure=ABSI_m$Index, sex="Men",
                             units.exposure="SD", pval.exposure=ABSI_m$gwasP, samplesize.exposure=186825, var=ABSI_m$var, 
                             effect_allele.exposure= ABSI_m$effect_allele, other_allele.exposure= ABSI_m$non_effect_allele, eaf.exposure= ABSI_m$MAF)
ABSI_m <- dat_exposure_m %>%   
  clump_data(., clump_r2 = 0.01) 


dat_exposure_w <- data.frame(SNP = WHI_w$SNP, beta.exposure= WHI_w$beta, se.exposure= WHI_w$se, id.exposure=3, exposure=WHI_w$Index, sex="Women",
                             units.exposure="SD", pval.exposure=WHI_w$gwasP, samplesize.exposure=219872, var=WHI_w$var, 
                             effect_allele.exposure= WHI_w$effect_allele, other_allele.exposure= WHI_w$non_effect_allele, eaf.exposure= WHI_w$MAF)
WHI_w <- dat_exposure_w %>%   
  clump_data(., clump_r2 = 0.01) 

# men
dat_exposure_m <- data.frame(SNP = WHI_m$SNP, beta.exposure= WHI_m$beta, se.exposure= WHI_m$se, id.exposure=4, exposure=WHI_m$Index, sex="Men",
                             units.exposure="SD", pval.exposure=WHI_m$gwasP, samplesize.exposure=186825, var=WHI_m$var, 
                             effect_allele.exposure= WHI_m$effect_allele, other_allele.exposure= WHI_m$non_effect_allele, eaf.exposure= WHI_m$MAF)
WHI_m <- dat_exposure_m %>%   
  clump_data(., clump_r2 = 0.01) 


dat_exposure_w <- data.frame(SNP = HI_w$SNP, beta.exposure= HI_w$beta, se.exposure= HI_w$se, id.exposure=5, exposure=HI_w$Index, sex="Women",
                             units.exposure="SD", pval.exposure=HI_w$gwasP, samplesize.exposure=219872, var=HI_w$var, 
                             effect_allele.exposure= HI_w$effect_allele, other_allele.exposure= HI_w$non_effect_allele, eaf.exposure= HI_w$MAF)
HI_w <- dat_exposure_w %>%   
  clump_data(., clump_r2 = 0.01) 

# men
dat_exposure_m <- data.frame(SNP = HI_m$SNP, beta.exposure= HI_m$beta, se.exposure= HI_m$se, id.exposure=6, exposure=HI_m$Index, sex="Men",
                             units.exposure="SD", pval.exposure=HI_m$gwasP, samplesize.exposure=186825, var=HI_m$var, 
                             effect_allele.exposure= HI_m$effect_allele, other_allele.exposure= HI_m$non_effect_allele, eaf.exposure= HI_m$MAF)
HI_m <- dat_exposure_m %>%   
  clump_data(., clump_r2 = 0.01) 

#################################################
#Combining EXPOSURES, EA
#################################################
#Generate sex-stratified exposures
dat_exposure_w <- bind_rows(ABSI_w, WHI_w, HI_w) %>%
  mutate(exposure=paste(exposure, "_w"))
dat_exposure_m <- bind_rows(ABSI_m, WHI_m, HI_m) %>%
  mutate(exposure=paste(exposure, "_m"))
dat_exposure <- bind_rows(dat_exposure_w, dat_exposure_m)

################################################
##OUTCOME=T2D, DIAGRAM
################################################
## WOMEN, 30,053 T2D cases; 434,336 controls
diagram_w <- vroom("G:/Lessard/T2D/Mahajan.NatGenet2018b.T2D.FEMALE.European.txt") %>% 
  mutate(var=paste(Chr, Pos, sep=":")) %>% 
  filter(var %in% dat_exposure_w$var) %>%
  distinct(., var, .keep_all = TRUE) %>%
  select(-SNP)
glimpse(diagram_w)

T2D_w <- diagram_w %>%
  rename(effect_allele.outcome=EA, other_allele.outcome=NEA, beta.outcome=Beta, se.outcome=SE, eaf.outcome=EAF, pval.outcome=Pvalue) %>% 
  mutate(ncase.outcome=30053, ncontrol.outcome=434336, prevalence.outcome=0.0647, units.outcome="log odds",
  id.outcome=1, outcome="T2D, women")

dat_outcome_w <- merge(x = T2D_w, y = dat_exposure_w[,c("var","SNP")], by = "var", all.x = TRUE)
dat_outcome_w <- distinct(dat_outcome_w, SNP, .keep_all = TRUE) #very important!

harm_dat_w <- harmonise_data(dat_exposure_w, dat_outcome_w) %>% 
  filter(mr_keep=="TRUE")

## MEN, 41,846 T2D cases; 383,767 controls
diagram_m <- vroom("G:/Lessard/T2D/Mahajan.NatGenet2018b.T2D.MALE.European.txt") %>% 
  mutate(var=paste(Chr, Pos, sep=":")) %>% 
  filter(var %in% dat_exposure_m$var) %>%
  distinct(., var, .keep_all = TRUE) %>%
  select(-SNP)
glimpse(diagram_m)

T2D_m <- diagram_m  %>% 
  rename(effect_allele.outcome=EA, other_allele.outcome=NEA, beta.outcome=Beta, se.outcome=SE, eaf.outcome=EAF, pval.outcome=Pvalue) %>% 
  mutate(ncase.outcome=41846, ncontrol.outcome=383767, prevalence.outcome=0.0983, units.outcome="log odds",
  id.outcome=2, outcome="T2D, men") 

dat_outcome_m <- merge(x = T2D_m, y = dat_exposure_m[,c("var","SNP")], by = "var", all.x = TRUE)
dat_outcome_m <- distinct(dat_outcome_m, SNP, .keep_all = TRUE) #very important!

harm_dat_m <- harmonise_data(dat_exposure_m, dat_outcome_m) %>% 
  filter(mr_keep=="TRUE")

harm_dat <- rbind(harm_dat_w, harm_dat_m)

#################################################
#MR analyses
#################################################
res <- mr(harm_dat, method_list=c( "mr_ivw", "mr_weighted_median", "mr_egger_regression"))


ivw <- subset_on_method(res,
                        single_snp_method = "Wald ratio",
                        multi_snp_method = "Inverse variance weighted"
)

#HETREOGENEITY
het<-mr_heterogeneity(harm_dat)
het <- het %>%
  mutate(I2=100*(Q-Q_df)/Q) %>%
  mutate(I2=if_else(I2<0, 0, I2)) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/
het <- subset(het, method=="Inverse variance weighted") 

#PLEIOTROPY
#Egger intercept
plt<-mr_pleiotropy_test(harm_dat)
plt <- plt %>%
  rename(se.egger=se, pval.egger=pval) %>%
  mutate(method="MR Egger")

#Egger I2GX
i2gxf <- function(x, y) {
  list(id.exposure=y$id.exposure, id.outcome=y$id.outcome, I2GX=mr_egger(mr_input(bx=x$beta.exposure, bxse=x$se.exposure, by=x$beta.outcome, byse=x$se.outcome))@I.sq)
}

I2GX <- harm_dat %>% 
  group_by(id.exposure, id.outcome) %>%
  group_map(i2gxf)

df <- data.frame(matrix(unlist(I2GX), nrow=length(I2GX), byrow=T),stringsAsFactors=FALSE)
df$X3 <- as.numeric(df$X3)
mregger.i2gx <- df %>% 
  rename(id.exposure=X1, id.outcome=X2, I2gx=X3) %>%
  mutate(I2gx=100*I2gx) %>%
  mutate(method="MR Egger")

##MRPRESSO
dat <- subset(harm_dat, mr_keep)
d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)),
            select = c(exposure, outcome, id.exposure, id.outcome))
mrpresso <- list()
attributes(mrpresso)$id.exposure <- d$id.exposure
attributes(mrpresso)$id.outcome <- d$id.outcome
attributes(mrpresso)$exposure <- d$exposure
attributes(mrpresso)$outcome <- d$outcome
for (j in 1:nrow(d)) {
  x <- subset(dat, exposure == d$exposure[j] & outcome ==
                d$outcome[j])
  message(x$exposure[1], " - ", x$outcome[1])
  mrpresso[[j]] <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                                       BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                                       SdExposure = "se.exposure", OUTLIERtest = TRUE,
                                       DISTORTIONtest = TRUE, data = x, NbDistribution = nrow(x) * 20,
                                       SignifThreshold = 0.05)
}

df <- data.frame(matrix(unlist(mrpresso), nrow=length(mrpresso), byrow=T),stringsAsFactors=FALSE)
mr.presso <- df %>% 
  rename(b=X1, se=X2, t=X3, p=X4, no_b=X5, no_se=X6, no_t=X7, no_p=X8, global.p=X9, bias.p=X10, snp.total=X11, snp.outliers=X12, snp.list=X13)

#mr.presso <- cbind(d, mr.presso)
d <- ivw %>%
  mutate(snp.total=nsnp) %>%
  select(b, snp.total, exposure, outcome, id.exposure, id.outcome)
mr.presso <- merge(d, mr.presso, all.x = TRUE, all.y = TRUE) 

mr.presso <- merge(d, mr.presso, all.x = TRUE, all.y = TRUE) 
mr.presso <- mr.presso %>%
  mutate(pval=no_p) %>%
  mutate(b=no_b) %>%
  mutate(se=no_se) %>%
  mutate(method="MR-PRESSO") %>%
  select(id.exposure, id.outcome, exposure, outcome, method, b, se, pval, global.p, bias.p, snp.total, snp.outliers, snp.list ) 

mr.presso[,6:12] <- lapply(mr.presso[,6:12], as.numeric) 
mr.presso <- mr.presso %>%
  mutate(nsnp=snp.total-snp.outliers)

#########################################################
#COMBINE ALL MR RESULTS
#########################################################
mr_res <- merge(res, plt, all.x=TRUE, all.y=TRUE)
mr_res <- merge(mr_res, mregger.i2gx, all.x=TRUE, all.y=TRUE)
mr_res <- merge(mr_res, het, all.x=TRUE, all.y=TRUE)
mr_res <- bind_rows(mr_res, mr.presso) 
mr_res <-mr_res[order(mr_res$id.exposure, mr_res$id.outcome, mr_res$method),]

mr_res <- mr_res %>%
  mutate(lci=b - 1.96 * se) %>%
  mutate(uci=b + 1.96 * se) %>%
  mutate(lci=exp(lci)) %>%
  mutate(uci=exp(uci)) %>%
  mutate(b=exp(b)) %>%
  rename(est=b) %>%
  mutate(sex=if_else(id.exposure==1| id.exposure==3 |id.exposure==5, "women", "men")) 

names(mr_res)
mr_res <- mr_res[,c(25, 5, 4, 1,2, 3, 7, 23, 24, 9, 14, 16, 17, 10, 12, 13, 6)]

#no binary exposures, only outcomes

binary.outcome <- harm_dat %>%
  mutate(r.outcome=NA) %>%
  mutate(pval.outcome=NA)%>%
  mutate(samplesize.outcome=NA) 

binary.outcome$r.outcome <- get_r_from_lor(
  lor=binary.outcome$beta.outcome,
  af=binary.outcome$eaf.outcome,
  ncase=binary.outcome$ncase.outcome,
  ncontrol=binary.outcome$ncontrol.outcome,
  prevalence=binary.outcome$prevalence.outcome,
  model = "logit",
  correction = FALSE
)

harm_dat <- binary.outcome

#I have pre-calculated r.outcome for binary traits (from get_r_from_lor), and deleted pval for outcome when binary
#so r2 is calculated from r for binary traits, and from pval and sample size for continuous traits: https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/steiger.R#sym-directionality_test

out <- directionality_test(harm_dat) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res<- merge(out, mr_res, all.x=TRUE, all.y=TRUE)
names(mr_res)

mr_res <- mr_res[,c(8, 3,4, 9:20, 5:7)]

mr_res$method <- factor(mr_res$method, order = TRUE, 
                        levels = c("Inverse variance weighted", "MR-PRESSO", "Weighted median", "MR Egger"))


unique(mr_res$exposure)
mr_res$exposure <- factor(mr_res$exposure, order = TRUE, 
                          levels = c("ABSIUKB _w", "ABSIUKB _m", "WHIUKB _w", "WHIUKB _m", "HIUKB _w", "HIUKB _m"))

mr_res$sex <- factor(mr_res$sex, order = TRUE, 
                     levels = c("women", "men"))

mr_res <-mr_res[order(mr_res$exposure, mr_res$sex, mr_res$method),] 

#########################################################
#Steiger
#########################################################
harm_dat.steiger<-steiger_filtering(harm_dat) %>% 
  filter(steiger_dir=="TRUE") 

res.steiger <- mr(harm_dat.steiger, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

#Steiger
ivw.steiger <- subset_on_method(
  res.steiger,
  single_snp_method = "Wald ratio",
  multi_snp_method = "Inverse variance weighted"
)

#Heterogeneity
het.steiger<-mr_heterogeneity(harm_dat.steiger)
het.steiger <- het.steiger %>%
  mutate(I2=100*(Q-Q_df)/Q) %>%
  mutate(I2=if_else(I2<0, 0, I2))
het.steiger <- subset(het.steiger, method=="Inverse variance weighted")

#Pleiotropy
plt.steiger<-mr_pleiotropy_test(harm_dat.steiger)
plt.steiger <- plt.steiger %>%
  rename(se.egger=se, pval.egger=pval) %>%
  mutate(method="MR Egger")

#Egger I2GX, Steiger
I2GX <- harm_dat.steiger %>%
  group_by(id.exposure, id.outcome) %>%
  group_map(i2gxf)

df <- data.frame(matrix(unlist(I2GX), nrow=length(I2GX), byrow=T),stringsAsFactors=FALSE)
df$X3 <- as.numeric(df$X3)
mreggeri2gx.steiger <- df %>%
  rename(id.exposure=X1, id.outcome=X2, I2gx=X3) %>%
  mutate(I2gx=100*I2gx) %>%
  mutate(method="MR Egger")

###############################################
#MR-PRESSO, Steiger
###############################################
dat <- subset(harm_dat.steiger, mr_keep)
d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)),
            select = c(exposure, outcome, id.exposure, id.outcome))

mrpresso.steiger <- list()
attributes(mrpresso.steiger)$id.exposure <- d$id.exposure
attributes(mrpresso.steiger)$id.outcome <- d$id.outcome
attributes(mrpresso.steiger)$exposure <- d$exposure
attributes(mrpresso.steiger)$outcome <- d$outcome
for (j in 1:nrow(d)) {
  x <- subset(dat, exposure == d$exposure[j] & outcome ==
                d$outcome[j])
  message(x$exposure[1], " - ", x$outcome[1])
  mrpresso.steiger[[j]] <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                                               BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
                                               SdExposure = "se.exposure", OUTLIERtest = TRUE,
                                               DISTORTIONtest = TRUE, data = x, NbDistribution = nrow(x) * 20,
                                               SignifThreshold = 0.05)
}

df <- data.frame(matrix(unlist(mrpresso.steiger), nrow=length(mrpresso.steiger), byrow=T),stringsAsFactors=FALSE)
mrpresso.steiger <- df %>% 
  rename(b=X1, se=X2, t=X3, p=X4, no_b=X5, no_se=X6, no_t=X7, no_p=X8, global.p=X9, bias.p=X10, snp.total=X11, snp.outliers=X12, snp.list=X13)

d <- ivw.steiger %>% 
  mutate(snp.total=nsnp) %>%
  select(b, snp.total, exposure, outcome, id.exposure, id.outcome) 

mrpresso.steiger  <- merge(d, mrpresso.steiger, all.x = TRUE, all.y = TRUE) 
mrpresso.steiger  <- mrpresso.steiger  %>%
  mutate(pval=no_p) %>%
  mutate(b=no_b) %>%
  mutate(se=no_se) %>%
  mutate(method="MR-PRESSO") %>%
  select(id.exposure, id.outcome, exposure, outcome, method, b, se, pval, global.p, bias.p, snp.total, snp.outliers, snp.list ) 

mrpresso.steiger[,6:12] <- lapply(mrpresso.steiger[,6:12], as.numeric) 
mrpresso.steiger <- mrpresso.steiger %>%
  mutate(nsnp=snp.total-snp.outliers)

###############################################
#COMBINE ALL MR RESULTS, Steiger
###############################################
mr_res.steiger <- merge(res.steiger, plt.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, mreggeri2gx.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, het.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- bind_rows(mr_res.steiger, mrpresso.steiger) 
mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$id.exposure, mr_res.steiger$id.outcome, mr_res.steiger$method),]

mr_res.steiger <- mr_res.steiger %>%
  mutate(lci=b - 1.96 * se) %>%
  mutate(uci=b + 1.96 * se) %>%
  mutate(lci=exp(lci)) %>%
  mutate(uci=exp(uci)) %>%
  mutate(b=exp(b)) %>%
  rename(est=b) %>%
  mutate(sex=if_else(id.exposure==1| id.exposure==3 |id.exposure==5, "women", "men")) 

names(mr_res.steiger)
mr_res.steiger <- mr_res.steiger[,c(25, 5, 4, 1,2, 3, 7, 23, 24, 9, 14, 16, 17, 10, 12, 13, 6)]

#I have pre-calculated r.outcome for binary traits (from get_r_from_lor), and deleted pval for outcome when binary
#so r2 is calculated from r for binary traits, and from pval and sample size for continuous traits: https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/steiger.R#sym-directionality_test

out.steiger <- directionality_test(harm_dat.steiger) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res.steiger<- merge(out.steiger, mr_res.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- mr_res.steiger[,c(8, 3,4, 9:20, 5:7)]

mr_res.steiger$method <- factor(mr_res.steiger$method, order = TRUE, 
                                levels = c("Inverse variance weighted", "MR-PRESSO", "Weighted median", "MR Egger"))

mr_res.steiger$exposure <- factor(mr_res.steiger$exposure, order = TRUE, 
                                  levels = c("ABSIUKB _w", "ABSIUKB _m", "WHIUKB _w", "WHIUKB _m", "HIUKB _w", "HIUKB _m"))
mr_res.steiger$sex <- factor(mr_res.steiger$sex, order = TRUE, 
                     levels = c("women", "men"))
mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$exposure, mr_res.steiger$sex, mr_res.steiger$method),] 

#############################################
#SUPPLEMENTARY SUMMARY
#############################################
OUT <- createWorkbook() # Create a blank workbook

# Add some sheets to the workbook
addWorksheet(OUT, "ABSI_T2D")
addWorksheet(OUT, "ABSI_T2D_sf")

# Write the data to the sheets
writeData(OUT, sheet = "ABSI_T2D", x = mr_res)
writeData(OUT, sheet = "ABSI_T2D_sf", x = mr_res.steiger)

# Export the file
saveWorkbook(OUT, "Supplementary/ABSI_T2D.xlsx", overwrite = TRUE)

#################################################
ivw_1 <- ivw
save(ivw_1, file = "ivw_1.RData")
#################################################