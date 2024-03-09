library(readr)
library(readxl)
library(openxlsx)
library(vroom)
library(dplyr)
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)

#################################################
#EXPOSURE: UACR EA
#################################################
dat <- read_excel("UACR_exposures.xlsx", 
                            sheet = "european") %>%
  mutate(var=paste(chr,`position (b37)`, sep=":"))

UACR <- data.frame(SNP =dat$SNP, beta.exposure= dat$Effect, se.exposure= dat$SE, id.exposure=1, exposure="UACR_EA",
                          units.exposure="SD", pval.exposure=dat$`P-value`, samplesize.exposure=dat$N, var=dat$var,
                          effect_allele.exposure= dat$EA, other_allele.exposure=dat$NEA,  eaf.exposure= dat$FreqEA) %>%
  clump_data(., clump_r2 = 0.01)


dat_exposure <- UACR   


################################################
##OUTCOME=T2D, DIAGRAM
################################################
diagram <- vroom("Mahajan.NatGenet2018b.T2D.European.txt") %>% 
  mutate(var=paste(Chr, Pos, sep=":"))   %>% 
  filter(var %in% UACR$var) %>% 
  distinct(., var, .keep_all = TRUE) %>%
  select(-SNP)

glimpse(diagram)

T2D <- diagram %>% 
  rename(effect_allele.outcome=EA, other_allele.outcome=NEA, beta.outcome=Beta, se.outcome=SE, eaf.outcome=EAF, pval.outcome=Pvalue) %>% 
  mutate(prevalence.outcome=0.0825, ncase.outcome=74124, ncontrol.outcome=824006, samplesize.outcome=898130, units.outcome="log odds",
         id.outcome="T2D", outcome="T2D") 

__________


##ABSI, WHI, HI
# T2D <- distinct(diagram, SNP, .keep_all = TRUE) %>% 
#   rename(var=SNP) %>% 
#   rename(effect_allele.outcome=EA) %>% 
#   rename(other_allele.outcome=NEA) %>% 
#   rename(beta.outcome=Beta) %>% 
#   rename(se.outcome=SE) %>% 
#   rename(eaf.outcome=EAF) %>% 
#   rename(pval.outcome=Pvalue) %>% 
#   mutate(ncase.outcome=74124) %>%
#   mutate(ncontrol.outcome=824006) %>%
#   mutate(prevalence.outcome=0.0825)%>%
#   mutate(units.outcome="log odds")%>%
#   mutate(id.outcome=1) %>% 
#   mutate(outcome="T2D") %>% 
#   filter(var %in% dat_exposure$var) 

dat_outcome <- merge(x = T2D, y = UACR[,c("var","SNP")], by = "var", all.x = TRUE)
glimpse(dat_outcome)
##################################################
#Harmonising
##################################################
harm_dat <- harmonise_data(UACR, dat_outcome) %>% 
  filter(mr_keep=="TRUE") 

##################################################
#MR analyses
##################################################
res <- mr(harm_dat, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

res_single <- mr_singlesnp(harm_dat) %>%
  select(outcome, SNP, b, se, p) 
res_single <- distinct(res_single) %>%
  mutate(ncase.outcome=74124) %>%
  mutate(ncontrol.outcome=824006) %>%
  mutate(effective.samplesize=231420)

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
mregger.i2gx <- df %>% 
  rename(id.exposure=X1, id.outcome=X2, I2gx=X3) %>%
  mutate(I2gx=as.numeric(I2gx)) %>%
  mutate(I2gx=100*I2gx) %>%
  mutate(method="MR Egger")

################################
#MR-PRESSO
##################################
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
  mutate(lci=exp(b - 1.96 * se)) %>%
  mutate(uci=exp(b + 1.96 * se)) %>%
  mutate(est=exp(b)) 

names(mr_res)
mr_res <- mr_res[,c(5, 4, 1,2, 3, 25, 23, 24, 9, 14, 16, 17, 10, 12, 13, 6)]

#MR Steiger, directionality test
#no binary exposures, only outcomes
binary.outcome <- harm_dat %>%
  filter(units.outcome=="log odds") %>%
  mutate(r.outcome=NA)

binary.outcome$r.outcome <- get_r_from_lor(
  lor=binary.outcome$beta.outcome,
  af=binary.outcome$eaf.outcome,
  ncase=binary.outcome$ncase.outcome,
  ncontrol=binary.outcome$ncontrol.outcome,
  prevalence=binary.outcome$prevalence.outcome,
  model = "logit",
  correction = FALSE
)

#I have pre-calculated r.outcome for binary traits (from get_r_from_lor), and deleted pval for outcome when binary
#so r2 is calculated from r for binary traits, and from pval and sample size for continuous traits: https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/steiger.R#sym-directionality_test
out <- directionality_test(binary.outcome) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res<- merge(out, mr_res, all.x=TRUE, all.y=TRUE)
names(mr_res)

mr_res <- mr_res[,c(3,4, 8:19, 5:7)]

mr_res$method <- factor(mr_res$method, order = TRUE, 
                        levels = c("Inverse variance weighted", "MR-PRESSO", "Weighted median", "MR Egger"))

unique(mr_res$exposure)
mr_res$exposure <- factor(mr_res$exposure, order = TRUE, 
                          levels = c("UACR_EA", "UACR_DM_EA"))

mr_res <-mr_res[order(mr_res$exposure, mr_res$method),] 
mr_res <- mr_res %>%
  filter(!is.na(outcome))

#########################################################
# After Steiger filtering
#########################################################
harm_dat.steiger<-steiger_filtering(binary.outcome) %>% 
  filter(steiger_dir=="TRUE")

res.steiger <- mr(harm_dat.steiger, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

#Steiger
ivw.steiger <- subset_on_method(
  res.steiger,
  single_snp_method = "Wald ratio",
  multi_snp_method = "Inverse variance weighted")

het.steiger<-mr_heterogeneity(harm_dat.steiger)
het.steiger <- het.steiger %>%
  mutate(I2=100*(Q-Q_df)/Q) %>%
  mutate(I2=if_else(I2<0, 0, I2))
het.steiger <- subset(het.steiger, method=="Inverse variance weighted")

plt.steiger<-mr_pleiotropy_test(harm_dat.steiger)
plt.steiger <- plt.steiger %>%
  rename(se.egger=se, pval.egger=pval) %>%
  mutate(method="MR Egger")

#Egger I2GX, Steiger
I2GX <- harm_dat.steiger %>%
  group_by(id.exposure, id.outcome) %>%
  group_map(i2gxf)

df <- data.frame(matrix(unlist(I2GX), nrow=length(I2GX), byrow=T),stringsAsFactors=FALSE)
mreggeri2gx.steiger <- df %>%
  rename(id.exposure=X1, id.outcome=X2, I2gx=X3) %>%
  mutate(I2gx=as.numeric(I2gx)) %>%
  mutate(I2gx=100*I2gx) %>%
  mutate(method="MR Egger")

######################################
#MR-PRESSO Steiger
######################################
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

######################################
#COMBINE ALL MR RESULTS Steiger
######################################
mr_res.steiger <- merge(res.steiger, plt.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, mreggeri2gx.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, het.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- bind_rows(mr_res.steiger, mrpresso.steiger) 
mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$id.exposure, mr_res.steiger$id.outcome, mr_res.steiger$method),]

mr_res.steiger <- mr_res.steiger %>%
  mutate(lci=exp(b - 1.96 * se)) %>%
  mutate(uci=exp(b + 1.96 * se)) %>%
  mutate(est=exp(b))

names(mr_res.steiger)
mr_res.steiger <- mr_res.steiger[,c(5, 4, 1,2, 3, 25, 23, 24, 9, 14, 16, 17, 10, 12, 13, 6)]

#I have pre-calculated r.outcome for binary traits (from get_r_from_lor) already!
#so r2 is calculated from r for binary traits, and from pval and sample size for continuous traits: https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/steiger.R#sym-directionality_test
out.steiger <- directionality_test(harm_dat.steiger) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res.steiger<- merge(out.steiger, mr_res.steiger, all.x=TRUE, all.y=TRUE)
names(mr_res.steiger)
mr_res.steiger <- mr_res.steiger[,c(3,4,8:19, 5:7)]

mr_res.steiger$method <- factor(mr_res.steiger$method, order = TRUE, 
                                levels = c("Inverse variance weighted", "MR-PRESSO", "Weighted median", "MR Egger"))

mr_res.steiger$exposure <- factor(mr_res.steiger$exposure, order = TRUE, 
                                  levels = c("UACR_EA", "UACR_DM_EA"))


mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$exposure, mr_res.steiger$method),] 
mr_res.steiger <- mr_res.steiger %>%
  filter(!is.na(outcome))

#############################################
#SUPPLEMENTARY SUMMARY
#############################################
OUT <- createWorkbook() # Create a blank workbook

# Add some sheets to the workbook
addWorksheet(OUT, "UACR_T2D")
addWorksheet(OUT, "UACR_T2D_sf")
addWorksheet(OUT, "single_SNPs")

# Write the data to the sheets
writeData(OUT, sheet = "UACR_T2D", x = mr_res)
writeData(OUT, sheet = "UACR_T2D_sf", x = mr_res.steiger)
writeData(OUT, sheet = "single_SNPs", x = res_single)

# Export the file
saveWorkbook(OUT, "Supplementary/UACR_T2D.xlsx", overwrite = TRUE)
#############################################