library(readxl)
library(TwoSampleMR)
library(dplyr)
library(tidyverse)
library(readr)
library(MendelianRandomization)
library(openxlsx)

`%notin%` <- Negate(`%in%`)

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
                             units.exposure="SD", pval.exposure=ABSI_w$gwasP, samplesize.exposure=219872, var=ABSI_w$var, #gene.exposure="Gene", mr
                             effect_allele.exposure= ABSI_w$effect_allele, other_allele.exposure= ABSI_w$non_effect_allele,  eaf.exposure= ABSI_w$MAF)
ABSI_w <- dat_exposure_w %>% 
  clump_data(., clump_r2 = 0.01) 

# men
dat_exposure_m <- data.frame(SNP = ABSI_m$SNP, beta.exposure= ABSI_m$beta, se.exposure= ABSI_m$se, id.exposure=2, exposure=ABSI_m$Index, sex="Men",
                             units.exposure="SD", pval.exposure=ABSI_m$gwasP, samplesize.exposure=186825, var=ABSI_m$var, gene.exposure="NEAREST_GENES", #mr
                             effect_allele.exposure= ABSI_m$effect_allele, other_allele.exposure= ABSI_m$non_effect_allele,  eaf.exposure= ABSI_m$MAF)
ABSI_m <- dat_exposure_m %>%   
  clump_data(., clump_r2 = 0.01) 


dat_exposure_w <- data.frame(SNP = WHI_w$SNP, beta.exposure= WHI_w$beta, se.exposure= WHI_w$se, id.exposure=3, exposure=WHI_w$Index, sex="Women",
                             units.exposure="SD", pval.exposure=WHI_w$gwasP, samplesize.exposure=219872, var=WHI_w$var, #gene.exposure="Gene", mr
                             effect_allele.exposure= WHI_w$effect_allele, other_allele.exposure= WHI_w$non_effect_allele,  eaf.exposure= WHI_w$MAF)
WHI_w <- dat_exposure_w %>%   
  clump_data(., clump_r2 = 0.01) 

# men
dat_exposure_m <- data.frame(SNP = WHI_m$SNP, beta.exposure= WHI_m$beta, se.exposure= WHI_m$se, id.exposure=4, exposure=WHI_m$Index, sex="Men",
                             units.exposure="SD", pval.exposure=WHI_m$gwasP, samplesize.exposure=186825, var=WHI_m$var, gene.exposure="NEAREST_GENES", #mr
                             effect_allele.exposure= WHI_m$effect_allele, other_allele.exposure= WHI_m$non_effect_allele,  eaf.exposure= WHI_m$MAF)
WHI_m <- dat_exposure_m %>%   
  clump_data(., clump_r2 = 0.01) 


dat_exposure_w <- data.frame(SNP = HI_w$SNP, beta.exposure= HI_w$beta, se.exposure= HI_w$se, id.exposure=5, exposure=HI_w$Index, sex="Women",
                             units.exposure="SD", pval.exposure=HI_w$gwasP, samplesize.exposure=219872, var=HI_w$var, #gene.exposure="Gene", mr
                             effect_allele.exposure= HI_w$effect_allele, other_allele.exposure= HI_w$non_effect_allele,  eaf.exposure= HI_w$MAF)
HI_w <- dat_exposure_w %>%   
  clump_data(., clump_r2 = 0.01) 

# men
dat_exposure_m <- data.frame(SNP = HI_m$SNP, beta.exposure= HI_m$beta, se.exposure= HI_m$se, id.exposure=6, exposure=HI_m$Index, sex="Men",
                             units.exposure="SD", pval.exposure=HI_m$gwasP, samplesize.exposure=186825, var=HI_m$var, gene.exposure="NEAREST_GENES", #mr
                             effect_allele.exposure= HI_m$effect_allele, other_allele.exposure= HI_m$non_effect_allele,  eaf.exposure= HI_m$MAF)
HI_m <- dat_exposure_m %>%   
  clump_data(., clump_r2 = 0.01) 

#############################################
################################################
#DIAGRAM
diagram_w <- vroom("G:/Lessard/T2D/Mahajan.NatGenet2018b.T2D.FEMALE.European.txt") # https://magicinvestigators.org/downloads/
diagram_m <- vroom("G:/Lessard/T2D/Mahajan.NatGenet2018b.T2D.MALE.European.txt") # https://magicinvestigators.org/downloads/
glimpse(diagram_w)

##ABSI
T2D_w <- diagram_w %>%
  mutate(var=paste(Chr, Pos, sep=":")) %>% 
  filter(var %in% ABSI_w$var) %>%
  distinct(., var, .keep_all = TRUE) 

T2D_w <- T2D_w %>%
  rename(effect_allele.outcome=EA, other_allele.outcome=NEA, beta.outcome=Beta, se.outcome=SE, eaf.outcome=EAF, pval.outcome=Pvalue, samplesize.outcome=Neff) %>% 
  mutate(id.outcome="T2D", outcome="T2D") %>%
  select(-SNP)

#merging SNPs
ABSI_T2D <- merge(x = T2D_w, y = ABSI_w[,c("var","SNP")], by = "var", all.x = TRUE) %>%
  select(SNP, !contains("exposure"))


##WHI
T2D_w <- diagram_w %>%
  mutate(var=paste(Chr, Pos, sep=":")) %>% 
  filter(var %in% WHI_w$var) %>%
  distinct(., var, .keep_all = TRUE)

T2D_w <- T2D_w %>%
  rename(effect_allele.outcome=EA, other_allele.outcome=NEA, beta.outcome=Beta, se.outcome=SE, eaf.outcome=EAF, pval.outcome=Pvalue, samplesize.outcome=Neff) %>% 
  mutate(id.outcome="T2D", outcome="T2D") %>%
  select(-SNP)

#merging SNPs
WHI_T2D <- merge(x = T2D_w, y = WHI_w[,c("var","SNP")], by = "var", all.x = TRUE)%>%
  select(SNP, !contains("exposure"))


##HI
T2D_w <- diagram_w %>%
  mutate(var=paste(Chr, Pos, sep=":")) %>% 
  filter(var %in% HI_w$var) %>%
  distinct(., var, .keep_all = TRUE)

T2D_w <- T2D_w %>%
  rename(effect_allele.outcome=EA, other_allele.outcome=NEA, beta.outcome=Beta, se.outcome=SE, eaf.outcome=EAF, pval.outcome=Pvalue, samplesize.outcome=Neff) %>% 
  mutate(id.outcome="T2D", outcome="T2D") %>%
  select(-SNP)

#merging SNPs
HI_T2D <- merge(x = T2D_w, y = HI_w[,c("var","SNP")], by = "var", all.x = TRUE) %>%
  select(SNP, !contains("exposure"))

##############################################################################################
#ABSI SNPs by concordance with T2D in DIAGRAM, only for women, as there was no assoc in men 
##############################################################################################
ABSI.T2D <- harmonise_data(ABSI_w, ABSI_T2D) %>% 
  filter(mr_keep=="TRUE") %>% 
  mutate(concordant=if_else(beta.exposure<0 & beta.outcome<0 | beta.exposure>0 & beta.outcome>0, 2, 0)) %>%
  mutate(concordant=if_else(pval.outcome>=0.05, 1, concordant)) %>% 
  filter(concordant==2) %>%
  mutate(id.exposure=21) %>%
  mutate(exposure="ABSI_T2D") %>% 
  select(SNP, contains("exposure"))

##################################################
#WHI SNPs by concordance with T2D in DIAGRAM
##################################################
WHI.T2D <- harmonise_data(WHI_w, WHI_T2D) %>% 
  filter(mr_keep=="TRUE") %>% 
  mutate(concordant=if_else(beta.exposure<0 & beta.outcome<0 | beta.exposure>0 & beta.outcome>0, 2, 0)) %>%
  mutate(concordant=if_else(pval.outcome>=0.05, 1, concordant)) %>% 
  filter(concordant==2) %>%
  mutate(id.exposure=22) %>%
  mutate(exposure="WHI_T2D") %>% 
  select(SNP, contains("exposure"))

##################################################
#HI SNPs by concordance with T2D in DIAGRAM
##################################################
HI.T2D <- harmonise_data(HI_w, HI_T2D) %>% 
  filter(mr_keep=="TRUE") %>% 
  mutate(concordant=if_else(beta.exposure<0 & beta.outcome<0 | beta.exposure>0 & beta.outcome>0, 2, 0)) %>%
  mutate(concordant=if_else(pval.outcome>=0.05, 1, concordant)) %>% 
  filter(concordant==2) %>%
  mutate(id.exposure=23) %>%
  mutate(exposure="HI_T2D") %>% 
  select(SNP, contains("exposure"))

#################################################
#Combining EXPOSURES, EA
#################################################
#Generate sex-stratified exposures
dat_exposure_w <- bind_rows(ABSI_w, WHI_w, HI_w, ABSI.T2D, WHI.T2D, HI.T2D) 
dat_exposure_m <- bind_rows(ABSI_m, WHI_m, HI_m) 
dat_exposure <- bind_rows(dat_exposure_w, dat_exposure_m)

#################################################
#OUTCOME:microalbuminuria in women, Neale lab 
#################################################
microalb_women <- read.delim("G:/BMI_eGFR/Nealelab/30500_irnt.gwas.imputed_v3.female.tsv")

microalb_w <- microalb_women %>%
  separate(variant, c("chr", "pos", "other_allele.outcome", "effect_allele.outcome"), ":", remove=F)

names(microalb_w)
microalb_w <- microalb_w %>%
  mutate(var=paste(chr, pos, sep=":"),
         se.outcome=se,
         id.outcome=1,
         outcome="alburia_w",
         samplesize.outcome=n_complete_samples,
         pval.outcome=pval,
         eaf.outcome=ifelse(minor_allele==effect_allele.outcome, minor_AF, 1-minor_AF),
         units.outcome="SD",
         beta.outcome=beta)
names(microalb_w)

microalb_w <- distinct(microalb_w, var, .keep_all = TRUE) %>% 
  filter(var %in% dat_exposure_w$var)

microalb_w <- microalb_w %>%
  select(var, other_allele.outcome, effect_allele.outcome, eaf.outcome, beta.outcome, se.outcome,
         pval.outcome, units.outcome, samplesize.outcome, id.outcome, outcome)

#merging SNPs
microalb_women <- merge(x = microalb_w, y = dat_exposure_w[,c("var","SNP")], by = "var") 
microalb_women <- distinct(microalb_women) 

#################################################
#OUTCOME:microalbuminuria in men, Neale lab 
#################################################
microalb_men <- read.delim("G:/BMI_eGFR/Nealelab/30500_irnt.gwas.imputed_v3.male.tsv")
microalb_m <- microalb_men %>%
  separate(variant, c("chr", "pos", "other_allele.outcome", "effect_allele.outcome"), ":", remove=F)

names(microalb_m)
microalb_m <- microalb_m %>%
  mutate(var=paste(chr, pos, sep=":"),
         se.outcome=se,
         id.outcome=2,
         outcome="alburia_m",
         samplesize.outcome=n_complete_samples,
         pval.outcome=pval,
         eaf.outcome=ifelse(minor_allele==effect_allele.outcome, minor_AF, 1-minor_AF),
         units.outcome="SD",
         beta.outcome=beta)
names(microalb_m)

microalb_m <- distinct(microalb_m, var, .keep_all = TRUE) %>% 
  filter(var %in% dat_exposure_m$var)

microalb_m <- microalb_m %>%
  select(var, other_allele.outcome, effect_allele.outcome, eaf.outcome, beta.outcome, se.outcome,
         pval.outcome, units.outcome, samplesize.outcome, id.outcome, outcome)

#merging SNPs
microalb_men <- merge(x = microalb_m, y = dat_exposure_m[,c("var","SNP")], by = "var") 
microalb_men <- distinct(microalb_men) 

#################################################
#OUTCOME:Microalbumin in spot urine,N=108,706
#https://gwas.mrcieu.ac.uk/datasets/ieu-a-1100/
#Dataset: ukb-d-30500_irnt sex-stratified
#################################################
microalbumin <- extract_outcome_data(
  snps = dat_exposure$SNP,
  outcomes = 'ukb-d-30500_irnt',
  proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3
)
microalbumin <- microalbumin %>%
  mutate(id.outcome="3") %>%
  mutate(id.outcome=as.numeric(id.outcome)) %>%
  mutate(outcome="microalbumin")%>%
  mutate(units.outcome="SD")

#################################################
#Combining all outcomes
#################################################
dat_outcome <- bind_rows(microalb_women, microalb_men)

################################################# 
#HARMONISING DATA, both sexes
#################################################
harm_dat <- harmonise_data(dat_exposure, dat_outcome)
harm_dat <- harm_dat %>%
  filter(!(id.exposure==1 & id.outcome==2 |
             id.exposure==2 & id.outcome==1 |
             id.exposure==3 & id.outcome==2 |
             id.exposure==4 & id.outcome==1 |
             id.exposure==5 & id.outcome==2 |
             id.exposure==6 & id.outcome==1 |
             id.exposure==21 & id.outcome==2 |
             id.exposure==22 & id.outcome==2 |
             id.exposure==23 & id.outcome==2))
save(harm_dat, file = "harm_ABSI_Neal.RData")

#################################################
#Steiger filtering
#################################################
harm_dat.steiger<-steiger_filtering(harm_dat) %>%
  filter(steiger_dir=="TRUE")
save(harm_dat, file = "harm_ABSI_Neal_sf.RData")

load("G:/BMI_eGFR/R BMI eGFR/harm_ABSI_Neal.RData")
load("G:/BMI_eGFR/R BMI eGFR/harm_ABSI_Neal_sf.RData")

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
het <- subset(het, method=="Inverse variance weighted") %>%
  mutate

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
  mutate(I2gx=100*I2gx) %>%
  mutate(method="MR Egger")


# dat <- subset(harm_dat, mr_keep)
# d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)),
#             select = c(exposure, outcome, id.exposure, id.outcome))
# mrpresso <- list()
# attributes(mrpresso)$id.exposure <- d$id.exposure
# attributes(mrpresso)$id.outcome <- d$id.outcome
# attributes(mrpresso)$exposure <- d$exposure
# attributes(mrpresso)$outcome <- d$outcome
# for (j in 1:nrow(d)) {
#   x <- subset(dat, exposure == d$exposure[j] & outcome ==
#                 d$outcome[j])
#   message(x$exposure[1], " - ", x$outcome[1])
#   mrpresso[[j]] <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
#                                        BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
#                                        SdExposure = "se.exposure", OUTLIERtest = TRUE,
#                                        DISTORTIONtest = TRUE, data = x, NbDistribution = nrow(x) * 20,
#                                        SignifThreshold = 0.05)
# }
# 
# 
# 
# save(mrpresso, file = "mrpresso_ABSI_Neal.RData")

load("G:/BMI_eGFR/R BMI eGFR/mrpresso_ABSI_Neal.RData")
df <- data.frame(matrix(unlist(mrpresso), nrow=length(mrpresso), byrow=T),stringsAsFactors=FALSE)
mr.presso <- df %>% 
  rename(b=X1, se=X2, t=X3, p=X4, no_b=X5, no_se=X6, no_t=X7, no_p=X8, global.p=X9, bias.p=X10, snp.total=X11, snp.outliers=X12, snp.list=X13)

d <- ivw %>% 
  mutate(snp.total=nsnp) %>%
  select(b, snp.total, exposure, outcome, id.exposure, id.outcome) 

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
  rename(est=b) %>%
  mutate(sex=if_else(id.exposure==2| id.exposure==4 |id.exposure==6, "men", "women")) 

names(mr_res)
mr_res <- mr_res[,c(25, 5, 4, 1,2, 3, 7, 23, 24,9, 14, 16, 17, 10, 12, 13, 6)]

out <- directionality_test(harm_dat) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res<- merge(out, mr_res, all.x=TRUE, all.y=TRUE)
names(mr_res)

mr_res <- mr_res[,c(8, 3,4, 9:20, 5:7)]

mr_res$method <- factor(mr_res$method, order = TRUE, 
                        levels = c("Inverse variance weighted", "MR-PRESSO", "Weighted median", "MR Egger"))


unique(mr_res$exposure)
mr_res$exposure <- factor(mr_res$exposure, order = TRUE, 
                          levels = c("ABSIUKB", "WHIUKB", "HIUKB", "ABSI_T2D", "WHI_T2D","HI_T2D"))
unique(mr_res$outcome)
mr_res$outcome <- factor(mr_res$outcome, order = TRUE, 
                         levels = c("alburia_w", "alburia_m"))

mr_res$sex <- factor(mr_res$sex, order = TRUE, 
                     levels = c("women", "men"))

mr_res <-mr_res[order(mr_res$sex, mr_res$outcome, mr_res$exposure, mr_res$method),] 
mr_res <- mr_res %>%
  filter(!is.na(outcome))


#########################################################
#COMBINE ALL MR RESULTS, Steiger
#########################################################
res.steiger <- mr(harm_dat.steiger, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

#Steiger
ivw.steiger <- subset_on_method(
  res.steiger,
  single_snp_method = "Wald ratio",
  multi_snp_method = "Inverse variance weighted"
)
wm.steiger <- subset_on_method(
  res.steiger,
  single_snp_method = "Wald ratio",
  multi_snp_method = "Weighted median"
)
egger.steiger <- subset_on_method(
  res.steiger,
  single_snp_method = "Wald ratio",
  multi_snp_method = "MR Egger"
)

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
i2gxf <- function(x, y) {
  list(id.exposure=y$id.exposure, id.outcome=y$id.outcome, I2GX=mr_egger(mr_input(bx=x$beta.exposure, bxse=x$se.exposure, by=x$beta.outcome, byse=x$se.outcome))@I.sq)
}
I2GX <- harm_dat.steiger %>%
  group_by(id.exposure, id.outcome) %>%
  group_map(i2gxf)


df <- data.frame(matrix(unlist(I2GX), nrow=length(I2GX), byrow=T),stringsAsFactors=FALSE)
mreggeri2gx.steiger <- df %>%
  rename(id.exposure=X1, id.outcome=X2, I2gx=X3) %>%
  mutate(I2gx=100*I2gx) %>%
  mutate(method="MR Egger")

#MR-PRESSO
# dat <- subset(harm_dat.steiger, mr_keep)
# d <- subset(dat, !duplicated(paste(id.exposure, " - ", id.outcome)),
#             select = c(exposure, outcome, id.exposure, id.outcome))
# 
# mrpresso.steiger <- list()
# attributes(mrpresso.steiger)$id.exposure <- d$id.exposure
# attributes(mrpresso.steiger)$id.outcome <- d$id.outcome
# attributes(mrpresso.steiger)$exposure <- d$exposure
# attributes(mrpresso.steiger)$outcome <- d$outcome
# for (j in 1:nrow(d)) {
#  x <- subset(dat, exposure == d$exposure[j] & outcome ==
#                d$outcome[j])
#  message(x$exposure[1], " - ", x$outcome[1])
#  mrpresso.steiger[[j]] <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
#                                       BetaExposure = "beta.exposure", SdOutcome = "se.outcome",
#                                       SdExposure = "se.exposure", OUTLIERtest = TRUE,
#                                       DISTORTIONtest = TRUE, data = x, NbDistribution = nrow(x) * 20,
#                                       SignifThreshold = 0.05)
# }
# 
# save(mrpresso.steiger, file = "mrpresso_ABSI_s.RData")


load("G:/BMI_eGFR/R BMI eGFR/mrpresso_ABSI_s.RData")

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

############################
#COMBINE ALL RESULTS
#############################
mr_res.steiger <- merge(res.steiger, plt.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, mreggeri2gx.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, het.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- bind_rows(mr_res.steiger, mrpresso.steiger) 
mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$id.exposure, mr_res.steiger$id.outcome, mr_res.steiger$method),]

mr_res.steiger <- mr_res.steiger %>%
  mutate(lci=b - 1.96 * se) %>%
  mutate(uci=b + 1.96 * se) %>%
  rename(est=b) %>%
  mutate(sex=if_else(id.exposure==2| id.exposure==4 |id.exposure==6, "men", "women")) 

names(mr_res.steiger)
mr_res.steiger <- mr_res.steiger[,c(25, 5, 4, 1,2, 3, 7, 23, 24,9, 14, 16, 17, 10, 12, 13, 6)]

out.steiger <- directionality_test(harm_dat.steiger) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res.steiger<- merge(out.steiger, mr_res.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- mr_res.steiger[,c(8, 3,4, 9:20, 5:7)]

mr_res.steiger$method <- factor(mr_res.steiger$method, order = TRUE, 
                                levels = c("Inverse variance weighted","MR-PRESSO",  "Weighted median", "MR Egger"))

mr_res.steiger$exposure <- factor(mr_res.steiger$exposure, order = TRUE, 
                                  levels = c("ABSIUKB", "WHIUKB","HIUKB", "ABSI_T2D", "WHI_T2D", "HI_T2D"))

mr_res.steiger$outcome <- factor(mr_res.steiger$outcome, order = TRUE, 
                                 levels = c("alburia_w", "alburia_m"))

mr_res.steiger$sex <- factor(mr_res.steiger$sex, order = TRUE, 
                             levels = c("women", "men"))

mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$sex,mr_res.steiger$outcome, mr_res.steiger$exposure,  mr_res.steiger$method),] 
mr_res.steiger <- mr_res.steiger %>%
  filter(!is.na(outcome))

#############################################
#SUPPLEMENTARY SUMMARY
#############################################
OUT <- createWorkbook() # Create a blank workbook

# Add some sheets to the workbook
addWorksheet(OUT, "ABSI_Neale")
addWorksheet(OUT, "ABSI_Neale_sf")

# Write the data to the sheets
writeData(OUT, sheet = "ABSI_Neale", x = mr_res)
writeData(OUT, sheet = "ABSI_Neale_sf", x = mr_res.steiger)

# Export the file
saveWorkbook(OUT, "Supplementary/ABSI_Neale.xlsx", overwrite = TRUE)
#############################################