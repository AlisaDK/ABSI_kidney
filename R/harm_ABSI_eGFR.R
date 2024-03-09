library(readxl)
library(TwoSampleMR)
library(dplyr)
library(tidyverse)
library(readr)
library(vroom)
library("writexl")

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
dat_exposure <- bind_rows(ABSI_w, ABSI_m, WHI_w, WHI_m, HI_m, HI_w)
names(dat_exposure)
snp_list <- dat_exposure %>%
  select(sex, exposure, SNP, beta.exposure, se.exposure, effect_allele.exposure, other_allele.exposure, eaf.exposure)

write_xlsx(snp_list,"Supplementary\\snp_list.xlsx")

#Generate sex-stratified exposures
dat_exposure_w <- bind_rows(ABSI_w, WHI_w, HI_w) %>%
  mutate(exposure=paste(exposure, "_w"))
dat_exposure_m <- bind_rows(ABSI_m, WHI_m, HI_m) %>%
  mutate(exposure=paste(exposure, "_m"))
dat_exposure <- bind_rows(dat_exposure_w, dat_exposure_m)
#################################################
#OUTCOME:eGFRcys, new download from CKDGen 
#################################################
dat <- vroom("outcomes/metal_eGFRcys_meta1.TBL.map.annot.gc.gz")
glimpse(dat)
dat <- dat %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE)
  
eGFRcys <-   data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=1, outcome="eGFRcys",
                      units.outcome="SD", pval.outcome=dat$P.value, samplesize.outcome=dat$n, 
                      effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#################################################
#OUTCOME:eGFRcrea_EA, new download from CKDGen
#################################################
dat <- vroom("outcomesEA/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz")
glimpse(dat)
dat <-  dat %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE)
  
eGFRcrea_EA <-  data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=2, outcome="eGFRcrea",
                          units.outcome="SD", pval.outcome=dat$P.value, samplesize.outcome=dat$n, 
                          effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)


#################################################
#OUTCOME:rapid eGFR decline, mostly EA
#################################################
dat <- vroom("outcomes/CKDGen_rapid3_overall.txt.gz")
glimpse(dat)
dat <-  dat %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 

Rapid3<- data.frame(SNP = dat$SNP, beta.outcome= dat$OR, se.outcome= dat$StdErr, id.outcome=3, outcome="Rapid3",
                    units.outcome="log odds", pval.outcome=dat$'P-value', samplesize.outcome=dat$N_total_sum, 
                    ncase.outcome=34874, ncontrol.outcome=107090, prevalence.outcome=0.25,
                    effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)
Rapid3 <- Rapid3 %>% 
  mutate(beta.outcome=log(beta.outcome)) 


# CKDi25
dat <- vroom("outcomes/CKDGen_ckdi25_overall.txt.gz")
glimpse(dat)

dat <-  dat %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 

CKDi25<- data.frame(SNP = dat$SNP, beta.outcome= dat$OR, se.outcome= dat$StdErr, id.outcome=4, outcome="CKDi25",
                    units.outcome="log odds", pval.outcome=dat$'P-value', samplesize.outcome=dat$N_total_sum, 
                    ncase.outcome=19901, ncontrol.outcome=175244, prevalence.outcome=0.10,
                    effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

CKDi25 <- CKDi25 %>% 
  mutate(beta.outcome=log(beta.outcome)) 

#################################################
#OUTCOME:BUN, new download from CKDGen
#################################################
dat <- vroom("outcomes/metal_bun_meta_all1.TBL.map.annot.gc.gz")
glimpse(dat)
dat <-  dat %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 

BUN <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=5, outcome="BUN_ALL",
                      units.outcome="SD", pval.outcome=dat$P.value, samplesize.outcome=dat$n, 
                      effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#################################################
#OUTCOME:BUN_EA, Wuttke 2019
#################################################
dat <- vroom("outcomesEA/BUN_CKDGen_EA.txt")
glimpse(dat)
dat <-  dat %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 

BUN_EA <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=6, outcome="BUN_EA",
                     units.outcome="SD", pval.outcome=dat$'P-value', samplesize.outcome=dat$n_total_sum,
                     effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#################################################
#OUTCOME:UACR2019_EA
#################################################
dat <- vroom("outcomesEA/UACR2019_EA.csv")
problems(dat)  #313 SNPs without rs ID do not have this column
glimpse(dat)
dat <-  dat %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 
UACR_EA <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=7, outcome="UACR_EA",
                      units.outcome="SD", pval.outcome=dat$`P-value`, samplesize.outcome=dat$n_total_sum, 
                      effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#################################################
#OUTCOME:UACR diabetes only, transethnic, mostly EA
#################################################
dat <- vroom("outcomes/formatted_20170020-UACR_DM-All-nstud_18-SumMac_400.tbl.rsid.gz") #Warning: 277 parsing failures #These are lacking rsids
glimpse(dat)
dat <-  dat  %>%
  filter(!is.na(RSID)) %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 

UACR_DM_ALL <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=8, outcome="UACR_DM_ALL",
                          units.outcome="SD", pval.outcome=dat$`P-value`, samplesize.outcome=dat$n_total_sum, 
                          effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#################################################
# OUTCOME:UACR2015_DM_EA,N=5,825, PMID:	26631737
# https://gwas.mrcieu.ac.uk/datasets/ieu-a-1100/
#################################################
UACR_DM_EA <- extract_outcome_data(
  snps = dat_exposure$SNP,
  outcomes = 'ieu-a-1100',
  proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3
)
UACR_DM_EA <- UACR_DM_EA %>%
  mutate(id.outcome="9") %>%
  mutate(id.outcome=as.numeric(id.outcome)) %>%
  mutate(outcome="UACR_DM_EA")%>%
  mutate(units.outcome="SD")

#################################################
#OUTCOME:UACR2015_noDM_EA
#WITHOUT DM: https://gwas.mrcieu.ac.uk/datasets/ieu-a-1101/
#Both with and without DM: https://gwas.mrcieu.ac.uk/datasets/ieu-a-1107/
#################################################
UACR_noDM_EA <- extract_outcome_data(
  snps = dat_exposure$SNP,
  outcomes = 'ieu-a-1101',
  proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3
)

UACR_noDM_EA <- UACR_noDM_EA %>%
  mutate(id.outcome="10") %>%
  mutate(id.outcome=as.numeric(id.outcome)) %>%
  mutate(outcome="UACR_noDM_EA")%>%
  mutate(units.outcome="SD")

#################################################
#OUTCOME: MA transethnic, mostly EA
#################################################
dat <- vroom("outcomes/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz") #Warning: 292 parsing failures #These are lacking rsids

dat <- dat %>%
  filter(!is.na(RSID)) %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 
glimpse(dat)

MA_ALL <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=11, outcome="MA_ALL",
                     units.outcome="log odds", pval.outcome=dat$`P-value`, samplesize.outcome=dat$n_total_sum, 
                     ncase.outcome=51861, ncontrol.outcome=297093, prevalence.outcome=0.15,
                     effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#################################################
#OUTCOME: MA2015_EA
#################################################
MA_EA <- extract_outcome_data(
  snps = dat_exposure$SNP,
  outcomes = 'ieu-a-1097',
  proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3
)

MA_EA <- MA_EA %>%
  mutate(id.outcome="12") %>%
  mutate(id.outcome=as.numeric(id.outcome)) %>%
  mutate(outcome="MA_EA")%>%
  mutate(units.outcome="log odds") %>%
  mutate(ncase.outcome=5400) %>%
  mutate(ncontrol.outcome=48600) %>%
  mutate(prevalence.outcome=0.10)

#################################################
#OUTCOME:CKD_EA, Wuttke 2019
#################################################
dat <- vroom("outcomesEA/CKD_overall_EA_nstud23.txt")
glimpse(dat)
dat <-  dat  %>%
  filter(!is.na(RSID)) %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 
CKD_EA <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=13, outcome="CKD_EA",
                     units.outcome="log odds", pval.outcome=dat$`P-value`, samplesize.outcome=as.numeric(dat$n_total_sum), 
                     ncase.outcome=41395, ncontrol.outcome=439303, prevalence.outcome=0.086,
                     effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#################################################
#Added February 8, 2023
# Annual eGFR decline, Gorski et al., Kidney Int. 2022 
#################################################
# CKDGen_eGFR-decline_overall_adjDM.txt.gz, though equivalent to those not adjusted for BMI
dat<- vroom("outcomesALL/CKDGen_eGFR-decline_overall_adjDM.txt.gz")
glimpse(dat)
dat <-  dat  %>%
  filter(!is.na(RSID)) %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 
egfr_decline_all <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=14, outcome="eGFR decline, all",
                     units.outcome="SD", pval.outcome=dat$`P-value`, samplesize.outcome=dat$N_total_sum, 
                     effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

# These GWAS summary statistics for eGFR-decline are based on 37,375 individuals WITH diabetes at baseline adjusted for age and sex. 
dat <- vroom("outcomesALL/CKDGen_eGFR-decline_DM.txt.gz")
glimpse(dat)
dat <-  dat  %>%
  filter(!is.na(RSID)) %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 
egfr_decline_dm <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=15, outcome="eGFR decline, DM",
                               units.outcome="SD", pval.outcome=dat$`P-value`, samplesize.outcome=dat$N_total_sum, 
                               effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#######################################################################################
#Added February 8, 2023
# eGFR in DM, Winkler TW. Commun Biol 2022;5(1):580. doi:10.1038/s42003-022-03448-z
#######################################################################################
dat <- vroom("outcomesEA/meta_egfr_dmstrat_stage1plus2_EURonly.txt.gz")
glimpse(dat)
dat <-  dat  %>%
  filter(!is.na(rsid)) %>%
  rename(SNP=rsid) %>%
  filter(SNP %in% dat_exposure$SNP)%>% 
  distinct(., SNP, .keep_all = TRUE) 
egfr_dm <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect.dm, se.outcome= dat$StdErr.dm, id.outcome=16, outcome="eGFR, DM",
                              units.outcome="SD", pval.outcome=dat$P.value.dm, samplesize.outcome=dat$n.dm, 
                              effect_allele.outcome= dat$Allele1, other_allele.outcome= dat$Allele2,  eaf.outcome= dat$Freq1.dm)

egfr_nodm <- data.frame(SNP = dat$SNP, beta.outcome= dat$Effect.nodm, se.outcome= dat$StdErr.nodm, id.outcome=17, outcome="eGFR, no DM",
                      units.outcome="SD", pval.outcome=dat$P.value.nodm, samplesize.outcome=dat$n.nodm, 
                      effect_allele.outcome= dat$Allele1, other_allele.outcome= dat$Allele2,  eaf.outcome= dat$Freq1.nodm)

#################################################
#ALL OUTCOMES combined
#################################################
dat_outcome <- bind_rows(eGFRcys, eGFRcrea_EA, Rapid3, CKDi25, BUN, BUN_EA, 
                         UACR_EA, UACR_DM_ALL, UACR_DM_EA, UACR_noDM_EA,
                         MA_ALL, MA_EA, CKD_EA,
                         egfr_decline_all, egfr_decline_dm,
                         egfr_dm, egfr_nodm) 

################################################# 
#HARMONISING DATA
#################################################
harm_dat <- harmonise_data(dat_exposure, dat_outcome) %>% 
  filter(mr_keep=="TRUE") #Removing SNPs for being palindromic with intermediate allele frequencies:

#preliminary results
res <- mr(harm_dat, method_list=c( "mr_ivw", "mr_weighted_median", "mr_egger_regression"))

save(harm_dat, file = "harm_ABSI_eGFR.RData")

#################################################
#Steiger filtering
#################################################
harm_dat.steiger<-steiger_filtering(harm_dat) %>% 
  filter(steiger_dir=="TRUE")
save(harm_dat.steiger, file = "harm_ABSI_eGFR_sf.RData")
##########################
#T2D concordance
##############################
dat_exposure <- bind_rows(ABSI.T2D, WHI.T2D, HI.T2D) %>%
  select(SNP, contains("exposure"))

dat_outcome <- bind_rows(UACR_EA, MA_ALL, MA_EA)

harm_dat_T2D <- harmonise_data(dat_exposure, dat_outcome) %>% 
  filter(mr_keep=="TRUE") #Removing SNPs for being palindromic with intermediate allele frequencies:

#preliminary results
res <- mr(harm_dat_T2D, method_list=c( "mr_ivw", "mr_weighted_median", "mr_egger_regression"))

save(harm_dat_T2D, file = "harm_ABSI_T2D.RData")