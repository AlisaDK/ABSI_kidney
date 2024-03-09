#################################################
#MR ANALYSES
#################################################
#devtools::install_github("MRCIEU/TwoSampleMR")
#devtools::install_github("thomasp85/patchwork")
# library(devtools)
# install_github("mnk/MR-PRESSO")
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(readr)
library(TwoSampleMR)
library(MendelianRandomization)
library(ggplot2)
library(forestplot)
library(ggplotify)
library(patchwork)

load("harm_ABSI_eGFR_EA.RData")
load("harm_ABSI_eGFR_sf.RData")

names(harm_dat)
names(harm_dat.steiger)

res_single <- mr_singlesnp(harm_dat) %>%
  select(outcome, samplesize, SNP, b, se, p) %>%
  filter(outcome!="urate")
res_single <- distinct(res_single) #, SNP, .keep_all = TRUE)

#F>29
F.strength <- harm_dat 
F.strength <- F.strength %>% 
  mutate(F.stat=(beta.exposure*beta.exposure)/(se.exposure*se.exposure)) %>%
  select(SNP, beta.exposure, se.exposure, F.stat)
  
#################################################
#IVW, WM and MR Egger
res <- mr(harm_dat, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

#Perform default MR analyses
#res <-mr_wrapper(dat, parameters = default_parameters())
#################################################
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



#save(mrpresso, file = "mrpresso_ABSI.RData")


load("G:/BMI_eGFR/R BMI eGFR/mrpresso_ABSI.RData")

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
  mutate(lci=if_else(outcome=="Rapid3"| outcome=="CKDi25"| outcome=="MA_ALL"| outcome=="MA_EA"| outcome=="CKD_EA", exp(lci), lci)) %>%
  mutate(uci=if_else(outcome=="Rapid3"| outcome=="CKDi25"| outcome=="MA_ALL"| outcome=="MA_EA"| outcome=="CKD_EA", exp(uci), uci)) %>%
  mutate(b=if_else(outcome=="Rapid3"| outcome=="CKDi25"| outcome=="MA_ALL"| outcome=="MA_EA"| outcome=="CKD_EA", exp(b), b)) %>%
  rename(est=b) %>%
  mutate(sex=if_else(id.exposure==1| id.exposure==3 |id.exposure==5, "women", "men")) 

names(mr_res)
mr_res <- mr_res[,c(25, 5, 4, 1,2, 3, 7, 23, 24,9, 14, 16, 17, 10, 12, 13, 6)]

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

binary.outcome <- binary.outcome %>%
  select(SNP, outcome, exposure, r.outcome)

harm_dat_new <- merge(x=harm_dat, y=binary.outcome, by=c("SNP", "outcome", "exposure"), all.x=TRUE)

harm_dat <- harm_dat_new %>%
  mutate(pval.outcome = if_else(units.outcome=="log odds", as.numeric(NA), pval.outcome)) %>%
  mutate(samplesize.outcome = if_else(units.outcome=="log odds", as.numeric(NA), samplesize.outcome)) 

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
unique(mr_res$outcome)
mr_res$outcome <- factor(mr_res$outcome, order = TRUE, 
                         levels = c("UACR_EA", "UACR_DM_EA", "UACR_noDM_EA", "MA_EA", 
                                    "eGFRcrea", "eGFR, DM", "eGFR, no DM",  "CKD_EA", 
                                    "BUN_EA", 
                                    "UACR_DM_ALL", "MA_ALL",  
                                    "eGFRcys",  "eGFR decline, all", "eGFR decline, DM", "Rapid3", "CKDi25", 
                                    "BUN_ALL"))
                         
mr_res$sex <- factor(mr_res$sex, order = TRUE, 
                         levels = c("women", "men"))
                         
mr_res <-mr_res[order(mr_res$outcome, mr_res$exposure, mr_res$sex, mr_res$method),] 
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
# library(devtools)
# install_github("mnk/MR-PRESSO")

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

#save(mrpresso.steiger, file = "mrpresso_ABSI_s.RData")


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
  mutate(lci=if_else(outcome=="Rapid3"| outcome=="CKDi25"| outcome=="MA_ALL"| outcome=="MA_EA"| outcome=="CKD_EA", exp(lci), lci)) %>%
  mutate(uci=if_else(outcome=="Rapid3"| outcome=="CKDi25"| outcome=="MA_ALL"| outcome=="MA_EA"| outcome=="CKD_EA", exp(uci), uci)) %>%
  mutate(b=if_else(outcome=="Rapid3"| outcome=="CKDi25"| outcome=="MA_ALL"| outcome=="MA_EA"| outcome=="CKD_EA", exp(b), b)) %>%
  rename(est=b) %>%
  mutate(sex=if_else(id.exposure==1| id.exposure==3 |id.exposure==5, "women", "men")) 


names(mr_res.steiger)
mr_res.steiger <- mr_res.steiger[,c(25, 5, 4, 1,2, 3, 7, 23, 24,9, 14, 16, 17, 10, 12, 13, 6)]
#I have pre-calculated r.outcome for binary traits (from get_r_from_lor) already!
#so r2 is calculated from r for binary traits, and from pval and sample size for continuous traits: https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/steiger.R#sym-directionality_test

out.steiger <- directionality_test(harm_dat.steiger) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res.steiger<- merge(out.steiger, mr_res.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- mr_res.steiger[,c(8, 3,4, 9:20, 5:7)]

mr_res.steiger$method <- factor(mr_res.steiger$method, order = TRUE, 
                        levels = c("Inverse variance weighted","MR-PRESSO",  "Weighted median", "MR Egger"))

mr_res.steiger$exposure <- factor(mr_res.steiger$exposure, order = TRUE, 
                                  levels = c("ABSIUKB _w", "ABSIUKB _m", "WHIUKB _w", "WHIUKB _m", "HIUKB _w", "HIUKB _m"))

mr_res.steiger$outcome <- factor(mr_res.steiger$outcome, order = TRUE, 
                                 levels = c("UACR_EA", "UACR_DM_EA", "UACR_noDM_EA", "MA_EA", 
                                            "eGFRcrea", "eGFR, DM", "eGFR, no DM",  "CKD_EA", 
                                            "BUN_EA", 
                                           "UACR_DM_ALL", "MA_ALL",  
                                            "eGFRcys",  "eGFR decline, all", "eGFR decline, DM", "Rapid3", "CKDi25", 
                                           "BUN_ALL"))

mr_res.steiger$sex <- factor(mr_res.steiger$sex, order = TRUE, 
                     levels = c("women", "men"))

mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$outcome, mr_res.steiger$exposure, mr_res.steiger$sex, mr_res.steiger$method),] 
mr_res.steiger <- mr_res.steiger %>%
  filter(!is.na(outcome))
#############################################
#SUPPLEMENTARY SUMMARY
#############################################
OUT <- createWorkbook() # Create a blank workbook

# Add some sheets to the workbook
addWorksheet(OUT, "ABSI_eGFR_EA")
addWorksheet(OUT, "ABSI_eGFR_EA_sf")
addWorksheet(OUT, "ABSI_single_SNPs")

# Write the data to the sheets
writeData(OUT, sheet = "ABSI_eGFR_EA", x = mr_res)
writeData(OUT, sheet = "ABSI_eGFR_EA_sf", x = mr_res.steiger)
writeData(OUT, sheet = "ABSI_single_SNPs", x = res_single)

# Export the file
saveWorkbook(OUT, "Supplementary/ABSI_eGFR.xlsx", overwrite = TRUE)
#############################################



###################################
#Fig 1: UACRC & MA, EA
#######################################
fig.1 <- ivw %>%
  filter(id.outcome==7 | id.outcome==9 | id.outcome==10 | id.outcome==12) %>%
  filter(id.exposure==1 | id.exposure==3 | id.exposure==5) #women

fig.1 <-fig.1[order(fig.1$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA,NA,NA, NA, NA, 1, 3, 5, NA, 1, 3, 5, NA, 1, 3, 5, NA, NA, 1, 3, 5) 
id.outcome <-  c(NA,NA, NA,NA, NA, 7, 7, 7, NA, 9, 9, 9, NA, 10, 10, 10, NA, NA, 12, 12, 12) 
id <-   c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
panel.fig1 <- data.frame(id.exposure, id.outcome, id)


fig.1 <-merge(panel.fig1, fig.1, all.x =  TRUE, all.y = TRUE)
fig.1 <-fig.1[order(fig.1$id),] %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b) #%>% 

#generating table for forestplot
text.fig1<-cbind(
  c("A. Women, European", NA, "Outcome", "UACR", "     Overall", NA, NA, NA, 
    "     Diabetes", NA, NA, NA,
    "     No diabetes", NA, NA, NA,
    "Microalbuminuria", "     Overall", NA, NA, NA),
  c(NA, NA, "Exposure",  NA, NA, "ABSI (201 SNPs)", "WHI (333 SNPs)", "HI (178 SNPs)",
    NA, "ABSI (147 SNPs)", "WHI (227 SNPs)", "HI (124 SNPs)",
    NA, "ABSI (148 SNPs)", "WHI (227 SNPs)", "HI (124 SNPs)",
    NA, NA, "ABSI (147 SNPs)", "WHI (228 SNPs)", "HI (124 SNPs)"),
  c(fig.1$est.size))
text.fig1[3, 3] <- "\u03b2 (95% CI)"

  
#creating plot
fig1a<- forestplot(text.fig1,
                   fn.ci_norm = fpDrawCircleCI,
                   txt_gp=fpTxtGp(ticks=gpar(cex=0.75),
                                  xlab=gpar(cex=0.75),
                                  label=gpar(cex=0.75)),
                   #title="Women, European", 
                   graph.pos = 3,
                   mean=cbind(fig.1$mean),
                   upper =cbind(fig.1$upper),
                   lower= cbind(fig.1$lower),
                   lwd.ci = 1,
                   vertices=TRUE,
                   boxsize = 0.25,
                   clip=c(-0.3, 0.3), 
                   is.summary=c(rep(TRUE,4), rep(FALSE,12),TRUE,rep(FALSE,12)),
                   col = fpColors(box = "black",
                                  line = "black"),
                   xticks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
                   xlab=expression(paste(beta, " (95% CI)")))

###################################
#Fig 1 B: UACRC & MA, EA MEN!
#######################################
fig.1 <- ivw %>%
  filter(id.outcome==7 | id.outcome==9 | id.outcome==10 | id.outcome==12) %>%
  filter(id.exposure==2 | id.exposure==4 | id.exposure==6) #Men

fig.1 <-fig.1[order(fig.1$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA,NA, NA, NA, 2, 4, 6, NA, 2, 4, 6, NA, 2, 4, 6, NA, NA, 2, 4, 6) 
id.outcome <-  c(NA, NA,NA, NA, NA, 7, 7, 7, NA, 9, 9, 9, NA, 10, 10, 10, NA, NA, 12, 12, 12) 
id <-   c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20,21)
panel.fig1 <- data.frame(id.exposure, id.outcome, id)


fig.1 <-merge(panel.fig1, fig.1, all.x =  TRUE, all.y = TRUE)
fig.1 <-fig.1[order(fig.1$id),] %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b) 

#generating table for forestplot
text.fig1<-cbind(
  c("B. Men, European", NA, "Outcome", "UACR", "     Overall", NA, NA, NA, 
    "     Diabetes", NA, NA, NA,
    "     No diabetes", NA, NA, NA,
    "Microalbuminuria", "     Overall", NA, NA, NA),
  c(NA, NA, "Exposure",  NA, NA, "ABSI (53 SNPs)", "WHI (88 SNPs)", "HI (66 SNPs)",
    NA, "ABSI (43 SNPs)", "WHI (71 SNPs)", "HI (51 SNPs)",
    NA, "ABSI (43 SNPs)", "WHI (71 SNPs)", "HI (52 SNPs)",
    NA, NA, "ABSI (43 SNPs)", "WHI (72 SNPs)", "HI (51 SNPs)"),
  c(fig.1$est.size))
text.fig1[3, 3] <- "\u03b2 (95% CI)"

#creating plot
fig1b <- forestplot(text.fig1,
                    fn.ci_norm = fpDrawCircleCI,
                    txt_gp=fpTxtGp(ticks=gpar(cex=0.75),
                                   xlab=gpar(cex=0.75),
                                   label=gpar(cex=0.75)),
                    #title="Men, European", 
                    graph.pos = 3,
                    mean=cbind(fig.1$mean),
                    upper =cbind(fig.1$upper),
                    lower= cbind(fig.1$lower),
                    lwd.ci = 1,
                    vertices=TRUE,
                    boxsize = 0.25,
                    clip=c(-0.3, 0.3), 
                    is.summary=c(rep(TRUE,4), rep(FALSE,12),TRUE,rep(FALSE,12)),
                    col = fpColors(box = "black",
                                   line = "black"),
                    xticks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
                    xlab=expression(paste(beta, " (95% CI)")))

# COMBINING FIGS 1A AND 1B
p1 <- grid2grob(print(fig1a))
p2 <- grid2grob(print(fig1b))
p_both <- wrap_elements(p1) | wrap_elements(p2)
p_both
#save as png or eps

###############################################################################
#Fig 2: UACRC_DM & MA, trans-ethnic
###############################################################################
fig.2 <-  ivw %>%
filter(id.outcome==8 | id.outcome==11) %>%
  filter(id.exposure==1 | id.exposure==3 | id.exposure==5) #women

fig.2 <-fig.2[order(fig.2$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA,NA, NA, 1, 3, 5, NA, 1, 3, 5)
id.outcome <- c(NA, NA,NA, NA, 8,8,8, NA, 11,11,11 )
id <- c(1, 2, 3, 4, 5, 6, 7, 8, 9,10,11)
panel.fig2 <- data.frame(id.exposure, id.outcome, id)

fig.2 <-merge(panel.fig2, fig.2, all.x = TRUE, all.y = TRUE)
fig.2 <-fig.2[order(fig.2$id),]  %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b)

#expression(paste(beta, " (95% CI)")))
#expression('title'[2])
#expression("A"['+'])
#generating table for forestplot
text.fig2<-cbind(
  c("A. Women, Transethnic", NA,"Outcome", "UACR, diabetes", NA, NA, NA, "Microalbuminuria",
    NA, NA, NA),
  c(NA, NA, "Exposure",  NA, "ABSI (189 SNPs)","WHI (319 SNPs)","HI (170 SNPs)",
    NA, "ABSI (202 SNPs)","WHI (338 SNPs)","HI (178 SNPs)"),
  c(fig.2$est.size))
text.fig2[3, 3] <- "\u03b2 (95% CI)" # replace value in matrix: row, column, new value

# #creating plot
fig2a <- forestplot(text.fig2,
               fn.ci_norm = fpDrawCircleCI,
               txt_gp=fpTxtGp(ticks=gpar(cex=0.8),
                              xlab=gpar(cex=0.8),
                              label=gpar(cex=0.8)),
               #"Women, Transethnic",
               graph.pos = 3,
               mean=cbind(fig.2$mean),
               upper =cbind(fig.2$upper),
               lower= cbind(fig.2$lower),
               lwd.ci = 1,
               vertices=TRUE,
               boxsize = 0.125,
               clip=c(-0.3, 0.3),
               is.summary=c(rep(TRUE,4),rep(FALSE,3),TRUE,rep(FALSE,3)),
               col = fpColors(box = "black",
                              line = "black"),
               xticks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
               xlab=expression(paste(beta, " (95% CI)")))

###############################################################################
#Fig 2B: UACRC_DM & MA, trans-ethnic MEN!
###############################################################################
fig.2 <-  ivw %>%
  filter(id.outcome==8 | id.outcome==11) %>%
  filter(id.exposure==2 | id.exposure==4 | id.exposure==6) #Men

fig.2 <-fig.2[order(fig.2$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA,NA, NA, 2,4,6, NA, 2,4,6)
id.outcome <- c(NA, NA, NA, NA,8,8,8, NA, 11,11,11 )
id <- c(1, 2, 3, 4, 5, 6, 7, 8, 9,10,11)
panel.fig2 <- data.frame(id.exposure, id.outcome, id)

fig.2 <-merge(panel.fig2, fig.2, all.x = TRUE, all.y = TRUE)
fig.2 <-fig.2[order(fig.2$id),]  %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b)

#expression(paste(beta, " (95% CI)")))
#expression('title'[2])
#expression("A"['+'])
#generating table for forestplot
text.fig2<-cbind(
  c("B. Men, Transethnic", NA, "Outcome", "UACR, diabetes", NA, NA, NA, "Microalbuminuria",
    NA, NA, NA),
  c(NA, NA,"Exposure",  NA, "ABSI (52 SNPs)","WHI (91 SNPs)","HI (68 SNPs)",
    NA, "ABSI (54 SNPs)","WHI (92 SNPs)","HI (70 SNPs)"),
  c(fig.2$est.size))
text.fig2[3, 3] <- "\u03b2 (95% CI)" # replace value in matrix: row, column, new value

# #creating plot
fig2b <- forestplot(text.fig2,
           fn.ci_norm = fpDrawCircleCI,
           txt_gp=fpTxtGp(ticks=gpar(cex=0.8),
                          xlab=gpar(cex=0.8),
                          label=gpar(cex=0.8)),
           #title="Men, Transethnic",
           graph.pos = 3,
           mean=cbind(fig.2$mean),
           upper =cbind(fig.2$upper),
           lower= cbind(fig.2$lower),
           lwd.ci = 1,
           vertices=TRUE,
           boxsize = 0.125,
           clip=c(-0.3, 0.3),
           is.summary=c(rep(TRUE,4),rep(FALSE,3),TRUE,rep(FALSE,3)),
           col = fpColors(box = "black",
                          line = "black"),
           xticks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
           xlab=expression(paste(beta, " (95% CI)")))

# COMBINING FIGS 2A AND 2B
p1 <- grid2grob(print(fig2a))
p2 <- grid2grob(print(fig2b))
p_both <- wrap_elements(p1) | wrap_elements(p2)
p_both
#save as png or eps

###############################################################################
#FigS1A: eGFRcrea, BUN & CKD, EA
###############################################################################
fig.3 <- ivw %>%
  filter(id.outcome==2 | id.outcome==6| id.outcome==13) %>%
  filter(id.exposure==1 | id.exposure==3 | id.exposure==5) #women

fig.3 <-fig.3[order(fig.3$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA,NA, NA, 1, 3,5, NA, 1, 3,5, NA, 1, 3,5)
id.outcome <- c(NA, NA,NA, NA, 2,2,2, NA, 13,13,13, NA, 6,6,6)
id <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14,15)
panel.fig3 <- data.frame(id.exposure, id.outcome, id)

fig.3 <-merge(panel.fig3, fig.3, all.x = TRUE, all.y = TRUE)
fig.3 <-fig.3[order(fig.3$id),]  %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b) 

#expression(paste(beta, " (95% CI)"))) 
#expression('title'[2])
#expression("A"['+'])
#generating table for forestplot
text.fig3<-cbind(
  c("A. Women, European", NA,"Outcome", "eGFRcrea", NA, NA, NA, 
    "CKD", NA, NA, NA,
    "BUN", NA, NA, NA),
  c(NA, NA, "Exposure",  NA, "ABSI (206 SNPs)","WHI (340 SNPs)","HI (182 SNPs)",
    NA, "ABSI (203 SNPs)","WHI (337 SNPs)","HI (182 SNPs)",
    NA, "ABSI (201 SNPs)","WHI (331 SNPs)","HI (177 SNPs)"),
  c(fig.3$est.size))
text.fig3[3, 3] <- "\u03b2 (95% CI)" # replace value in matrix: row, column, new value

#creating plot
figs1a <- forestplot(text.fig3,
           fn.ci_norm = fpDrawCircleCI,
           txt_gp=fpTxtGp(ticks=gpar(cex=0.8),
                          xlab=gpar(cex=0.8),
                          label=gpar(cex=0.8)),
           #title="Women, European", 
           graph.pos = 3,
           mean=cbind(fig.3$mean),
           upper =cbind(fig.3$upper),
           lower= cbind(fig.3$lower),
           lwd.ci = 1,
           vertices=TRUE,
           boxsize = 0.175, 
           clip=c(-0.2, 0.1), 
           is.summary=c(rep(TRUE,4),rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
           col = fpColors(box = "black",
                          line = "black"),
           xticks = c(-0.2, -0.1, 0, 0.1),
           xlab=expression(paste(beta, " (95% CI)"))) 

###############################################################################
#FigS1B: eGFRcrea, BUN & CKD, EA MEN!
###############################################################################
fig.3 <- ivw %>%
  filter(id.outcome==2 | id.outcome==6| id.outcome==13) %>%
  filter(id.exposure==2 | id.exposure==4 | id.exposure==6) #Men

fig.3 <-fig.3[order(fig.3$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA,NA, NA, 2,4,6, NA, 2,4,6, NA, 2,4,6)
id.outcome <- c(NA, NA, NA, NA,2,2,2, NA,13,13,13, NA, 6,6,6)
id <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14,15)
panel.fig3 <- data.frame(id.exposure, id.outcome, id)

fig.3 <-merge(panel.fig3, fig.3, all.x = TRUE, all.y = TRUE)
fig.3 <-fig.3[order(fig.3$id),]  %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b) 

#expression(paste(beta, " (95% CI)"))) 
#expression('title'[2])
#expression("A"['+'])
#generating table for forestplot
text.fig3<-cbind(
  c("B. Men, European", NA, "Outcome", "eGFRcrea", NA, NA, NA, 
    "CKD", NA, NA, NA,
    "BUN", NA, NA, NA),
  c(NA, NA,"Exposure",  NA, "ABSI (55 SNPs)","WHI (92 SNPs)","HI (71 SNPs)",
    NA, "ABSI (54 SNPs)","WHI (91 SNPs)","HI (71 SNPs)",
    NA, "ABSI (53 SNPs)","WHI (87 SNPs)","HI (65 SNPs)"),
  c(fig.3$est.size))
text.fig3[3, 3] <- "\u03b2 (95% CI)" # replace value in matrix: row, column, new value

#creating plot
figs1b <- forestplot(text.fig3,
           fn.ci_norm = fpDrawCircleCI,
           txt_gp=fpTxtGp(ticks=gpar(cex=0.8),
                          xlab=gpar(cex=0.8),
                          label=gpar(cex=0.8)),
           #title="Men, European", 
           graph.pos = 3,
           mean=cbind(fig.3$mean),
           upper =cbind(fig.3$upper),
           lower= cbind(fig.3$lower),
           lwd.ci = 1,
           vertices=TRUE,
           boxsize = 0.175, 
           clip=c(-0.2, 0.1), 
           is.summary=c(rep(TRUE,4),rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
           col = fpColors(box = "black",
                          line = "black"),
           xticks = c(-0.2, -0.1, 0, 0.1),
           xlab=expression(paste(beta, " (95% CI)"))) 

# COMBINING FIGS S1A AND S1B
p1 <- grid2grob(print(figs1a))
p2 <- grid2grob(print(figs1b))
p_both <- wrap_elements(p1) | wrap_elements(p2)
p_both
#save as png or eps

###################################
#FigS2A: eGFRcys, Rapid3, CKDi25 & BUN, trans-ethnic
#######################################
fig.4 <- ivw %>%
  filter(id.outcome==1 | id.outcome==3 | id.outcome==4 | id.outcome==5) %>%
  filter(id.exposure==1 | id.exposure==3 | id.exposure==5) #women

fig.4 <-fig.4[order(fig.4$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA,NA, NA, 1, 3, 5, NA, 1, 3, 5, NA, 1, 3, 5, NA, 1, 3, 5) 
id.outcome <-  c(NA, NA,NA, NA, 1,1,1, NA, 3,3,3, NA, 4,4,4, NA, 5,5,5) 
id <-   c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,18,19)
panel.fig4 <- data.frame(id.exposure, id.outcome, id)


fig.4 <-merge(panel.fig4, fig.4, all.x =  TRUE, all.y = TRUE)
fig.4 <-fig.4[order(fig.4$id),] %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b) #%>% 

#generating table for forestplot
text.fig4<-cbind(
  c("A. Women, Transethnic", NA, "Outcome", "eGFRcys", NA, NA, NA, 
    "Rapid3", NA, NA, NA,
    " CKDi25", NA, NA, NA,
    "BUN", NA, NA, NA),
  c(NA, NA, "Exposure",  NA, "ABSI (197 SNPs)", "WHI (332 SNPs)", "HI (179 SNPs)",
    NA, "ABSI (200 SNPs)", "WHI (338 SNPs)", "HI (179 SNPs)",
    NA, "ABSI (207 SNPs)", "WHI (344 SNPs)", "HI (183 SNPs)",
    NA, "ABSI (202 SNPs)", "WHI (333 SNPs)", "HI (177 SNPs)"),
  c(fig.4$est.size))
text.fig4[3, 3] <- "\u03b2 (95% CI)"

#creating plot
figs2a <- forestplot(text.fig4,
           fn.ci_norm = fpDrawCircleCI,
           txt_gp=fpTxtGp(ticks=gpar(cex=0.75),
                          xlab=gpar(cex=0.75),
                          label=gpar(cex=0.75)),
           #title="Women, Transethnic", 
           graph.pos = 3,
           mean=cbind(fig.4$mean),
           upper =cbind(fig.4$upper),
           lower= cbind(fig.4$lower),
           lwd.ci = 1,
           vertices=TRUE,
           boxsize = 0.25,
           clip=c(-0.2, 0.2), 
           is.summary=c(rep(TRUE,4), rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
           col = fpColors(box = "black",
                          line = "black"),
           xticks = c(-0.2, -0.1, 0, 0.1, 0.2),
           xlab=expression(paste(beta, " (95% CI)"))) 

###################################
#FigS2B: eGFRcys, Rapid3, CKDi25 & BUN, trans-ethnic MEN!
#######################################
fig.4 <- ivw %>%
  filter(id.outcome==1 | id.outcome==3 | id.outcome==4 | id.outcome==5) %>%
  filter(id.exposure==2 | id.exposure==4 | id.exposure==6) #men

fig.4 <-fig.4[order(fig.4$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA, NA, NA, 2,4,6, NA, 2,4,6, NA, 2,4,6, NA, 2,4,6) 
id.outcome <-  c(NA, NA, NA, NA, 1,1,1, NA, 3,3,3, NA, 4,4,4, NA, 5,5,5) 
id <-   c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,18,19)
panel.fig4 <- data.frame(id.exposure, id.outcome, id)


fig.4 <-merge(panel.fig4, fig.4, all.x =  TRUE, all.y = TRUE)
fig.4 <-fig.4[order(fig.4$id),] %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b) 

#generating table for forestplot
text.fig4<-cbind(
  c("B. Men, Transethnic", NA, "Outcome", "eGFRcys", NA, NA, NA, 
    "Rapid3", NA, NA, NA,
    " CKDi25", NA, NA, NA,
    "BUN", NA, NA, NA),
  c(NA, NA,"Exposure",  NA, "ABSI (53 SNPs)", "WHI (91 SNPs)", "HI (70 SNPs)",
    NA, "ABSI (54 SNPs)", "WHI (93 SNPs)", "HI (70 SNPs)",
    NA, "ABSI (55 SNPs)", "WHI (94 SNPs)", "HI (71 SNPs)",
    NA, "ABSI (54 SNPs)", "WHI (91 SNPs)", "HI (70 SNPs)"),
  c(fig.4$est.size))
text.fig4[3, 3] <- "\u03b2 (95% CI)"

#creating plot
figs2b <- forestplot(text.fig4,
           fn.ci_norm = fpDrawCircleCI,
           txt_gp=fpTxtGp(ticks=gpar(cex=0.75),
                          xlab=gpar(cex=0.75),
                          label=gpar(cex=0.75)),
           #title="Men, Transethnic", 
           graph.pos = 3,
           mean=cbind(fig.4$mean),
           upper =cbind(fig.4$upper),
           lower= cbind(fig.4$lower),
           lwd.ci = 1,
           vertices=TRUE,
           boxsize = 0.25,
           clip=c(-0.2, 0.2), 
           is.summary=c(rep(TRUE,4), rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
           col = fpColors(box = "black",
                          line = "black"),
           xticks = c(-0.2, -0.1, 0, 0.1, 0.2),
           xlab=expression(paste(beta, " (95% CI)"))) 

# COMBINING FIGS S2A AND S2B
p1 <- grid2grob(print(figs2a))
p2 <- grid2grob(print(figs2b))
p_both <- wrap_elements(p1) | wrap_elements(p2)
p_both
#save as png or eps