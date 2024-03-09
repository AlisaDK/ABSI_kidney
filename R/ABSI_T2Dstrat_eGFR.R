library(openxlsx)
library(readxl)
library(readr)
library(dplyr)
library(tidyverse)
library(TwoSampleMR)
library(MendelianRandomization)
library(forestplot)
library(ggplot2)
library(patchwork)

load("G:/BMI_eGFR/R BMI eGFR/harm_ABSI_T2D.RData") # saved in harm_ABSI_eGFR.R
load("G:/BMI_eGFR/R BMI eGFR/harm_ABSI_T2D_sf.RData")
harm_dat <- harm_dat_T2D

res <- mr(harm_dat, method_list=c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))

#################################################
ivw <- subset_on_method(res,
                        single_snp_method = "Wald ratio",
                        multi_snp_method = "Inverse variance weighted")

# HETREOGENEITY
het<-mr_heterogeneity(harm_dat)
het <- het %>%
  mutate(I2=100*(Q-Q_df)/Q) %>%
  mutate(I2=if_else(I2<0, 0, I2)) #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC192859/
het <- subset(het, method=="Inverse variance weighted")

# PLEIOTROPY: Egger intercept
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
  mutate(lci=if_else(outcome=="MA_ALL"| outcome=="MA_EA", exp(lci), lci)) %>%
  mutate(uci=if_else(outcome=="MA_ALL"| outcome=="MA_EA", exp(uci), uci)) %>%
  mutate(b=if_else(outcome=="MA_ALL"| outcome=="MA_EA", exp(b), b)) %>%
  rename(est=b) 
names(mr_res)
mr_res <- mr_res[,c(5, 4, 1,2, 3, 7, 23, 24,9, 14, 16, 17, 10, 12, 13, 6)]

#MR Steiger
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

mr_res <- mr_res[,c(8, 3,4, 9:19, 5:7)]

mr_res$method <- factor(mr_res$method, order = TRUE, 
                        levels = c("Inverse variance weighted", "MR-PRESSO", "Weighted median", "MR Egger"))

unique(mr_res$exposure)
mr_res$exposure <- factor(mr_res$exposure, order = TRUE, 
                          levels = c("ABSI_T2D", "WHI_T2D", "HI_T2D"))
unique(mr_res$outcome)
mr_res$outcome <- factor(mr_res$outcome, order = TRUE, 
                         levels = c("UACR_EA", "MA_EA", "MA_ALL"))

mr_res <-mr_res[order(mr_res$outcome, mr_res$exposure, mr_res$method),] 

#########################################################
# After Steiger filtering
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
mreggeri2gx.steiger <- df %>%
  rename(id.exposure=X1, id.outcome=X2, I2gx=X3) %>%
  mutate(I2gx=100*I2gx) %>%
  mutate(method="MR Egger")

##MRPRESSO STEIGER
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

##################################################
##MERGE ALL MR RESULTS STEIGER
#####################################################################
mr_res.steiger <- merge(res.steiger, plt.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, mreggeri2gx.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, het.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- bind_rows(mr_res.steiger, mrpresso.steiger) 
mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$id.exposure, mr_res.steiger$id.outcome, mr_res.steiger$method),]

mr_res.steiger <- mr_res.steiger %>%
  mutate(lci=b - 1.96 * se) %>%
  mutate(uci=b + 1.96 * se) %>%
  mutate(lci=if_else(outcome=="MA_ALL"| outcome=="MA_EA", exp(lci), lci)) %>%
  mutate(uci=if_else(outcome=="MA_ALL"| outcome=="MA_EA", exp(uci), uci)) %>%
  mutate(b=if_else(outcome=="MA_ALL"| outcome=="MA_EA", exp(b), b)) %>%
  rename(est=b) 

names(mr_res.steiger)
mr_res.steiger <- mr_res.steiger[,c(5, 4, 1,2, 3, 7, 23, 24,9, 14, 16, 17, 10, 12, 13, 6)]

#I have pre-calculated r.outcome for binary traits (from get_r_from_lor), and deleted pval for outcome when binary
#so r2 is calculated from r for binary traits, and from pval and sample size for continuous traits: https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/steiger.R#sym-directionality_test

out.steiger <- directionality_test(harm_dat.steiger) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res.steiger<- merge(out.steiger, mr_res.steiger, all.x=TRUE, all.y=TRUE)
names(mr_res.steiger)
mr_res.steiger <- mr_res.steiger[,c(3,4, 8:19, 5:7)]

mr_res.steiger$method <- factor(mr_res.steiger$method, order = TRUE, 
                                levels = c("Inverse variance weighted","MR-PRESSO",  "Weighted median", "MR Egger"))

mr_res.steiger$exposure <- factor(mr_res.steiger$exposure, order = TRUE, 
                                  levels = c("ABSI_T2D", "WHI_T2D", "HI_T2D"))

mr_res.steiger$outcome <- factor(mr_res.steiger$outcome, order = TRUE, 
                                 levels = c("UACR_EA", "MA_EA", "MA_ALL"))

mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$outcome, mr_res.steiger$exposure, mr_res.steiger$method),] 

#############################################
#SUPPLEMENTARY SUMMARY
#############################################
OUT <- createWorkbook() # Create a blank workbook

# Add some sheets to the workbook
addWorksheet(OUT, "ABSI_T2Dstrat_eGFR")
addWorksheet(OUT, "ABSI_T2Dstrat_eGFR_sf")

# Write the data to the sheets
writeData(OUT, sheet = "ABSI_T2Dstrat_eGFR", x = mr_res)
writeData(OUT, sheet = "ABSI_T2Dstrat_eGFR_sf", x = mr_res.steiger)

# Export the file
saveWorkbook(OUT, "Supplementary/ABSI_T2Dstrat_eGFR.xlsx", overwrite = TRUE)
###############################################################################

###############################################################################
#Fig 3: UACR EA, MA EA+trans-ethnic, for women only!
###############################################################################
fig.3 <- ivw 

fig.3 <-fig.3[order(fig.3$id.outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA,NA, NA, 21, 22,23, NA, 21, 22,23, NA, 21,22,23)
id.outcome <- c(NA, NA,NA, NA, 7,7,7, NA, 12,12,12, NA, 11,11,11)
id <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14,15)
panel.fig3 <- data.frame(id.exposure, id.outcome, id)

fig.3 <-merge(panel.fig3, fig.3, all.x = TRUE, all.y = TRUE)
fig.3 <-fig.3[order(fig.3$id),]  %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(mean=b) 

text.fig3<-cbind(
  c("Women", NA, "Outcome", "UACR, European", NA, NA, NA, 
    "MA, European", NA, NA, NA,
    "MA, Transethnic", NA, NA, NA),
  c(NA, NA, "Exposure",  NA, "ABSI, T2D (49 SNPs)","WHI, T2D (88 SNPs)","HI, T2D (9 SNPs)",
    NA, "ABSI, T2D (35 SNPs)","WHI, T2D (56 SNPs)","HI, T2D (8 SNPs)",
    NA, "ABSI, T2D (47 SNPs)","WHI, T2D (87 SNPs)","HI, T2D (9 SNPs)"),
  c(fig.3$est.size))
text.fig3[3, 3] <- "\u03b2 (95% CI)" # replace value in matrix: row, column, new value

#creating plot
library(forestplot)
pdf(file = 'Fig3new.pdf', onefile=F) 
forestplot(text.fig3,
           fn.ci_norm = fpDrawCircleCI,
           txt_gp=fpTxtGp(ticks=gpar(cex=0.8),
                          xlab=gpar(cex=0.8),
                          label=gpar(cex=0.8)),
           #title="Women", 
           graph.pos = 3,
           mean=cbind(fig.3$mean),
           upper =cbind(fig.3$upper),
           lower= cbind(fig.3$lower),
           lwd.ci = 1,
           vertices=TRUE,
           boxsize = 0.175, 
           clip=c(-0.25, 0.75), 
           is.summary=c(rep(TRUE,4),rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
           col = fpColors(box = "black",
                          line = "black"),
           xticks = c(-0.25, 0, 0.25, 0.5, 0.75),
           xlab=expression(paste(beta, " (95% CI)"))) 
dev.off() 