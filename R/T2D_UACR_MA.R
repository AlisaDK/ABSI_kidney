library(openxlsx)
library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(vroom)
library(TwoSampleMR)
library(MendelianRandomization)
library(forestplot)
library(ggplot2)
library(patchwork)

#################################################
#EXPOSURE: T2D
#Genetic instruments: T2D EA
#Mahajan Tables S2 (DIAGRAM)
#################################################
Mahajan2018T2D <- read_excel("Mahajan2018T2D.xlsx", 
                             sheet = "ST 2") 
DM <- Mahajan2018T2D %>%
  rename(Pos='Position (Build 37 bp)') %>%
  mutate(var=paste(Chromosome, Pos, sep=":")) %>%
  rename(SNP='Index variant') %>% 
  rename(pval.exposure='p-value') %>% 
  rename(cases.exposure=Cases) 
DM$pval.exposure <- as.numeric(DM$pval.exposure) 
DM <- DM %>%
  filter(pval.exposure<5e-8) %>% # 298 SNPs
  select(SNP, var, pval.exposure)  %>%
  clump_data(., clump_r2 = 0.01) # 246 SNPs


diagram <- vroom("Mahajan.NatGenet2018b.T2D.European.txt") %>% 
  mutate(var=paste(Chr, Pos, sep=":")) %>% 
  distinct(., var, .keep_all = TRUE) %>%
  select(-SNP)

glimpse(diagram)

T2D <- diagram %>% 
  rename(effect_allele.exposure=EA, other_allele.exposure=NEA, beta.exposure=Beta, se.exposure=SE, eaf.exposure=EAF, pval.exposure=Pvalue) %>% 
  mutate(prevalence.exposure=0.0825, ncase.exposure=74124, ncontrol.exposure=824006, samplesize.exposure=898130, units.exposure="log odds",
  id.exposure="T2D", exposure="T2D") %>% 
  filter(var %in% DM$var) 

dat_exposure <- merge(x = T2D, y = DM[,c("var","SNP")], by = "var", all.x = TRUE) #missing SNP column

#################################################
#OUTCOME:UACR2019_EA
#################################################
dat <- vroom("outcomesEA/UACR2019_EA.csv") %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)
  
UACR_EA <- data.frame(SNP=dat$SNP, beta.outcome= dat$Effect, se.outcome= dat$StdErr, id.outcome=7, outcome="UACR_EA",
                      units.outcome="SD", pval.outcome=dat$`P-value`, samplesize.outcome=dat$n_total_sum, 
                      effect_allele.outcome= toupper(dat$Allele1), other_allele.outcome= toupper(dat$Allele2),  eaf.outcome= dat$Freq1)

#################################################
#OUTCOME: MA transethnic, mostly EA
#################################################
dat <- vroom("outcomes/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz") %>%
  filter(!is.na(RSID)) %>%
  rename(SNP=RSID) %>%
  filter(SNP %in% dat_exposure$SNP)
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
  mutate(samplesize.outcome=54000) %>%
  mutate(prevalence.outcome=0.10)

dat_outcome<- bind_rows(UACR_EA, MA_EA, MA_ALL)

############# HARMONISE DATA #########################
harm_dat <- harmonise_data(dat_exposure, dat_outcome) %>% 
  filter(mr_keep=="TRUE") 

######################################################
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

############################
##MR-PRESSO
#############################
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
  mutate(lci=if_else(outcome=="MA_EA" | outcome=="MA_ALL", exp(lci), lci)) %>%
  mutate(uci=if_else(outcome=="MA_EA" | outcome=="MA_ALL", exp(uci), uci)) %>%
  mutate(b=if_else(outcome=="MA_EA" | outcome=="MA_ALL", exp(b), b)) %>%
  rename(est=b) 

names(mr_res)
mr_res <- mr_res[,c(5, 4, 1,2, 3, 7, 23,24, 9, 14, 16, 17, 10, 12, 13, 6)]

#some binary outcomes (MA_EA and MA_ALL), and binary exposure (T2D)
binary.outcome <- harm_dat %>%
  filter(units.outcome=="log odds") %>%
  mutate(r.outcome=NA)%>%
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
binary.outcome <- binary.outcome %>%
  select(SNP, outcome, exposure, r.outcome)

binary.exposure <- harm_dat %>%
  filter(units.exposure=="log odds") %>%
  mutate(r.exposure=NA)%>%
  mutate(pval.exposure=NA)%>%
  mutate(samplesize.exposure=NA) 

binary.exposure$r.exposure <- get_r_from_lor(
  lor=binary.exposure$beta.exposure,
  af=binary.exposure$eaf.exposure,
  ncase=binary.exposure$ncase.exposure,
  ncontrol=binary.exposure$ncontrol.exposure,
  prevalence=binary.exposure$prevalence.exposure,
  model = "logit",
  correction = FALSE
)


harm_dat_new <- binary.exposure 
harm_dat_new_new <- merge(x=harm_dat_new, y=binary.outcome, by=c("SNP", "outcome", "exposure"), all.x=TRUE)

harm_dat <- harm_dat_new_new %>%
  mutate(pval.outcome = if_else(units.outcome=="log odds", as.numeric(NA), pval.outcome)) %>%
  mutate(samplesize.outcome = if_else(units.outcome=="log odds", as.numeric(NA), samplesize.outcome)) 

#I have pre-calculated r.exposure and r.outcome for binary traits (from get_r_from_lor), and deleted pval for outcome when binary
#so r2 is calculated from r for binary traits, and from pval and sample size for continuous traits: https://rdrr.io/github/MRCIEU/TwoSampleMR/src/R/steiger.R#sym-directionality_test

out <- directionality_test(harm_dat) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res<- merge(out, mr_res, all.x=TRUE, all.y=TRUE)
names(mr_res)

mr_res <- mr_res[,c(3,4, 8:19, 5:7)]

mr_res$method <- factor(mr_res$method, order = TRUE, 
                        levels = c("Inverse variance weighted", "MR-PRESSO", "Weighted median", "MR Egger"))

unique(mr_res$outcome)
mr_res$outcome <- factor(mr_res$outcome, order = TRUE, 
                          levels = c("UACR_EA", "MA_EA", "MA_ALL"))


mr_res <-mr_res[order(mr_res$outcome, mr_res$method),] 

#########################################################
#MR Steiger
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
df$X3 <- as.numeric(df$X3)
mreggeri2gx.steiger <- df %>%
  rename(id.exposure=X1, id.outcome=X2, I2gx=X3) %>%
  mutate(I2gx=100*I2gx) %>%
  mutate(method="MR Egger")

#####################################
#MR-PRESSO Steiger
###########################################
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

#####################################
#Merge all MR results, Steiger
###########################################
mr_res.steiger <- merge(res.steiger, plt.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, mreggeri2gx.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- merge(mr_res.steiger, het.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- bind_rows(mr_res.steiger, mrpresso.steiger) 
mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$id.exposure, mr_res.steiger$id.outcome, mr_res.steiger$method),]

mr_res.steiger <- mr_res.steiger %>%
  mutate(lci=b - 1.96 * se) %>%
  mutate(uci=b + 1.96 * se) %>%
  mutate(lci=if_else(outcome=="MA_EA" | outcome=="MA_ALL", exp(lci), lci)) %>%
  mutate(uci=if_else(outcome=="MA_EA" | outcome=="MA_ALL", exp(uci), uci)) %>%
  mutate(b=if_else(outcome=="MA_EA" | outcome=="MA_ALL", exp(b), b)) %>%
  rename(est=b)   

names(mr_res.steiger)
mr_res.steiger <- mr_res.steiger[,c(5, 4, 1,2, 3, 7, 23, 24, 9, 14, 16, 17, 10, 12, 13, 6)]

#I have already pre-calculated r.exposure and r.outcome for harm_dat, and as harm_dat.steiger is based on calculations after this calculation, ready to ho
out.steiger <- directionality_test(harm_dat.steiger) %>%
  select(id.exposure, id.outcome, exposure, outcome, snp_r2.exposure, snp_r2.outcome, correct_causal_direction) #, steiger_pval  

mr_res.steiger<- merge(out.steiger, mr_res.steiger, all.x=TRUE, all.y=TRUE)
mr_res.steiger <- mr_res.steiger[,c(3,4, 8:19, 5:7)]

mr_res.steiger$method <- factor(mr_res.steiger$method, order = TRUE, 
                                levels = c("Inverse variance weighted", "MR-PRESSO", "Weighted median", "MR Egger"))

mr_res.steiger$outcome <- factor(mr_res.steiger$outcome, order = TRUE, 
                                  levels = c("UACR_EA", "MA_EA", "MA_ALL"))
mr_res.steiger <-mr_res.steiger[order(mr_res.steiger$outcome, mr_res.steiger$method),] 

#############################################
#SUPPLEMENTARY SUMMARY
#############################################
OUT <- createWorkbook() # Create a blank workbook

# Add some sheets to the workbook
addWorksheet(OUT, "T2D_UACR")
addWorksheet(OUT, "T2D_UACR_sf")

# Write the data to the sheets
writeData(OUT, sheet = "T2D_UACR", x = mr_res)
writeData(OUT, sheet = "T2D_UACR_sf", x = mr_res.steiger)

# Export the file
saveWorkbook(OUT, "Supplementary/T2D_UACR.xlsx", overwrite = TRUE)
#################################################

#################################################
#Fig4: ABSI ->T2D, T2D -> UACR/MA, 
#################################################
load("G:/BMI_eGFR/R BMI eGFR/ivw_1.RData")

part_two <- ivw %>%
  mutate(id.exposure=as.numeric(7)) %>%
  mutate(sex="Both sexes")

part_one <- ivw_1 %>% 
  mutate(sex=if_else(id.exposure==1| id.exposure==3 |id.exposure==5, "Women", "Men")) 

fig.3 <- bind_rows(part_one, part_two) 

fig.3 <-fig.3[order(fig.3$outcome),] %>%
  mutate(beta=sprintf("%5.3f", b)) %>%
  mutate(lower=b-1.96*se) %>%
  mutate(upper=b+1.96*se) %>%
  mutate(lci=sprintf("%5.3f", lower)) %>%
  mutate(uci=sprintf("%5.3f", upper)) %>%
  mutate(est.size=paste(beta, "(",lci,",",uci,")"))

#align with table (number of rows to account for line jumps)
id.exposure <- c(NA, NA, 1, 3, 5, NA, 2,4,6, NA, 7,7,7)
id.outcome <- c(NA, NA, 1,1,1, NA, 2,2,2, NA, 7,12,11)
id <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
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
  c("Outcome", "T2D", NA, NA, NA, 
    "T2D", NA, NA, NA,
    "Albuminuria", "UACR, European", "MA, European", "MA, Transethnic"),
  c("Exposure", "Women", "ABSI (201 SNPs)","WHI (333 SNPs)","HI (177 SNPs)",
    "Men", "ABSI (53 SNPs)","WHI (87)","HI (65 SNPs)",
    "Both sexes", "T2D (234 SNPs)","T2D (158 SNPs)","T2D (233 SNPs)"),
  c(fig.3$est.size))
text.fig3[1, 3] <- "\u03b2 (95% CI)" # replace value in matrix: row, column, new value

#creating plot
pdf(file = 'newFig4.pdf', onefile=F) 
forestplot(text.fig3,
           fn.ci_norm = fpDrawCircleCI,
           txt_gp=fpTxtGp(ticks=gpar(cex=0.8),
                          xlab=gpar(cex=0.8),
                          label=gpar(cex=0.8)),
           title="Two-step MR: 1) allometric body shape indices -> T2D and 2) T2D -> albuminuria", 
           graph.pos = 3,
           mean=cbind(fig.3$mean),
           upper =cbind(fig.3$upper),
           lower= cbind(fig.3$lower),
           lwd.ci = 1,
           vertices=TRUE,
           boxsize = 0.175, 
           clip=c(-0.8, 0.8), 
           is.summary=c(TRUE,TRUE,rep(FALSE,3),TRUE,rep(FALSE,3),TRUE,rep(FALSE,3)),
           col = fpColors(box = "black",
                          line = "black"),
           xticks = c(-0.8, -0.4, 0, 0.4, 0.8),
           xlab=expression(paste(beta, " (95% CI)"))) 
dev.off() 