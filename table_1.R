##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: table_1.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Descriptive characteristics of study population, stratified by race (Table 1)

# Statistical Analyses:
    # Descriptive statistics

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####

rm(list = ls())
main.dir = ""
setwd(main.dir)

# Load packages
library(Biobase)
library(tableone)
library(tidyverse)

# Load data
load(file.path(main.dir,"data/known.eset.final.RData"))

# Extract phenotype data
pdata <- pData(known.eset.final)

# Numbers in the text
# % non-Hispanic
tab1 <- table(pdata$nhs1_hisp)
tab1
round(prop.table(tab1),2)

# Age - median [IQR]
round(summary(pdata$knn_nhs1_agebld),1)

# Create data labels for race and for categorical variables with >2 categories

# Race
pdata <- mutate(pdata, Race = factor(race_bw),
                Race = fct_recode(Race,"Black" = "1","White" = "0"),
                Race = fct_relevel(Race,"Black","White"))
summary(pdata$Race)

# Non-hispanic indicator
pdata$nhs1_nonhisp = ifelse(pdata$nhs1_hisp==1, 0, 1)
pdata <- mutate(pdata, nhs1_nonhisp = factor(nhs1_nonhisp),
                nhs1_nonhisp = fct_recode(nhs1_nonhisp,"Non-Hispanic" = "1","Hispanic" = "0"))
summary(pdata$nhs1_nonhisp)

# Fasting
pdata <- mutate(pdata, knn_nhs1_fast = factor(knn_nhs1_fast),
                knn_nhs1_fast = fct_relevel(knn_nhs1_fast,
                                "0","1"))
summary(pdata$knn_nhs1_fast)

# Menopausal status in questionnaire cycle prior to blood collection
pdata <- mutate(pdata, knn_nhs1_menopmh_88_der = factor(knn_nhs1_menopmh_88_der),
                knn_nhs1_menopmh_88_der = fct_recode(knn_nhs1_menopmh_88_der,
                                              "Premenopausal" = "1",
                                              "Postmenopausal, no HRT" = "2",
                                              "Postmenopausal, HRT" = "3"),
                knn_nhs1_menopmh_88_der = fct_relevel(knn_nhs1_menopmh_88_der,
                                                   "Premenopausal",
                                                   "Postmenopausal, no HRT",
                                                   "Postmenopausal, HRT"))
summary(pdata$knn_nhs1_menopmh_88_der)

# Current smoker (binary for the table)
table(pdata$knn_nhs1_smk_status)
pdata <- mutate(pdata, knn_nhs1_smk_status = factor(knn_nhs1_smk_status),
                knn_nhs1_smk_status = fct_recode(knn_nhs1_smk_status,
                                              "Never" = "1",
                                              "Former" = "2",
                                              "Current" = "3"),
                knn_nhs1_smk_status = fct_relevel(knn_nhs1_smk_status,
                                  "Never",
                                  "Former",
                                  "Current"))
table(pdata$knn_nhs1_smk_status)

# Cigarettes smoked per day
pdata <- mutate(pdata, knn_nhs1_smk_per_day = factor(knn_nhs1_smk_per_day),
                knn_nhs1_smk_per_day = fct_recode(knn_nhs1_smk_per_day,
                                             "0 cigarettes/day" = "0",
                                             "1-14 cigarettes/day" = "1",
                                             "15+ cigarettes/day" = "2"))
summary(pdata$knn_nhs1_smk_per_day)

# History of hypertension (excluding treatment)
pdata <- mutate(pdata, knn_nhs1_hypthx_med = factor(knn_nhs1_hypthx_med),
                knn_nhs1_hypthx_med = fct_recode(knn_nhs1_hypthx_med,
                                             "No" = "0",
                                             "Yes" = "1",
                                             "Yes" = "2"),
                knn_nhs1_hypthx_med = fct_relevel(knn_nhs1_hypthx_med,
                                              "No",
                                              "Yes"))
summary(pdata$knn_nhs1_hypthx_med)

# History of high cholesterol and treatment
pdata <- mutate(pdata, knn_nhs1_cholhx_med = factor(knn_nhs1_cholhx_med),
                knn_nhs1_cholhx_med = fct_recode(knn_nhs1_cholhx_med,
                                             "No" = "0",
                                             "Yes" = "1",
                                             "Yes" = "2"),
                knn_nhs1_cholhx_med = fct_relevel(knn_nhs1_cholhx_med,
                                              "No",
                                              "Yes"))
summary(pdata$knn_nhs1_cholhx_med)

# History of diabetes and treatment
pdata <- mutate(pdata, knn_nhs1_dbhx_med = factor(knn_nhs1_dbhx_med),
                knn_nhs1_dbhx_med = fct_recode(knn_nhs1_dbhx_med,
                                             "No" = "0",
                                             "Yes" = "1",
                                             "Yes" = "2"),
                knn_nhs1_dbhx_med = fct_relevel(knn_nhs1_dbhx_med,
                                              "No",
                                              "Yes"))
summary(pdata$knn_nhs1_dbhx_med)

# Convert census tract % variables from decimal to % for table presentation
pdata$knn_nhs1_percent_college_cens <- pdata$knn_nhs1_percent_college_cens*100
pdata$knn_nhs1_percent_white_cens <- pdata$knn_nhs1_percent_white_cens*100

# Convert logged income to dollar scale
pdata$knn_nhs1_income_cens <- exp(pdata$knn_nhs1_ln_income_cens)
summary(pdata$knn_nhs1_income_cens)

# Region
summary(pdata$nhs1_region, useNA="always")

# Birth cohort
summary(pdata$nhs1_birth_cohort)

# Create list of variables for table 1
all.vars = c("nhs1_nonhisp","knn_nhs1_agebld","knn_nhs1_fast","knn_nhs1_menopmh_88_der",
             "knn_nhs1_hysthx","knn_nhs1_oophhx",
             "knn_nhs1_hypthx_med","knn_nhs1_cholhx_med","knn_nhs1_dbhx_med",
             "knn_nhs1_smk_status","knn_nhs1_smk_per_day","knn_nhs1_alc","knn_nhs1_ahei2010_noALC","knn_nhs1_act",
             "knn_nhs1_wtchg",
             "knn_nhs1_income_cens","knn_nhs1_percent_college_cens","knn_nhs1_percent_white_cens",
             "nhs1_region") 
cat.vars = c("nhs1_nonhisp","knn_nhs1_fast","knn_nhs1_menopmh_88_der",
             "knn_nhs1_hysthx","knn_nhs1_oophhx",
             "knn_nhs1_hypthx_med","knn_nhs1_cholhx_med","knn_nhs1_dbhx_med",
             "knn_nhs1_smk_status","knn_nhs1_smk_per_day",
             "nhs1_region")
cont.vars = setdiff(all.vars, cat.vars)

# Check distributions of categorical variables
apply(pdata[,cat.vars], 2, table, useNA="always")

# Check distributions of continuous variables
pdata[,cont.vars] %>%
  gather() %>% 
  ggplot(aes(value)) +
  facet_wrap(~ key, scales = "free") +
  geom_histogram(bins=20, na.rm=TRUE) 
# Use median [IQR] in table 1 given some skewed distributions


#### 2 - CREATE TABLE ####

# Create table 1
t1 <- CreateTableOne(vars = all.vars, strata = "Race", data = pdata, 
                    factorVars = cat.vars, test=FALSE)
t1.mat <- print(t1, nonnormal = cont.vars, format='fp', noSpaces=TRUE, contDigits=1, catDigits=1)
t1.mat

# Export table 1
write.csv(as.data.frame(t1.mat), file = file.path(main.dir,"final_results/tables/table_1.csv"), row.names=TRUE)

# Export table 1 footnotes - calculate n (%) missing
# Create list of variables for table 1
all.vars.not.imputed = c("nhs1_nonhisp","nhs1_agebld","nhs1_fast","nhs1_menopmh_88_der",
                                    "nhs1_hysthx","nhs1_oophhx",
                                    "nhs1_hypthx_med","nhs1_cholhx_med","nhs1_dbhx_med",
                                    "nhs1_smk_status","nhs1_smk_per_day","nhs1_alc","nhs1_ahei2010_noALC","nhs1_act",
                                    "nhs1_income_cens","nhs1_percent_college_cens","nhs1_percent_white_cens",
                                    "nhs1_wtchg") 
miss.percent <- apply(pdata[,all.vars.not.imputed], 2, function(x) round(length(which(is.na(x)))/length(x)*100,1))
miss.n <- apply(pdata[,all.vars.not.imputed], 2, function(x) round(length(which(is.na(x))),1))
miss.dat <- cbind(miss.n, miss.percent)
write.csv(as.data.frame(miss.dat), file = file.path(main.dir,"final_results/tables/table_1_footnote.csv"), row.names=TRUE)


# % missing data among Black vs. white women
pdata.black <- pdata[which(pdata$Race=="Black"),]
miss.percent <- apply(pdata.black[,all.vars.not.imputed], 2, function(x) round(length(which(is.na(x)))/length(x)*100,1))
miss.n <- apply(pdata.black[,all.vars.not.imputed], 2, function(x) round(length(which(is.na(x))),1))
miss.dat <- cbind(miss.n, miss.percent)
miss.dat

pdata.white <- pdata[which(pdata$Race=="White"),]
miss.percent <- apply(pdata.white[,all.vars.not.imputed], 2, function(x) round(length(which(is.na(x)))/length(x)*100,1))
miss.n <- apply(pdata.white[,all.vars.not.imputed], 2, function(x) round(length(which(is.na(x))),1))
miss.dat <- cbind(miss.n, miss.percent)
miss.dat


