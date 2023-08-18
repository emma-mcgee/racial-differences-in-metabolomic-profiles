##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_5.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Descriptive characteristics of study population, stratified by race, WHI replication dataset (Supplemental Table 5)

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

# Load WHI data
source(file = file.path(main.dir,"whi_data.R"))
pdata <- whi.pdata

# Numbers in the text
# % non-Hispanic

# Age - median [IQR]
round(summary(pdata$age),1)

# Create data labels for race and for categorical variables with >2 categories

# Race
table(pdata$race_bw, useNA="always")
pdata <- mutate(pdata, Race = factor(race_bw),
                Race = fct_recode(Race,"Black" = "1","White" = "0"),
                Race = fct_relevel(Race,"Black","White"))
summary(pdata$Race)

# Age
summary(pdata$age, useNA="always")

# Fasting
table(pdata$fast, useNA="always")
pdata <- mutate(pdata, fast = factor(fast),
                fast = fct_relevel(fast,
                                "0","1"))
summary(pdata$fast)

# All women are menopausal
pdata$meno <- "1"

# Menopausal status in questionnaire cycle prior to blood collection
table(pdata$ht, useNA="always")

# Hysterectomy
table(pdata$hyst, useNA="always")

# Current smoker (binary for the table)
table(pdata$smoking, useNA="always")
pdata <- mutate(pdata, smoking_ever = factor(smoking),
                smoking_ever = fct_recode(smoking_ever,
                                              "Never" = "Never Smoked",
                                              "Former" = "Past Smoker",
                                              "Current" = "Current Smoker"),
                smoking_ever = fct_relevel(smoking_ever,
                                  "Never",
                                  "Former",
                                  "Current"))
table(pdata$smoking_ever)

# Cigarettes smoked per day
table(pdata$cigsday_new, useNA="always")
pdata <- mutate(pdata, cigsday = factor(cigsday_new),
                cigsday = fct_recode(cigsday,
                                             "0 cigarettes/day" = "0",
                                             "1-14 cigarettes/day" = "Less than 1",
                                             "1-14 cigarettes/day" = "1-4",
                                             "1-14 cigarettes/day" = "5-14",
                                             "15+ cigarettes/day" = "15-24",
                                             "15+ cigarettes/day" = "25-34",
                                             "15+ cigarettes/day" = "35-44",
                                             "15+ cigarettes/day" = "45 or more"))
table(pdata$cigsday, useNA="always")

# History of hypertension (excluding treatment)
table(pdata$hypt_trt, useNA="always")
pdata <- mutate(pdata, hypt = factor(hypt_trt),
                hypt = fct_recode(hypt,
                                             "No" = "0",
                                             "Yes" = "1",
                                             "Yes" = "2"),
                hypt = fct_relevel(hypt,
                                              "No",
                                              "Yes"))
table(pdata$hypt, useNA="always")

# Huperlipidemia treatment
table(pdata$hyperlipid, useNA="always")


# History of diabetes and treatment
table(pdata$diab_trt, useNA="always")
pdata <- mutate(pdata, diab = factor(diab_trt),
                diab = fct_recode(diab,
                                             "No" = "0",
                                             "Yes" = "1",
                                             "Yes" = "2"),
                diab = fct_relevel(diab,
                                              "No",
                                              "Yes"))
table(pdata$diab, useNA="always")

# Alcohol
summary(pdata$alcohol)

# Diet
summary(pdata$ahei)

# Physical activity
summary(pdata$phyact)

# BMI
summary(pdata$bmi)

# Income
table(pdata$income, useNA="always")
pdata <- mutate(pdata, income_cat = factor(income),
                income_cat = fct_recode(income_cat,
                                     "<$20,000" = "Less than $10,000",
                                     "<$20,000" = "$10,000 to $19,999",
                                     "$20,000 - <$35,000" = "$20,000 to $34,999",
                                     "$35,000 - <$50,000" = "$35,000 to $49,999",
                                     "$50,000 - <$100,000" = "$50,000 to $74,999",
                                     "$50,000 - <$100,000" = "$75,000 to $99,999",
                                     ">$100,000" = "$100,000 to $149,999",
                                     ">$100,000" = "$150,000 or more"),
                income_cat = fct_relevel(income_cat,
                                   "<$20,000",
                                   "$20,000 - <$35,000",
                                   "$35,000 - <$50,000",
                                   "$50,000 - <$100,000",
                                   ">$100,000"))
table(pdata$income_cat, useNA="always")

# Education
table(pdata$educ, useNA="always")
pdata <- mutate(pdata, educ_cat = factor(educ),
                educ_cat = fct_recode(educ_cat,
                                        "1" = "College graduate or Baccalaureate Degree",
                                        "1" = "Doctoral Degree (Ph.D,M.D.,J.D.,etc.)",
                                        "0" = "Grade school (1-4 years)",
                                        "0" = "Grade school (5-8 years)" ,
                                        "0" = "High school diploma or GED",
                                        "1" = "Master's Degree",
                                        "0" = "Some college or Associate Degree",
                                        "0" = "Some high school (9-11 years)",
                                        "1" = "Some post-graduate or professional",
                                        "0" = "Vocational or training school"))
table(pdata$educ_cat, useNA="always")

# Region
table(pdata$region, useNA="always")
pdata <- mutate(pdata, region = factor(region),
                region = fct_relevel(region,
                                   "Northeast",
                                   "Midwest",
                                   "South",
                                   "West"))
table(pdata$region, useNA="always")


# Create list of variables for table 1
all.vars = c("age","fast","meno","ht",
             "hyst",
             "hypt","diab","hyperlipid",
             "smoking_ever","cigsday","alcohol","ahei","phyact",
             "bmi",
             "income_cat",
             "educ_cat",
             "region") 
cat.vars = c("fast","meno","ht",
             "hyst",
             "hypt","diab","hyperlipid",
             "smoking_ever","cigsday",
             "income_cat",
             "educ_cat",
             "region")
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

# Create table
t1 <- CreateTableOne(vars = all.vars, strata = "Race", data = pdata, 
                    factorVars = cat.vars, test=FALSE)
t1.mat <- print(t1, nonnormal = cont.vars, format='fp', noSpaces=TRUE, contDigits=1, catDigits=1)

# Export table
write.csv(as.data.frame(t1.mat), file = file.path(main.dir,"final_results/supplemental_tables/supp_table_5.csv"), row.names=TRUE)
