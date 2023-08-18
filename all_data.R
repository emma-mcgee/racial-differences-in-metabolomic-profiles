##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

#########################################################################

# Program: all_data.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Study Aims:
  # Aim 1: quantify observed differences in individual metabolites and groups of molecularly or biologically similar metabolites between Black and White women
  # Aim 2: estimate residual differences in metabolomic profiles in a counterfactual population in which the distributions of established disease 
           # risk factors are equalized across racial groups

# Purpose of Program: Create datasets of 1) known metabolites (primary analyses) and 2) unknown metabolites/peaks (sensitivity analyses) for NHS

# Study Design: Longitudinal

# Study Population: NHS

# Population Selection:
# Eligibility Criteria
  # 1) blood sample provided at first blood collection (1989-1990)
  # 2) cancer-free at blood collection
  # 3) underwent metabolomic profiling (selected into a breast cancer metabolomics case-control study 
     # OR Black woman potentially eligible to be selected into such a study)
  # 4) self-identified as Black or White race 

# Descriptive variable for differences analyses:
  # RACE_BW: self-reported race, Black or White
    # Timing of assessment: derived variable based on questionnaire responses
    # Missing data: none
    # Categorical: 1:Black, 0:White

# Risk factors for counterfactual analysis:
  # KNN_NHS1_SMK_PER_DAY: average number of cigarettes smoked per day
    # Timing of assessment: 1986
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: 0, 1-14, 15+ cigs/day
  # KNN_NHS1_ALC: average alcohol consumption
    # Timing of assessment: 1986
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: grams/day, linear term for primary analyses, quadratic evaluated in sensitivity analysis
  # KNN_NHS1_AHEI2010_NOALC: AHEI diet score excluding alcohol (range: 0-100)
    # Timing of assessment: 1986
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: score, linear term for primary analyses
  # KNN_NHS1_ALC: average total physical activity per week
    # Timing of assessment: 1986
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: met-hours/week, linear term for primary analyses, quadratic evaluated in sensitivity analysis
  # KNN_NHS1_HYSTHX: history of hysterectomy 
    # Timing of assessment: prior to blood collection
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: yes, no
  # KNN_NHS1_OOPHHX: history of oophorectomy
    # Timing of assessment: prior to blood collection
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: yes, no
  # KNN_NHS1_MENOPMH_88_DER: menopausal status and PMH use prior to blood collection
    # Timing of assessment: 1988
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: premenopausal, postmenopausal no HRT, postmenopausal HRT
  # KNN_NHS1_LN_INCOME_CENS: census-tract median family income
    # Timing of assessment: residential address in 1986 + data from nearest census
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: ln(dollars), linear term for primary analyses, quadratic evaluated in sensitivity analysis
  # KNN_NHS1_PERCENT_COLLEGE_CENS: census-tract % of population over 25 with a college-level education
    # Timing of assessment: residential address in 1986 + data from nearest census
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: %, linear term for primary analyses, quadratic evaluated in sensitivity analysis
  # KNN_NHS1_PERCENT_WHITE_CENS: census-tract % White
    # Timing of assessment: residential address in 1986 + data from nearest census
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: %, linear term for primary analyses, quadratic evaluated in sensitivity analysis
  # KNN_NHS1_WTCHG: adult weight change (difference in weight between age 18 and blood collection)
    # Timing of assessment: age 18 to prior to blood collection
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: kg, linear term for primary analyses, quadratic evaluated in sensitivity analysis
  # KNN_NHS1_HYPTHX_MED: hypertension treatment and diagnosis history
    # Timing of assessment: 1988
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: diagnosed + treated, diagnosed + not treated, not diagnosed
  # KNN_NHS1_DBHX_MED: diabetes treatment and diagnosis history
    # Timing of assessment: 1988
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: diagnosed + treated, diagnosed + not treated, not diagnosed
  # KNN_NHS1_CHOLHX_MED: high cholesterol treatment and diagnosis history
    # Timing of assessment: 1988
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: diagnosed + treated, diagnosed + not treated, not diagnosed

# Outcome variable: 
  # Plasma metabolites measured on the C8-positive and HILIC-positive platforms
    # Exclusion criteria for main analyses:
       # 1) known metabolites that fail processing delay reproducibility pilot (rho <0.75 AND ICC <0.75)
       # 2) known and unknown metabolites that have >=10% missing values
    # Timing of assessment: first blood collection (1989-1990)
    # Missing data: imputed using k-Nearest Neighbors for metabolites with <10% missing
                    # Metabolites with >=10% missing data are analyzed in a sensitivity analysis
    # Continuous: natural log and Z-score transformed

# Stratifiers:
  # NHS1_REGION: region of residence in the US
    # Timing of assessment: 1988
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: Northeast, Midwest, South, West
  # NHS1_BIRTH_COHORT: birth cohort in decades
    # Timing of assessment: birth year
    # Missing data: N/A
    # Categorical: 20s, 30s, 40s

# Matching factors from original case-control study (accounted for via inverse probability weighting):
  # KNN_NHS1_AGEBLD: age at blood collection
    # Timing of assessment: at first blood collection (1989-1990)
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: years, linear + quadratic term
  # KNN_NHS1_MENOPMH_BLD: menopausal status and HRT use at blood collection
    # Timing of assessment: at first blood collection (1989-1990)
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: premenopausal, postmenopausal no HRT, postmenopausal HRT, unknown
  # KNN_NHS1_BLDMONTH: month of blood collection
    # Timing of assessment: at first blood collection (1989-1990)
    # Missing data: imputed using k-Nearest Neighbors
    # Continuous: months since 1900, linear term
  # KNN_NHS1_BTIMELR: time of day of blood collection
    # Timing of assessment: at first blood collection (1989-1990)
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: midnight-7am, 8am-noon, 1pm-11pm
  # KNN_NHS1_FAST: fasting status at blood collection
    # Timing of assessment: at first blood collection (1989-1990)
    # Missing data: imputed using k-Nearest Neighbors
    # Categorical: >=8 hours, <8 hours


# Statistical Analyses:
  # Descriptive statistics used to assess data distributions, missing values, etc.

##########################################################################


###### 1 - CREATE MERGED METABOLOMICS DATASET OF KNOWN AND UNKNOWN METABOLITES ######

rm(list = ls())
main.dir = ""
setwd(main.dir)

# Source merge_metab_data function and other packages
#devtools::install("/proj/nhmngs/nhmng00/R/chanmetab", build_vignettes = TRUE, upgrade = FALSE)
library(chanmetab)
library(Biobase)
library(tidyverse)

# Package documentation
#vignette("chanmetab")

# Create merged metabolomics dataset: nhs1 breast endpoint (labcodes 10114 and 10116) and nhs1 racial differences endpoint (labcode 10122)
# Include both known and unknown metabolites, which will later be subset for primary analyses (known metabolites) and secondary analyses (unknown metabolites)
# Using the row-name merge function because all samples were analyzed together in a single batch
merged.data <- row_name_merge(endpoints = c("breast", "racial.diff"),
                              labcodes = c("10114", "10116", "10122"),
                              cohorts = "nhs1",
                              methods = c("C8-pos", "HILIC-pos"),
                              collection_to_use = "first",
                              combine_cohorts = FALSE,
                              transformation = "transform_none",
                              impute_cutoff = 0.0,
                              controls_only = FALSE,
                              merge_type = "union",
                              extended_annotation = FALSE,
                              keep_failed_pm_metabolites = TRUE,
                              debug = TRUE,
                              keep_unknown_metabolites = TRUE)
# NOTE: Imputations, transformations, and dropping of PM failed metabolites will be done in later steps

# Check merged dataset
merged.data$summary
merged.eset <- merged.data$expr_set$nhs1

# Remove 6 metabolites that are duplicated across methods + 1 that should not be included (metformin) due to data quality issues
# + 2 metabolites which are compounds used by the lab for quality control purposes, not compounds we are measuring (valine-d8 and phenylalinine-d8)
# Duplicated on metabolite name
# Keeping duplicates that were measured on the preferred method (C8-positive platform, since all are lipid metabolites)
# Also removing metformin (HMDB001921, not reliably measured based on high #s of non-0 values)
source(file.path(main.dir,"function_removeMetabolitesFromESet.R"))
dup.names <- c("X701_5586mz_7_32","X764_6749mz_1_49","X790_6906mz_1_48","X792_7061mz_1_48","HMDB01921","HMDB07973_1","HMDB10382","X126_1367mz_7_55","X174_1365mz_6_86") 
test <-fData(merged.eset)
test[dup.names,c("hmdb_id","method","metabolite_name","mean_cv","mean_icc")] 
merged.eset = removeMetabolitesFromESet(eset = merged.eset, metabolites.to.remove = dup.names)
test <-fData(merged.eset) 
test[c("HMDB10419","X701_5602mz_7_93","X769_6328mz_10_85","X795_6481mz_10_96","HMDB07973","X518_3227mz_4_90")
     ,c("hmdb_id","method","metabolite_name","mean_cv","mean_icc")] # confirmed that correct metabolites were kept

# Remove 3 metabolites that are duplicated within method and were not removed from merged eset
# Duplicated on both HMDB ID and metabolite name
# Keeping duplicate with better laboratory CV
test[c("HMDB05385","HMDB05385_1","HMDB05391","HMDB05391_1","HMDB08039","HMDB08039_1"),c("hmdb_id","method","metabolite_name","mean_cv","mean_icc")] 
dup.names.2 <- c("HMDB05385","HMDB05391","HMDB08039_1")
merged.eset = removeMetabolitesFromESet(eset = merged.eset, metabolites.to.remove = dup.names.2)
test <-fData(merged.eset) 
test[c("HMDB05385_1","HMDB05391_1","HMDB08039"),c("hmdb_id","method","metabolite_name","mean_cv","mean_icc")] # confirmed that metabolites with lower CV were kept


# Save
save(merged.eset, file = file.path(main.dir, "data/merged.eset.RData"))



###### 2 - ADD COVARIATE DATA ######

# Source addCovariateData function
source(file.path(main.dir,"function_addCovariateData.R"))

# Import covariates
nhs1.cov <- read.csv(file = file.path(main.dir, "data/NHS1_covariate_data_blood_cohort.csv"), sep = ",", header = T, na.strings = "")

# Create new id variable with 6 digits - merged metabolomics dataset uses 7 digit ids but 6 digits are used elsewhere
merged.eset$old_id <- merged.eset$id
merged.eset$id <- substr(merged.eset$id, 1, 6)

# Check for non-overlapping ids between merged metabolomics dataset and covariate data
pdata <- pData(merged.eset)
non.overlapping.samples <- setdiff(pdata$id, as.character(nhs1.cov$id))
length(non.overlapping.samples)

# Add covariate data (NHS1 covariates from 1st blood collection)
merged.eset.raw <- addCovariateData(covariate.data = nhs1.cov, eset.data = merged.eset, prefix = "nhs1_",verbose = T)
save(merged.eset.raw, file = file.path(main.dir, "data/merged.eset.raw.RData"))


###### 3 - SUBSET PARTICIPANTS ######

# Check for participants with a history of cancer (excluding non-melanoma skin cancer) prior to blood draw
table(merged.eset.raw$nhs1_canhx, merged.eset.raw$nhs1_caco, useNA="always")
# 3 cases who were originally included were later recorded to have a prior cancer at blood
# Because this discrepancy is due to the dynamic nature of the NHS dataset, these three cases will be retained to retain the integrity of the *original* case-control selection criteria

# Subset to Black (nhs1_race=2) and White (nhs1_race=1) women, using derived race variable from CACO matrix
merged.eset.raw <- merged.eset.raw[ , (is.na(merged.eset.raw$nhs1_race)==FALSE) & (merged.eset.raw$nhs1_race == "1" | merged.eset.raw$nhs1_race == "2")]
dim(merged.eset.raw)

# Check final data
table(merged.eset.raw$nhs1_race, useNA="always")
table(merged.eset.raw$nhs1_caco, merged.eset.raw$race, useNA="always")

# Check for participants who are ineligible based on yobf (missing or <20)
summary(merged.eset.raw$nhs1_yobf, useNA="always") # no ineligible participants


###### 4 - LN Z-SCORE TRANSFORMATION OF METABOLITES & MISSING DATA IMPUTATION ######

# Impute missing data using K-nearest neighbors with Gower distances for mixed type datasets 
# Both missing metabolite and missing covariate data is imputed using K-nearest neighbors
# All covariates and metabolites are used in the calculation of nearest neighbors
# All metabolite values are ln z-score transformed to ensure a common scale and decrease the influence of extreme values


#### a) - PREPARE DATA FOR IMPUTATION ####

# Load packages and set seed
library(VIM)
library(laeken)
set.seed(8499023)

# Create outcome variable (black=1 vs. white=0) - using nhs1 race data from CACO matrix file
merged.eset.raw$race_bw <- recode(merged.eset.raw$nhs1_race, '1'="0", '2'="1")
table(merged.eset.raw$race_bw, useNA="always")

# Do natural log + z-score transformation of metabolites
# Need to scale all metabolites to same scale before imputation, since KNN works better when data is normalized/scaled
source(file.path(main.dir, "function_transformMetabolitesToLnZScores.R"))
merged.eset.raw2 <- transformMetabolitesToLnZScores(merged.eset.raw)

# Restrict to known metabolites for KNN imputation (unknowns will be imputed separately)
merged.eset.fdata <- fData(merged.eset.raw2)
table(merged.eset.fdata$known)
unknown.names <- rownames(merged.eset.fdata)[which(merged.eset.fdata$known==0)]
dim(merged.eset.raw2)
merged.eset.raw3 <- removeMetabolitesFromESet(eset = merged.eset.raw2, metabolites.to.remove = unknown.names)

# Extract metabolites for imputation
met.vals <- exprs(merged.eset.raw3)

# Restrict to metabolites with <10% missing
missing <- apply(met.vals, 1, FUN = function(x) return(length(which(is.na(x)))/length(x)*100))
table(missing)
missing.names <- rownames(met.vals)[which(missing < 10)]
length(missing.names)
met.vals.lt10miss <- met.vals[row.names(met.vals) %in% missing.names,]
dim(met.vals.lt10miss)

# Transpose so that metabolites are in columns and individuals are in rows
met.vals.t <- t(met.vals.lt10miss)
rownames(met.vals.t) <- colnames(met.vals)

# Add epidemiologic/demographic covariates to be used during imputation
pdata <- pData(merged.eset.raw3)

# Subset to variables to be used in imputation algorithm and covariates to be imputed
imput.var = c("id","race_bw",
              "nhs1_agebld","nhs1_menopmh_bld","nhs1_fast","nhs1_bldmonth","nhs1_btimelr",
              "nhs1_smk_per_day","nhs1_alc","nhs1_ahei2010_noALC","nhs1_act",
              "nhs1_oophhx","nhs1_hysthx",
              "nhs1_ln_income_cens","nhs1_percent_college_cens","nhs1_percent_white_cens",
              "nhs1_wtchg",
              "nhs1_menopmh_88_der",
              "nhs1_hypthx_med","nhs1_dbhx_med","nhs1_cholhx_med",
              "nhs1_smk_status","nhs1_smk_status_84","nhs1_smk_per_day_84","nhs1_alc_84","nhs1_ahei2010_noALC_84","nhs1_acth82",
              "nhs1_percent_pov_cens_76","nhs1_percent_unem_cens_76","nhs1_percent_college_cens_76","nhs1_percent_white_cens_76","nhs1_percent_int_cens_76",
              "nhs1_bmi18","nhs1_bmibld",
              "nhs1_ahei2010_frtI86","nhs1_ahei2010_naI86","nhs1_ahei2010_nutI86","nhs1_ahei2010_omegaI86",
              "nhs1_ahei2010_polyI86","nhs1_ahei2010_ptranI86","nhs1_ahei2010_rmtI86","nhs1_ahei2010_ssbI86",
              "nhs1_ahei2010_vegI86","nhs1_ahei2010_whgrnI86",
              "nhs1_menarc","nhs1_afbpar","nhs1_brfdhx",
              "nhs1_ses_score_cens",
              "nhs1_pred_fat_mass",
              "nhs1_bcfamhx","nhs1_mifamhx",
              "nhs1_state_88")
pdata <- pdata[,imput.var]

# Check missing values of covariates
map(pdata, ~mean(is.na(.))) 

# Combine covariate and metabolite data
dim(met.vals.t)
met.vals.t.rf <- cbind(pdata,met.vals.t)

# Convert all categorical variables to factors
cat.vars = c("race_bw",
             "nhs1_menopmh_bld","nhs1_fast","nhs1_btimelr",
             "nhs1_smk_per_day",
             "nhs1_oophhx","nhs1_hysthx",
             "nhs1_menopmh_88_der",
             "nhs1_hypthx_med","nhs1_dbhx_med","nhs1_cholhx_med",
             "nhs1_smk_status","nhs1_smk_status_84","nhs1_smk_per_day_84","nhs1_acth82",
             "nhs1_ahei2010_naI86", # treated as factor for imputation only, since this variable is discrete ordinal
             "nhs1_menarc", # treated as factor for imputation only, since this variable is discrete ordinal
             "nhs1_afbpar","nhs1_brfdhx",
             "nhs1_bcfamhx","nhs1_mifamhx",
             "nhs1_state_88")
met.vals.t.rf[,cat.vars] <- lapply(met.vals.t.rf[,cat.vars], factor) 
sapply(met.vals.t.rf[,cat.vars], class)


#### b) - IMPUTE MISSING DATA USING K-NEAREST NEIGHBORS ####

# Set number of nearest neighbors (k) to be the square root of the total number of observations per previous simulation study findings
k_val <- floor(sqrt(dim(met.vals.t.rf)[1])) # round down to ensure an odd number, which helps with breaking ties
k_val

# Do KNN imputation - variation of Gower distances used for mixed continuous and categorical data
met.vals.t.imputed <-kNN(met.vals.t.rf,
                         dist_var = colnames(met.vals.t.rf[,-1]), # Use all metabolites and risk factors for imputation algorithm, excluding ID
                         numFun = weightedMean, # Impute weighted mean of neighbor values for continuous data
                         catFun = maxCat, # Impute mode of neighbor values for categorical data
                         useImputedDist = TRUE, # Allow an imputed value to be used for distance calculation for imputing another variable (when =TRUE, this creates a dependency on ordering of variables)
                         weightDist=TRUE, # Use weighted distances in the aggregation step
                         k=k_val, # Set number of neighbors to square root of number of observations
                         imp_var=FALSE, # Do not create new variables with _imp suffixes
                         trace=FALSE)

# Check data
anyNA(met.vals.t.imputed)

#### c) - MERGE IMPUTED DATA WITH ORIGINAL DATA ####

# Add imputed covariates to the original dataset - excluding outcome variable (not imputed)
risk.factor.imputed <-met.vals.t.imputed[,c(1, 3:length(imput.var))]

# Check for non-overlapping ids between merged metabolomics dataset and covariate data
pdata <- pData(merged.eset.raw)
non.overlapping.samples <- setdiff(pdata$id, as.character(risk.factor.imputed$id))
length(non.overlapping.samples)

# Add imputed covariate data
dim(merged.eset.raw)
merged.eset.final <- addCovariateData(covariate.data = risk.factor.imputed, eset.data = merged.eset.raw, prefix = "knn_",verbose = T)
dim(merged.eset.final)
pdata <- pData(merged.eset.final)

# Add imputed metabolite data to the original eset

# Transpose metabolites and drop covariate data to prepare for merging
met.vals.imputed <- t(met.vals.t.imputed)
dim(met.vals.imputed)
met.vals.imputed <-met.vals.imputed[-c(1:length(imput.var)),]
dim(met.vals.imputed)
colnames(met.vals.imputed) <- rownames(met.vals.t)

# Subset known metabolites that were missing >=10% or were unknown (not imputed) to prepare for merging
# Knowns with >=10% missing data
met.vals.orig <- exprs(merged.eset.raw3)
missing <- apply(met.vals.orig, 1, FUN = function(x) return(length(which(is.na(x)))/length(x)*100))
missing.names <- rownames(met.vals)[which(missing >= 10)]
length(missing.names)
met.vals.gt10miss <- met.vals[row.names(met.vals) %in% missing.names,]
dim(met.vals.gt10miss)
# Unknowns
known.names <- rownames(merged.eset.fdata)[which(merged.eset.fdata$known==1)]
dim(merged.eset.final)
unknown.eset <- removeMetabolitesFromESet(eset = merged.eset.raw2, metabolites.to.remove = known.names)
met.vals.unknown <- exprs(unknown.eset)

# Combine metabolites with <10% missing and metabolites with >=10% missing
met.vals.imputed.final <- rbind(met.vals.imputed, met.vals.gt10miss, met.vals.unknown)
dim(met.vals.imputed.final)
class(met.vals.imputed.final) <- "numeric"

# Assign met.vals to dataset
exprs(merged.eset.final) <- met.vals.imputed.final
dim(merged.eset.final)
pdata <- pData(merged.eset.final)



#### d) - Z-SCORE TRANSFORM METABOLITES ####

# Final Z-score transformation after all data pre-processing, subsetting, and imputation - ensure mean=0 and SD = 1 after imputation in final sample
source(file.path(main.dir, "function_transformMetabolitesToZScores.R"))
merged.eset.final <- transformMetabolitesToZScores(merged.eset.final)
dim(merged.eset.final)

# Save imputed data
save(merged.eset.final, file = file.path(main.dir, "data/merged.eset.final.RData"))



###### 5 - CREATE DERIVED COVARIATES ######

# Create derived covariates for multivariate models, including quadratic terms for continuous variables and derived scores/prediction variables
# Relevel categorical variables so that category with largest N is the reference category for multivariate models, where appropriate


#### a) - VARIABLES USED AS MATCHING FACTORS IN ORIGINAL CASE-CONTROL STUDY ####
# Age at blood draw - create quadratic term
merged.eset.final$knn_nhs1_agebld_sqr <- merged.eset.final$knn_nhs1_agebld^2
summary(merged.eset.final$knn_nhs1_agebld, useNA="always")
summary(merged.eset.final$knn_nhs1_agebld_sqr, useNA="always")

# Fasting status (1=fasting, 0=non-fasting)
table(merged.eset.final$knn_nhs1_fast, useNA="always")
merged.eset.final$knn_nhs1_fast <- factor(merged.eset.final$knn_nhs1_fast, levels = c("1", "0"))
table(merged.eset.final$knn_nhs1_fast, useNA="always")

# Menopausal status at blood - relevel
table(merged.eset.final$knn_nhs1_menopmh_bld, useNA="always")
merged.eset.final$knn_nhs1_menopmh_bld <- factor(merged.eset.final$knn_nhs1_menopmh_bld, levels = c("2", "1", "3","4"))
summary(merged.eset.final$knn_nhs1_menopmh_bld)

# Date of blood draw (continuous, months since 1900) - create quadratic term
table(merged.eset.final$knn_nhs1_bldmonth, useNA="always")
merged.eset.final$knn_nhs1_bldmonth_sqr <- merged.eset.final$knn_nhs1_bldmonth^2
summary(merged.eset.final$knn_nhs1_bldmonth)
summary(merged.eset.final$knn_nhs1_bldmonth_sqr)

# Time of day of blood draw (1=midnight-7am, 2=8am-noon,3=1pm-11pm) - relevel
table(merged.eset.final$knn_nhs1_btimelr, useNA="always")
merged.eset.final$knn_nhs1_btimelr <- factor(merged.eset.final$knn_nhs1_btimelr, levels = c("2", "1", "3")) # change reference to category with largest N
table(merged.eset.final$knn_nhs1_btimelr)



#### b) - LIFESTYLE FACTORS ####

# Smoke status
table(merged.eset.final$knn_nhs1_smk_status, useNA="always")
table(merged.eset.final$knn_nhs1_smk_status_84, useNA="always")

# Smoke cigarettes per day - collapse categories due to small numbers for counterfactual analysis
table(merged.eset.final$knn_nhs1_smk_per_day, useNA="always")
pData(merged.eset.final)  <- mutate(pData(merged.eset.final), knn_nhs1_smk_per_day = factor(knn_nhs1_smk_per_day),
                knn_nhs1_smk_per_day = fct_recode(knn_nhs1_smk_per_day,
                                              "0" = "0",
                                              "1" = "1",
                                              "1" = "2",
                                              "2" = "3",
                                              "2" = "4",
                                              "2" = "5",
                                              "2" = "6"))
table(merged.eset.final$knn_nhs1_smk_per_day, useNA="always")
# Smoke cigaretters per day in 1984 does not need to be collapsed
table(merged.eset.final$knn_nhs1_smk_per_day_84, useNA="always")

# Alcohol consumption - create quadratic term
summary(merged.eset.final$knn_nhs1_alc, useNA="always")
merged.eset.final$knn_nhs1_alc_sqr <- merged.eset.final$knn_nhs1_alc^2
summary(merged.eset.final$knn_nhs1_alc_sqr)

# AHEI 2010 (excluding alcohol)
# AHEI total score (excluding alcohol) - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_noALC, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_noALC_sqr <- merged.eset.final$knn_nhs1_ahei2010_noALC^2
summary(merged.eset.final$knn_nhs1_ahei2010_noALC_sqr)

summary(merged.eset.final$knn_nhs1_ahei2010_noALC_84, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_noALC_84_sqr <- merged.eset.final$knn_nhs1_ahei2010_noALC_84^2
summary(merged.eset.final$knn_nhs1_ahei2010_noALC_84_sqr)

# AHEI-2010 (excluding alcohol) COMPONENTS
# Fruit - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_frtI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_frtI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_frtI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_frtI86_sqr)
# Sodium  - convert to numeric, create quadratic term
merged.eset.final$knn_nhs1_ahei2010_naI86 = as.numeric(as.character(merged.eset.final$knn_nhs1_ahei2010_naI86))
summary(merged.eset.final$knn_nhs1_ahei2010_naI86)
summary(merged.eset.final$knn_nhs1_ahei2010_naI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_naI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_naI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_naI86_sqr)
# Nuts - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_nutI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_nutI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_nutI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_nutI86_sqr)
# Omega - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_omegaI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_omegaI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_omegaI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_omegaI86_sqr)
# Poly - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_polyI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_polyI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_polyI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_polyI86_sqr)
# Ptran - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_ptranI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_ptranI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_ptranI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_ptranI86_sqr)
# RMT - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_rmtI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_rmtI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_rmtI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_rmtI86_sqr)
# SSB - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_ssbI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_ssbI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_ssbI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_ssbI86_sqr)
# Veg - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_vegI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_vegI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_vegI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_vegI86_sqr)
# Whole grain - create quadratic term
summary(merged.eset.final$knn_nhs1_ahei2010_whgrnI86, useNA="always")
merged.eset.final$knn_nhs1_ahei2010_whgrnI86_sqr <- merged.eset.final$knn_nhs1_ahei2010_whgrnI86^2
summary(merged.eset.final$knn_nhs1_ahei2010_whgrnI86_sqr)

# Physical activity - create quadratic term and relevel 1982 variable
summary(merged.eset.final$knn_nhs1_act, useNA="always")
merged.eset.final$knn_nhs1_act_sqr <- merged.eset.final$knn_nhs1_act^2
summary(merged.eset.final$knn_nhs1_act_sqr)

table(merged.eset.final$knn_nhs1_acth82, useNA="always")
merged.eset.final$knn_nhs1_acth82 <- factor(merged.eset.final$knn_nhs1_acth82, levels = c("less than 1 hour/week", "1 hour", "2-3 hours", "4-6 hours", "7 or more hours per week")) 
table(merged.eset.final$knn_nhs1_acth82)



#### c) - REPRODUCTIVE FACTORS ####

# History of oophorectomy
summary(merged.eset.final$knn_nhs1_oophhx)

# History of hysterectomy
summary(merged.eset.final$knn_nhs1_hysthx)

# Menopausal status and PMH
table(merged.eset.final$knn_nhs1_menopmh_88_der, useNA="always")
merged.eset.final$knn_nhs1_menopmh_88_der<- factor(merged.eset.final$knn_nhs1_menopmh_88_der, levels = c("2", "1", "3"))
summary(merged.eset.final$knn_nhs1_menopmh_88_der)

# Age at menarche
merged.eset.final$knn_nhs1_menarc = as.numeric(as.character(merged.eset.final$knn_nhs1_menarc))
summary(merged.eset.final$knn_nhs1_menarc)

# Parity/age at first birth - relevel
table(merged.eset.final$knn_nhs1_afbpar, useNA="always")
merged.eset.final$knn_nhs1_afbpar<- factor(merged.eset.final$knn_nhs1_afbpar, levels = c("4", "1", "2", "3", "5"))
summary(merged.eset.final$knn_nhs1_afbpar)

# History of breast feeding - relevel
table(merged.eset.final$knn_nhs1_brfdhx, useNA="always")
merged.eset.final$knn_nhs1_brfdhx <- factor(merged.eset.final$knn_nhs1_brfdhx, levels = c("1", "0"))
table(merged.eset.final$knn_nhs1_brfdhx, useNA="always")


#### d) - CONTEXTUAL FACTORS ####

# Median family income (natural log transformed) - create quadratic term
summary(merged.eset.final$knn_nhs1_ln_income_cens, useNA="always")
merged.eset.final$knn_nhs1_ln_income_cens_sqr <- merged.eset.final$knn_nhs1_ln_income_cens^2
summary(merged.eset.final$knn_nhs1_ln_income_cens_sqr)

# % over 25 w/ college education - create quadratic term
summary(merged.eset.final$knn_nhs1_percent_college_cens, useNA="always")
merged.eset.final$knn_nhs1_percent_college_cens_sqr <- merged.eset.final$knn_nhs1_percent_college_cens^2
summary(merged.eset.final$knn_nhs1_percent_college_cens_sqr)

# % White, census tract - create quadratic term
summary(merged.eset.final$knn_nhs1_percent_white_cens, useNA="always")
merged.eset.final$knn_nhs1_percent_white_cens_sqr <- merged.eset.final$knn_nhs1_percent_white_cens^2
summary(merged.eset.final$knn_nhs1_percent_white_cens_sqr)

# SES score, census tract - score representing census tract socioeconomic advantage and affluence, defined based on prior literature
summary(merged.eset.final$knn_nhs1_ses_score_cens, useNA="always")
merged.eset.final$knn_nhs1_ses_score_cens_sqr <- merged.eset.final$knn_nhs1_ses_score_cens^2
summary(merged.eset.final$knn_nhs1_ses_score_cens_sqr)


##### e) - ANTHROPOMETRIC FACTORS ####

# Weight change from age 18 to blood - create quadratic term
summary(merged.eset.final$knn_nhs1_wtchg, useNA="always")
merged.eset.final$knn_nhs1_wtchg_sqr <- merged.eset.final$knn_nhs1_wtchg^2
summary(merged.eset.final$knn_nhs1_wtchg_sqr)

# BMI at age 18 - create quadratic term
summary(merged.eset.final$knn_nhs1_bmi18, useNA="always")
merged.eset.final$knn_nhs1_bmi18_sqr <- merged.eset.final$knn_nhs1_bmi18^2
summary(merged.eset.final$knn_nhs1_bmi18_sqr)

# BMI at blood - create quadratic term
summary(merged.eset.final$knn_nhs1_bmibld, useNA="always")
merged.eset.final$knn_nhs1_bmibld_sqr <- merged.eset.final$knn_nhs1_bmibld^2
summary(merged.eset.final$knn_nhs1_bmibld_sqr)

# Predicted fat mass based on an anthropometric prediction equation developed and validated in NHANES: https://pubmed.ncbi.nlm.nih.gov/29110742/
summary(merged.eset.final$knn_nhs1_pred_fat_mass)
merged.eset.final$knn_nhs1_pred_fat_mass_sqr <- merged.eset.final$knn_nhs1_pred_fat_mass^2
summary(merged.eset.final$knn_nhs1_pred_fat_mass_sqr)

#### f) - HYPERTENSION HISTORY AND MEDICATIONS ####
summary(merged.eset.final$knn_nhs1_hypthx_med)

#### g) - DIABETES HISTORY AND MEDICATIONS ####
summary(merged.eset.final$knn_nhs1_dbhx_med)

#### h) - CHOLESTEROL HISTORY AND MEDICATIONS ####
summary(merged.eset.final$knn_nhs1_cholhx_med)

#### i) - FAMILY HEALTH HISTORY ####

# Family history of breast cancer
summary(merged.eset.final$knn_nhs1_bcfamhx)

# Family history of MI
summary(merged.eset.final$knn_nhs1_mifamhx)


#### j) - STRATIFIERS ####

# Create geographic region variable based on US census bureau designations and state of residence in 1986
table(merged.eset.final$knn_nhs1_state_88, useNA="always")
pData(merged.eset.final) <- mutate(pData(merged.eset.final), nhs1_region = factor(knn_nhs1_state_88),
                                   nhs1_region = fct_recode(nhs1_region,
                                                            "1=Northeast" = "CT","1=Northeast" = "ME","1=Northeast" = "MA",
                                                            "1=Northeast" = "NH","1=Northeast" = "RI","1=Northeast" = "VT",
                                                            "1=Northeast" = "NJ","1=Northeast" = "NY","1=Northeast" = "PA",
                                                            "2=Midwest" = "IN","2=Midwest" = "IL","2=Midwest" = "MI",
                                                            "2=Midwest" = "OH","2=Midwest" = "IA","2=Midwest" = "MN",
                                                            "2=Midwest" = "MO","2=Midwest" = "NE","2=Midwest" = "ND",
                                                            #"2=Midwest" = "WI","2=Midwest" = "KS","2=Midwest" = "SD",
                                                            "3=South" = "DE","3=South" = "DC","3=South" = "FL",
                                                            "3=South" = "GA","3=South" = "MD","3=South" = "NC",
                                                            "3=South" = "SC","3=South" = "VA","3=South" = "AL",
                                                            "3=South" = "KY","3=South" = "MS","3=South" = "TN",
                                                            "3=South" = "AR","3=South" = "OK","3=South" = "TX",
                                                            #"3=South" = "WV","3=South" = "LA",
                                                            "4=West" = "AZ","4=West" = "CO","4=West" = "ID",
                                                            "4=West" = "UT", "4=West" = "CA","4=West" = "OR",
                                                            #"4=West" = "NM", #"4=West" = "MT","4=West" = "NV",
                                                            #"4=West" = "WY","4=West" = "AK","4=West" = "HI",
                                                            "4=West" = "WA"),
                                   nhs1_region = fct_relevel(nhs1_region,"1=Northeast","2=Midwest","3=South","4=West"))
table(merged.eset.final$nhs1_region, useNA="always")


# Create variable for birth cohort in decades (based on year of birth, recorded in years since 1900)
table(merged.eset.final$nhs1_yobf, useNA="always")
pData(merged.eset.final) <- mutate(pData(merged.eset.final), nhs1_birth_cohort = factor(nhs1_yobf),
                                   nhs1_birth_cohort = fct_recode(nhs1_birth_cohort,
                                                                  "20s" = "21","20s" = "22","20s" = "23","20s" = "24","20s" = "25",
                                                                  "20s" = "26","20s" = "27","20s" = "28","20s" = "29",
                                                                  "30s" = "30","30s" = "31","30s" = "32","30s" = "33","30s" = "34",
                                                                  "30s" = "35","30s" = "36","30s" = "37","30s" = "38","30s" = "39",
                                                                  "40s" = "40","40s" = "41","40s" = "42","40s" = "43","40s" = "44",
                                                                  "40s" = "45","40s" = "46",),
                                   nhs1_birth_cohort = fct_relevel(nhs1_birth_cohort,"20s","30s","40s"))
table(merged.eset.final$nhs1_birth_cohort, useNA="always")


###### 6 - CREATE IP WEIGHTS TO ACCOUNT FOR CASE-CONTROL SELECTION ######

# Merge entire blood cohort with racial differences data
pdata <- pData(merged.eset.final)
not.included  <- setdiff(as.character(nhs1.cov$id),pdata$id)
length(not.included)
nhs1.cov.not <- nhs1.cov[which(nhs1.cov$id %in% not.included),]
nhs1.cov.include <- pdata

# Create inclusion indicator
nhs1.cov.not$Included <- 0
nhs1.cov.include$Included <- 1

# Subset women who are not included to those who are potentially eligible for this study (no prior history of cancer at blood collection, Black or White race)
nhs1.cov.not <- nhs1.cov.not[(nhs1.cov.not$canhx == "0"),]
nhs1.cov.not <- nhs1.cov.not[(is.na(nhs1.cov.not$race_bldq)==FALSE & (nhs1.cov.not$race_bldq == "1" | nhs1.cov.not$race_bldq == "2")),]

# Prepare data for merging
imput.var = c("id","race_bw",
              "knn_nhs1_agebld","knn_nhs1_agebld_sqr","knn_nhs1_menopmh_bld","knn_nhs1_fast","knn_nhs1_bldmonth","knn_nhs1_btimelr",
              "nhs1_lasteat",
              "Included","nhs1_caco","nhs1_canhx","nhs1_dxmonth","nhs1_bldmonth",
              "nhs1_can90","nhs1_can92","nhs1_can94","nhs1_can96","nhs1_can98",
              "nhs1_can00","nhs1_can02","nhs1_can04","nhs1_can06","nhs1_can08",
              "nhs1_can10","nhs1_can12","nhs1_can14")
nhs1.cov.include <- nhs1.cov.include[,imput.var]
nhs1.cov.not <- nhs1.cov.not %>% rename_with( ~ paste("nhs1", .x, sep = "_"))
nhs1.cov.not$knn_nhs1_agebld <- nhs1.cov.not$nhs1_agebld
nhs1.cov.not$knn_nhs1_agebld_sqr <- nhs1.cov.not$nhs1_agebld^2
nhs1.cov.not$knn_nhs1_bldmonth <- nhs1.cov.not$nhs1_bldmonth
table(nhs1.cov.not$nhs1_race_bldq, useNA="always")
nhs1.cov.not$race_bw <- nhs1.cov.not$nhs1_race_bldq - 1
table(nhs1.cov.not$race_bw, useNA="always")
nhs1.cov.not$id <- nhs1.cov.not$nhs1_id
nhs1.cov.not$Included <- nhs1.cov.not$nhs1_Included
nhs1.cov.not$knn_nhs1_menopmh_bld <- nhs1.cov.not$nhs1_menopmh_bld_der
nhs1.cov.not$knn_nhs1_btimelr <- nhs1.cov.not$nhs1_btimelr_bldq
nhs1.cov.not$knn_nhs1_fast <- ifelse(nhs1.cov.not$nhs1_lasteat>=8,1,
                                     ifelse(nhs1.cov.not$nhs1_lasteat<8,0,NA))
nhs1.cov.not <- nhs1.cov.not[,imput.var]


# Merge data
all_data <- rbind(nhs1.cov.include, nhs1.cov.not)


#### a) - IP weights for cases ####

# Subset to women who were eligible for selection as cases
# I.e., women who developed breast cancer after blood collection and by 2016 at the latest

# Create new dx year variable for entire blood cohort
all_data <- mutate(all_data, dxyear = nhs1_dxmonth/12)
summary(all_data$dxyear)

# Subset to eligible women with a date of diagnosis of breast cancer after blood collection but no later than 2016
potential_cases_final <- filter(all_data, nhs1_dxmonth < 1404 & nhs1_dxmonth > knn_nhs1_bldmonth)

# Fit a model for the log odds of being selected as a case (c_case==0)
potential_cases_final <- mutate(potential_cases_final, c_case = case_when(nhs1_caco== 1 ~ 0,
                                                                          TRUE ~ 1)) #all other cases
summary(as.factor(potential_cases_final$c_case), useNA="always")
case_denom <- glm(c_case==0~ 1, data=potential_cases_final, family=binomial(link='logit'))
summary(case_denom)

# Create predicted probabilities from log odds
potential_cases_final$ipw_denom <- predict(case_denom, potential_cases_final, type="response")

# Create weights
potential_cases_final$ipw_c_case <- (1/potential_cases_final$ipw_denom)
summary(potential_cases_final$ipw_c_case)



#### b) - IP weights for controls ####

# Parametrically estimate IP weights for controls to improve efficiency
potential_controls <- all_data

# Fit a model for the log odds of being selected as a control (c_control==0), conditional on the matching factors, which includes age and month of blood draw as functions of time (to account for risk-set sampling)
potential_controls <- mutate(potential_controls, c_control = case_when(nhs1_caco== 2 ~ 0,
                                                                       TRUE ~ 1))
table(potential_controls$c_control)
control_denom <- glm(c_control==0 ~ knn_nhs1_agebld + knn_nhs1_agebld_sqr + knn_nhs1_menopmh_bld + knn_nhs1_bldmonth + knn_nhs1_btimelr + knn_nhs1_fast,
                     data=potential_controls, family=binomial(link='logit'))
summary(control_denom)

# Create predicted probabilities from log odds
potential_controls$ipw_denom <- predict(control_denom, potential_controls, type="response")

# Create weights
potential_controls$ipw_c_control <- (1/potential_controls$ipw_denom)

# Check weights for controls included in this analysis
controls <- filter(potential_controls, nhs1_caco == 2)
summary(controls$ipw_c_control)



#### c) - IP weights for Black women not selected in the case-control study ####
potential_black <- filter(all_data, race_bw == 1 & is.na(nhs1_caco)==TRUE) # restrict to Black women not selected as a case or control

# Fit model for log odds of being selected as a Black woman (Included ==1)
black_denom <- glm(Included==1~ 1, data=potential_black, family=binomial(link='logit'))
summary(black_denom)

# Create predicted probabilities from log odds
potential_black$ipw_denom <- predict(black_denom, potential_black, type="response")

# Create weights
potential_black$ipw_c_black <- (1/potential_black$ipw_denom)
summary(potential_black$ipw_c_black)



#### d) -  Assign weights in final dataset ####
potential_cases_final <- potential_cases_final[,c("id","ipw_c_case")]
potential_controls <- potential_controls[,c("id","ipw_c_control")]
potential_black <- potential_black[,c("id","ipw_c_black")]
merged <- pdata %>% 
  left_join(potential_cases_final, by = "id") %>% 
  left_join(potential_controls, by="id") %>%
  left_join(potential_black, by = "id")
merged$ipw_c <- ifelse(merged$nhs1_caco==1, merged$ipw_c_case,
                       ifelse(merged$nhs1_caco==2, merged$ipw_c_control, NA))
merged$ipw_c <- ifelse(merged$race_bw=="1" & is.na(merged$ipw_c)==TRUE, merged$ipw_c_black, merged$ipw_c)
summary(merged$ipw_c)
keep.var = c("id","ipw_c")

# Add IPW to pdata
merged.eset.final <- addCovariateData(covariate.data = merged[,keep.var], eset.data = merged.eset.final, prefix = "",verbose = T)
summary(merged.eset.final$ipw_c) # Could consider truncation of extreme weights to improve efficiency


###### 7 - SUBSET TO KNOWN METABOLITES ######

# Subset to known metabolites (metabolites with names)
merged.eset.fdata <- fData(merged.eset.final)
unknown.names <- rownames(merged.eset.fdata)[which(merged.eset.fdata$known==0)]
dim(merged.eset.final)
known.eset.final <- removeMetabolitesFromESet(eset = merged.eset.final, metabolites.to.remove = unknown.names)

# Subset known metabolites to those that pass the processing delay pilot (heparin PM flag NOT = 1)
known.eset.fdata <- fData(known.eset.final)
table(known.eset.fdata$pm_pilot_fail_heparin, useNA="always")
# Ensure that PM flag agrees with manual assessment (rho <0.75 AND ICC <0.75)
pm.fail.check <- filter(known.eset.fdata, pm_pilot_rho_heparin <0.75, pm_pilot_icc_heparin <0.75)
length(pm.fail.check$metabolite_name)
# Check names of metabolites that don't pass processing delay
pm.fail <- known.eset.fdata[ which(known.eset.fdata$pm_pilot_fail_heparin=='1'), ]
table(pm.fail$metabolite_name, useNA="always")
# Exclude metabolites that don't pass processing delay pilot
fail.process.delay.names = rownames(known.eset.fdata)[which(known.eset.fdata$pm_pilot_fail_heparin=='1')]
dim(known.eset.final)
known.eset.final = removeMetabolitesFromESet(eset = known.eset.final, metabolites.to.remove = fail.process.delay.names)

# Check for any metabolites that were not measured in the protocol data
merged.data.prot = protocolData(known.eset.final)
measured <- t(pData(merged.data.prot))
not_measured <- measured[rowSums("not-measured" == measured) > 0,]
length(row.names(not_measured)) # none are not-measured

# SAVE FINAL DATA - KNOWNS
save(known.eset.final, file = file.path(main.dir, "data/known.eset.final.RData"))

###### 8 - SUBSET TO UNKNOWN METABOLITES ######

# Subset to unknown metabolites (metabolites without names)
known.names <- rownames(merged.eset.fdata)[which(merged.eset.fdata$known==1)]
dim(merged.eset.final)
unknown.eset.final <- removeMetabolitesFromESet(eset = merged.eset.final, metabolites.to.remove = known.names)

# SAVE FINAL DATA - UNKNOWNS
save(unknown.eset.final, file = file.path(main.dir, "data/unknown.eset.final.RData"))
