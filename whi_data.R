##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: whi_data.R

# Programmer: Raji Balasubramanian & Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Create WHI dataset

# Statistical Analyses:
    # Descriptive statistics used to assess data distributions, missing values, etc.

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - LOAD AND FORMAT WHI DATA ####
rm(list = ls())
main.dir <- ""
results.dir <- file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(stats)
library(readxl)
library(Biobase)
library(tidyverse)
library(openxlsx)

# Create a function to load WHI data

##################################################################################################
# Function: selects all baseline data (visit = BASE) for 2306 participants in the WHI-OS and HT 
## VISIT = "Base" 
## ARM = > table(mx.dat$ARM)
#
#  E.CONTR   E.INTER E.P.CONTR E.P.INTER   OS 
#      369       363       278       352   944 

## OS=observational study
## E.CONTR = placebo arm in the estrogen only trial; E.INTER= active treatment (estrogen) in the estrogen only trial
## E.P.CONTR = placebo arm in the estrogen+progestin trial; E.P.INTER= active treatment (estrogen+progestin) in the estrogen+progestintrial

## CACO: 1= CHD case; 0 = control
##################################################################################################

select.base <- function(dset.dir, dt, by.x, by.y)
{
  select.base.partic <- function(dset.dir)
  {
    partic    <- read.csv(file=paste(dset.dir, "metab_particip_map.csv", sep=""), sep=",", as.is=T, header=T)
    partic.os.disc <- partic[partic$VISIT == "Base",]      
    return(subset(partic.os.disc, select=c("WHI.ID", "Complete.Draw.ID", "CACO", "VISIT", "ARM")))
  }
  
  os.disc.part <- select.base.partic(dset.dir)
  dt <- merge(os.disc.part[, c("Complete.Draw.ID", "WHI.ID", "CACO")], dt, by.x=by.x, by.y=by.y) 
  
  return(dt)
}

# Set dataset directories 
dset.dir <- ""
dset.dir2 <- ""

# Read metabolomics datasets 
newmx.dat <- read.csv(file = "MergedMxData-01252019.txt", sep="\t", header=T)

# Read dataset with metabolite annotations -- metabolite name, platform, hmdbid (if available)
mxnames3 <- read_excel("MergedMxInfo-10232018.xlsx")
mxnames4 <- as.data.frame(mxnames3)

# Read datasets with covariates:  in /proj/ncdats/ncdat02/WHI/derived
oth.dat     <- read.csv(file = paste(dset.dir,   "whibsln_all.csv", sep=""),   sep=",", header=T, as.is=T)
markers.dat <- read.csv(file = paste(dset.dir,   "biomarkers_bsln.csv", sep=""), sep=",", header=T, as.is=T)
new.dat     <- read.csv(file = paste(dset.dir,   "meds_all.csv", sep=""),      sep=",", header=T, as.is=T)
draw <- read.csv(file = paste(dset.dir, "draws_info.csv", sep=""), sep=",", header=T, as.is=T)

# Select 2306 participants with baseline metabolomics (visit==base)
# use these data for analyses (mx.dat = metabolite data; all other covariates are in the other datasets)
mx.dat        <- select.base(dset.dir2, newmx.dat,      "Complete.Draw.ID","Complete.Draw.ID")
oth.dat       <- select.base(dset.dir2, oth.dat,     "WHI.ID", "ID")
markers.dat   <- select.base(dset.dir2, markers.dat, "WHI.ID", "WHI.ID")
new.dat       <- select.base(dset.dir2, new.dat,     "WHI.ID", "ID")
draw <- select.base(dset.dir2, draw, "Complete.Draw.ID","Complete.Draw.ID")

# Put the datasets in the same order as the mx data 
whi.id   <- mx.dat[["WHI.ID.x"]]

oth.dat2 <- c()
markers.dat2 <- c()
new.dat2 <- c()
draw2 <- c()
for (i in 1:length(whi.id)){
  ind <- which(oth.dat[["WHI.ID"]] %in% whi.id[i])
  ind4 <- which(markers.dat[["WHI.ID"]] %in% whi.id[i])
  newind <- which(new.dat[["WHI.ID"]] %in% whi.id[i])
  ind6 <- which(draw[["WHI.ID.x"]] %in% whi.id[i])
  
  oth.dat2 <- rbind(oth.dat2, oth.dat[ind,]) 
  markers.dat2 <- rbind(markers.dat2, markers.dat[ind4,])
  new.dat2 <- rbind(new.dat2, new.dat[newind,])
  draw2 <- rbind(draw2, draw[ind6,])
}

# Check for any non-overlapping IDs
sum(oth.dat2[,1] != whi.id)
sum(markers.dat2[,1] != whi.id)
sum(new.dat2[,1] != whi.id) 
sum(draw2[,2] != whi.id)

# Drop participants who do not identify as White or Black

# Extract eligible IDs who identify as Black or White
oth.dat.wb <- filter(oth.dat2, ETHNIC == "Black or African-American" | ETHNIC == "White (not of Hispanic origin)")
race.id <- oth.dat.wb[["WHI.ID"]]
length(race.id)
# Total eligible participants = 2128

# Subset each dataset
oth.dat2 <- filter(oth.dat2, WHI.ID %in% race.id)
new.dat2 <- filter(new.dat2, WHI.ID %in% race.id)
draw2 <- filter(draw2, WHI.ID.x %in% race.id)
mx.dat$WHI.ID <- mx.dat$WHI.ID.x
mx.dat <- filter(mx.dat, WHI.ID %in% race.id)


#### 2 - CREATE WHI FDATA WITH METABOLITES THAT WERE ANALYZED IN NHS ####

# Load eset of known metabolites from NHS
load(file.path(main.dir, "data/known.eset.final.RData"))

# Remove metabolites that are missing >10% of all values
# Main analyses are restricted to metbaolites with <10% missing
# Only main analyses will be replicated in WHI
met.vals <- exprs(known.eset.final)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
miss.names <- row.names(met.vals[which(missing > 10),])
source(file.path(main.dir,"function_removeMetabolitesFromESet.R"))
known.eset.final <- removeMetabolitesFromESet(eset = known.eset.final, metabolites.to.remove = miss.names)

# Prepare data for merging

# NHS metabolite data
met.vals <- exprs(known.eset.final)
fdata <- fData(known.eset.final)
myvars <- c("class_broad","hmdb_id","metabolite_name","method")
metab.fdata <- fdata[myvars]
# Recode PC plasmalogen 
metab.fdata$class_broad <- ifelse(metab.fdata$class_broad=="PC plasmalogens", "Phosphatidylcholine plasmalogens", metab.fdata$class_broad)
metab.fdata$HMDB.ID <- row.names(metab.fdata)

# Check for duplicate metabolite names in NHS
metab.fdata[duplicated(metab.fdata$metabolite_name),]
metab.fdata[metab.fdata$metabolite_name %in% metab.fdata[duplicated(metab.fdata$metabolite_name),]$metabolite_name,] 
# No duplicates

# WHI metabolite data
mxnames4$metabolite_name <- mxnames4$Metabolite
table(mxnames4$Method)
# total number of metabolites from C8-pos and HILIC-pos platforms 164 + 244 = 408

# Subset to two platforms included in NHS (C8-pos and HILIC-pos)
mxnames.final <- filter(mxnames4, Method == "C8-pos" | Method == "HILIC-pos")
table(mxnames.final$Method)

# Check for duplicate metabolite names
mxnames.final[duplicated(mxnames.final$metabolite_name),] # none

# Check for duplicate HMDBIDs
# mxnames.final[duplicated(mxnames.final$HMDB.ID),] # some duplicated HMDBIDs, mostly NA

# Match NHS and WHI data by metabolite names
# Merging on names first since names were more uniformly documented in the early data than HMDB IDs were
merged.names <- merge(metab.fdata, mxnames.final, by="metabolite_name")
dim(merged.names)[1]
# 267 metabolites

# Subset to metabolites not matched on name
metab.fdata.id <- filter(metab.fdata, !metabolite_name %in% merged.names$metabolite_name)
mxnames.final.id <- filter(mxnames.final, !metabolite_name %in% merged.names$metabolite_name)

# Add metabolites that match on HMDB ID but not on name
mxnames.final.id$HMDB.ID <- substr(mxnames.final.id$HMDB.ID, 1, 9)
merged.names.id <- merge(metab.fdata.id, mxnames.final.id, by="HMDB.ID")

# Remove HMDB ID matches which are duplicated (n=3 after removing duplicates)
# There are A/B versions of these metabolites in WHI, but not in NHS - we do not know which is the appropriate match
merged.names.id.final <- filter(merged.names.id, !HMDB.ID %in% merged.names.id[duplicated(merged.names.id$HMDB.ID),]$HMDB.ID)

# Check metabolites not matched on name or HMDB ID - checking for synonyms that were not matched
mxnames.final.no.match <- filter(mxnames.final.id, !HMDB.ID %in% merged.names.id.final$HMDB.ID)
write.xlsx(mxnames.final.no.match, file = file.path(main.dir,"final_results/data_checks/unmatched_whi_metabolites.xlsx"),
           rowNames=FALSE, overwrite = TRUE)

# Add metabolites that match on synonyms (n=2)
mxnames.final.id["metabolite_name"][ mxnames.final.id["metabolite_name"] == "C16:0 LPC"] <- "Na_C16:0 LPC"
mxnames.final.id["metabolite_name"][ mxnames.final.id["metabolite_name"] == "guanidoacetic acid"] <- "guanidinoacetic acid"
# mxnames.final.id["metabolite_name"][ mxnames.final.id["metabolite_name"] == "C18:1 LPC plasmalogen"] <- "C18:1 LPC plasmalogen_minor"
# the above metabolite cannot be confirmed to be a match (also note this metab was measured on different platforms)

merged.names.syn.final <- merge(metab.fdata.id, mxnames.final.id, by="metabolite_name")

# Combine metabolites matched on name, ID, or synonym

# Format merged.syn.final
merged.names.syn.final$metabolite_name <- merged.names.syn.final$Metabolite
merged.names.syn.final$row.names <- merged.names.syn.final$HMDB.ID.x
merged.names.syn.f = subset(merged.names.syn.final, select = c("row.names","hmdb_id","class_broad","metabolite_name","method","MetabID") )

# Format merged.names.id.final
merged.names.id.final$metabolite_name <- merged.names.id.final$metabolite_name.x
merged.names.id.final$row.names <- merged.names.id.final$HMDB.ID
merged.names.id.f = subset(merged.names.id.final, select = c("row.names","hmdb_id","class_broad","metabolite_name","method","MetabID") )

# Format merged.names
merged.names$row.names <- merged.names$HMDB.ID.x
merged.names.f <- subset(merged.names, select = c("row.names","hmdb_id","class_broad","metabolite_name","method","MetabID") )

# Merge datasets to create metabolite fdata for WHI - 272 metabolites in total matched among those with <10% missing data
whi.fdata <- rbind(merged.names.id.f, merged.names.f, merged.names.syn.f)
dim(whi.fdata)[1]


#### 3 - CREATE MERGED WHI DATASET WITH COVARIATES AND METABOLITES ####

# Subset to metabolites overlapping with NHS
ids <- whi.fdata$MetabID  
mx.dat.final <- mx.dat[ids]
rownames(mx.dat.final) <- mx.dat$WHI.ID

# Create final pdata dataset
whi.pdata.1 <- cbind(oth.dat2, new.dat2[,4:length(new.dat2)])
whi.pdata <- cbind(whi.pdata.1, draw2$FASTHRS)

# Create final dataset of metabolite values
whi.met.vals <- mx.dat.final

#### RECODE PDATA VARIABLES ####

# Race
table(whi.pdata$ETHNIC)
whi.pdata <- mutate(whi.pdata, race_bw = factor(ETHNIC),
                    race_bw = fct_recode(race_bw,"1" = "Black or African-American","0" = "White (not of Hispanic origin)"),
                    race_bw = fct_relevel(race_bw,"0","1"))
table(whi.pdata$race_bw, useNA="always")

# Age
summary(whi.pdata$AGE, useNA="always")
whi.pdata$age <- whi.pdata$AGE
summary(whi.pdata$age)

# Hormone therapy - ever use
# Note that all women in WHI were postmenopausal
whi.pdata$ht <- rep(0, dim(whi.pdata)[1])
whi.pdata$ht[whi.pdata$estprob == 1 | whi.pdata$testosteroneb == 1] <- 1
table(whi.pdata$ht, useNA="always")

# Fast <8 vs. >=8 hours
summary(whi.pdata$'draw2$FASTHRS')
whi.pdata$fast <- ifelse(whi.pdata$'draw2$FASTHRS'<8, 0, ifelse(whi.pdata$'draw2$FASTHRS'>=8, 1, NA))
table(whi.pdata$fast, useNA="always")

# Smoking status
whi.pdata$smoking <- ifelse(whi.pdata$SMOKING=="", NA, whi.pdata$SMOKING)
table(whi.pdata$smoking, useNA="always")

# Cigarettes per day
whi.pdata$cigsday <- ifelse(whi.pdata$CIGSDAY=="", NA, whi.pdata$CIGSDAY)
table(whi.pdata$cigsday, useNA="always")
table(whi.pdata$cigsday, whi.pdata$smoking)
# Set cigarettes per day to 0 if never smoker (past smokers can have values other than 0 for this variable)
whi.pdata$cigsday_new <- ifelse(is.na(whi.pdata$cigsday)==TRUE & whi.pdata$smoking=="Never Smoked",
                                "0", whi.pdata$cigsday)
table(whi.pdata$cigsday_new, useNA="always")

# Alcohol
whi.pdata$alcohol <- whi.pdata$ALCSWK
summary(whi.pdata$alcohol)

# AHEI diet score - 2005 version
whi.pdata$ahei <- whi.pdata$HEI2005
summary(whi.pdata$ahei)

# Physical activity
whi.pdata$phyact <- whi.pdata$TEXPWK
summary(whi.pdata$phyact)

# Oophorectomy is not available in WHI

# Hysterectomy 
whi.pdata$hyst <- whi.pdata$HYST
table(whi.pdata$hyst, useNA="always")

# Census tract variables are not available but we will include individual income and education

# Income
whi.pdata$income <- ifelse((whi.pdata$INCOME=="" | whi.pdata$INCOME=="Don't know"), NA, whi.pdata$INCOME)
table(whi.pdata$income, useNA="always")

# Education
whi.pdata$educ <- ifelse(whi.pdata$EDUC=="", NA, whi.pdata$EDUC)
table(whi.pdata$educ, useNA="always")

# Neither weight change nor weight at age 18 are available

# Hypertension ever + treatment
table(whi.pdata$HYPT,whi.pdata$HTNTRT, useNA="always")
table(whi.pdata$HYPT,whi.pdata$hypertensionb, useNA="always")
whi.pdata$hypt_trt <- ifelse((whi.pdata$HTNTRT=="Treated hypertensive" | whi.pdata$hypertensionb==1), 2,
                             ifelse(whi.pdata$HYPT=="Yes", 1,
                                    ifelse(whi.pdata$HYPT=="No", 0, NA)))
table(whi.pdata$hypt_trt, useNA="always")

# Diabetes ever + treatment
table(whi.pdata$DIAB,whi.pdata$DIABTRT, useNA="always")
table(whi.pdata$DIAB,whi.pdata$diabetesb, useNA="always")
whi.pdata$diab_trt <- ifelse((whi.pdata$DIABTRT=="Yes" | whi.pdata$diabetesb==1), 2,
                             ifelse(whi.pdata$DIAB=="Yes", 1,
                                    ifelse(whi.pdata$DIAB=="No", 0, NA)))
table(whi.pdata$diab_trt, useNA="always")


# Hyperlipidemia treatment (ever diagnosed is not available): combine 2 variables (STATIN, hilipidb)
statin.ind <- which(colnames(whi.pdata) %in% "statinb")
statin <- whi.pdata[,statin.ind] 
 
hilipidb <- whi.pdata$hilipidb
table(statin, hilipidb)
 
whi.pdata$hyperlipid <- rep(0, length(hilipidb))
whi.pdata$hyperlipid[statin == 1 | hilipidb == 1] <- 1
table(whi.pdata$hyperlipid, useNA="always")


# BMI
whi.pdata$bmi <- whi.pdata$BMI
summary(whi.pdata$bmi)

# Age at menarche, age at first birth / parity, and breastfeeding history will not be included (not required for main analyses)

# Family history of breast cancer and MI will not be included (not required for main analyses)

# Region
whi.pdata$region <- whi.pdata$REGION
table(whi.pdata$region, useNA="always")


# No data is available prior to baseline (e.g., AHEI score in the qx cycle before baseline)


#### 3 - IMPUTE AND TRANSFORM WHI METABOLITES ####

# Load packages and set seed
# install.packages("VIM")
# remotes::install_github("statistikat/VIM")
library(VIM)
library(laeken)
set.seed(02586010)

# Do natural log + z-score transformation of metabolites
# Need to scale all metabolites to same scale before imputation, since KNN works better when data is normalized/scaled
whi.met.vals.ln <- log(whi.met.vals)
whi.met.vals.2 <-scale(whi.met.vals.ln, center = TRUE, scale = TRUE)

# Check histogram
hist(whi.met.vals.2[,3])

# Check for missing data
anyNA(whi.met.vals.2)

# Add demographic variables for imputation
# Subset to variables to be used in imputation algorithm and covariates to be imputed
imput.var = c("WHI.ID","race_bw",
              "age","ht","fast",
              "cigsday_new","alcohol","ahei","phyact",
              "hyst",
              "income","educ",
              "hypt_trt","hyperlipid","diab_trt",
              "smoking",
              "bmi",
              "region") 
whi.pdata.subset <- whi.pdata[,imput.var]

# Check missing values of covariates
map(whi.pdata.subset, ~mean(is.na(.))) # very small % missing

# Combine covariate and metabolite data
dim(whi.met.vals.2)
met.vals.t.rf <- cbind(whi.pdata.subset, whi.met.vals.2)


# Convert all categorical variables to factors
cat.vars = c("race_bw",
             "ht","fast",
             "cigsday_new",
             "hyst",
             "income","educ",
             "hypt_trt","hyperlipid","diab_trt",
             "smoking",
             "region")
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


# Replace met.vals with imputed and transformed data
whi.met.vals <- met.vals.t.imputed[, (length(imput.var)+1):length(met.vals.t.imputed)]

# Set rownames to fdata row.names that matched to NHS
colnames(whi.met.vals) <- whi.fdata$row.names

# Replace pdata with imputed data
whi.pdata <- met.vals.t.imputed[, 1:length(imput.var)]

# Create age-squared variable
whi.pdata$age_sqr <- whi.pdata$age^2

