##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_1.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Characteristics of included metabolites (Supplemental Table 1)

# Statistical Analyses:
    # Descriptive statistics
    # The following statistics are included from the merged metabolomics dataset
            # Coefficients of variation (CVs)
            # Intraclass correlation coefficients (ICCs)
            # Spearman rho

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####
rm(list = ls())
main.dir = ""
results.dir = file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(tidyverse)
library(Biobase)
library(openxlsx)

# Load expression set
load(file.path(main.dir,"data/known.eset.final.RData"))
load(file.path(main.dir, "data/merged.eset.raw.RData"))

# Extract data
fdata <- fData(known.eset.final)

# Create percent missing variable - original unimputed data
table(merged.eset.raw$nhs1_canhx, merged.eset.raw$nhs1_caco, useNA="always")
merged.eset.raw <- merged.eset.raw[ , (is.na(merged.eset.raw$nhs1_race)==FALSE) & (merged.eset.raw$nhs1_race == "1" | merged.eset.raw$nhs1_race == "2")]
dim(merged.eset.raw)
met.vals <- exprs(merged.eset.raw)
met.vals.df <- as.data.frame(met.vals)
met.vals.df$percent_missing <- apply(met.vals.df, 1, function(x) paste(format(round(length(which(is.na(x)))/length(x)*100,1),nsmall=1)))
met.vals.df$percent_missing <- ifelse(met.vals.df$percent_missing == "0.0", "0", met.vals.df$percent_missing)
miss <- met.vals.df[c("percent_missing")]
table.dat <- merge(fdata, miss, by="row.names")

# Susbet to <10% missing (primary analysis - using 11 for subsetting due to rounding)
table.dat <- table.dat[which(as.numeric(table.dat$percent_missing)<11),]
summary(as.numeric(table.dat$percent_missing))

# Subset and format variables for table
keep.vars <- c("class_broad","hmdb_id","metabolite_name","method",
              "percent_missing",
              "mean_cv","mean_icc",
              "wps_pilot_icc_p658","wps_pilot_rho_p658")
final.dat <- table.dat[,keep.vars]
# Recode PC plasmalogen 
final.dat$class_broad <- ifelse(final.dat$class_broad=="PC plasmalogens", "Phosphatidylcholine plasmalogens", final.dat$class_broad)
# 1 digit after the decimal for CVs
format.vars.1 <- c("mean_cv")
final.dat$mean_cv <- as.numeric(final.dat$mean_cv)
final.dat[format.vars.1] <- format(final.dat[format.vars.1], digits=1, nsmall=1)
head(final.dat$mean_cv)
any(is.na(final.dat[format.vars.1]))
# 2 digits after the decimal for all other numeric values
format.vars.2 <- c("mean_icc",
                "wps_pilot_icc_p658","wps_pilot_rho_p658")
final.dat[format.vars.2] <- format(final.dat[format.vars.2], nsmall=2)
final.dat[format.vars.2] <- mutate_all(final.dat[format.vars.2],
                                    str_replace_all, "NA", "NM")


# Arrange by metabolite class and hmdb_id
final.dat <- arrange(final.dat, class_broad, hmdb_id)

# Format missing HMDB ID and class information
final.dat$hmdb_id <- ifelse(is.na(final.dat$hmdb_id) == TRUE, "Unavailable", final.dat$hmdb_id)
final.dat$class_broad <- ifelse(is.na(final.dat$class_broad) == TRUE, "Unclassified", final.dat$class_broad)

# Format methods
final.dat <- mutate(final.dat, method = factor(method),
                    method = fct_recode(method, "HILIC-positive" = "HILIC-pos", "C8-positive" = "C8-pos"))


#### 2 - CREATE TABLE ####

# Export table
write.xlsx(final.dat, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_1.xlsx"), overwrite=TRUE)

