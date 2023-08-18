##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_3.R

# Programmer: Emma McGee 

# Date Last Updated: February 8, 2023

# Purpose of Program: Estimate residual differences in metabolomic profiles in a hypothetical counterfactual population
                      # in which distributions of risk factors in the Black population are set equal to the distributions in the White
                      # population using a counterfactual analysis based on causal mediation approaches (Supplemental Table 3)

# Statistical Analyses:
    # Inverse odds ratio weighting based analysis for multiple mediators implemented using the CMAverse package
    # Residual differences are estimated, along with model-based observed differences
    # Confidence intervals and p-values are estimated by bootstrap sampling (n=1000 bootstrap samples)

# Additional Study Information: See all_data.R

##########################################################################


#### 1 - PREPARE DATA ####

rm(list = ls())
main.dir = ""
results.dir = file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
# devtools::install_github("BS1125/CMAverse")
library(CMAverse)
library(tidyverse)
library(Biobase)
library(geepack)
library(stats)
library(openxlsx)

# Load eset of known metabolites
load(file.path(main.dir, "data/known.eset.final.RData"))

# Extract pdata and metabolite data
pdata <- pData(known.eset.final)
pdata$race_bw_f <- as.factor(pdata$race_bw)
summary(pdata$race_bw_f)
met.vals <- exprs(known.eset.final)

# Convert hypercholestolemia, hypertension, and diabetes to binary (due to small numbers)
pdata$knn_nhs1_hypthx <- ifelse(pdata$knn_nhs1_hypthx_med=="0", "0", 
                                ifelse(pdata$knn_nhs1_hypthx_med=="1", "1",
                                       ifelse(pdata$knn_nhs1_hypthx_med=="2", "1",NA)))
table(pdata$knn_nhs1_hypthx)
pdata$knn_nhs1_dbhx <- ifelse(pdata$knn_nhs1_dbhx_med=="0", "0", 
                                ifelse(pdata$knn_nhs1_dbhx_med=="1", "1",
                                       ifelse(pdata$knn_nhs1_dbhx_med=="2", "1",NA)))
table(pdata$knn_nhs1_dbhx)
pdata$knn_nhs1_cholhx <- ifelse(pdata$knn_nhs1_cholhx_med=="0", "0", 
                              ifelse(pdata$knn_nhs1_cholhx_med=="1", "1",
                                     ifelse(pdata$knn_nhs1_cholhx_med=="2", "1",NA)))
table(pdata$knn_nhs1_cholhx)

# Restrict to metabolites with 0% missing data after imputation (<10% missing in raw data)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
met.vals <- met.vals[which(missing == 0),]
dim(met.vals)
#met.vals <- met.vals[1:2,]

# Prepare regression output files
nde.model.iorw <- matrix(0,ncol = 4,nrow = dim(met.vals)[1]) #NDE
colnames(nde.model.iorw) <- c("effect.estimate","L_CI","U_CI","p.value")
rownames(nde.model.iorw) <- rownames(met.vals)
nie.model.iorw <- nde.model.iorw #NIE
te.model.iorw <- nde.model.iorw #TE
pm.model.iorw <- nde.model.iorw #PM

#### 2 - RUN MEDIATION ANALYSIS ####

# The inverse odds ratio weighting approach is used because it is agnostic with regards to exposure-mediator interactions and allows for mediator-mediator interactions

# Run all mediation analyses
for(i in 1:dim(met.vals)[1]){
  curr.data <- cbind(pdata[,c("race_bw_f",
                             "id",
                             "knn_nhs1_agebld","knn_nhs1_agebld_sqr","knn_nhs1_menopmh_bld","knn_nhs1_fast","knn_nhs1_bldmonth",
                                              "knn_nhs1_bldmonth_sqr","knn_nhs1_btimelr",
                             "knn_nhs1_smk_per_day","knn_nhs1_alc","knn_nhs1_alc_sqr","knn_nhs1_ahei2010_noALC","knn_nhs1_ahei2010_noALC_sqr",
                                              "knn_nhs1_act","knn_nhs1_act_sqr",
                             "knn_nhs1_oophhx","knn_nhs1_hysthx",
                             "knn_nhs1_ln_income_cens","knn_nhs1_percent_college_cens","knn_nhs1_percent_white_cens",
                             "knn_nhs1_wtchg","knn_nhs1_wtchg_sqr",
                             "knn_nhs1_menopmh_88_der",
                             "knn_nhs1_hypthx","knn_nhs1_dbhx","knn_nhs1_cholhx",
                             "knn_nhs1_percent_pov_cens_76","knn_nhs1_percent_unem_cens_76","knn_nhs1_percent_college_cens_76","knn_nhs1_percent_white_cens_76",
                             "knn_nhs1_bmi18","knn_nhs1_bmi18_sqr"
                             )],met = met.vals[i,])
  
   ### Mediation Model Using Inverse Odds Ratio Weighting: all mediators (regardless of how well-defined the intervention is)
  
   # Model for the outcome (conditional on A and mediator-outcome confounders)
   yreg1 <- glm(met ~ race_bw_f  # A
               # mediator-outcome confounders
               + knn_nhs1_agebld + knn_nhs1_agebld_sqr 
               + knn_nhs1_percent_pov_cens_76 + knn_nhs1_percent_unem_cens_76 + knn_nhs1_percent_college_cens_76 + knn_nhs1_percent_white_cens_76
               + knn_nhs1_bmi18,
               data=curr.data, family=gaussian)
  
   # Model for A (conditional on mediators and mediator-outcome confounders)
   ereg1 <- glm(race_bw_f ~
                 # mediators
                 knn_nhs1_smk_per_day + knn_nhs1_alc + knn_nhs1_ahei2010_noALC + knn_nhs1_act 
               + knn_nhs1_ln_income_cens + knn_nhs1_percent_college_cens + knn_nhs1_percent_white_cens
               + knn_nhs1_hysthx + knn_nhs1_oophhx
               + knn_nhs1_menopmh_88_der
               + knn_nhs1_wtchg
               + knn_nhs1_hypthx + knn_nhs1_dbhx + knn_nhs1_cholhx
               # mediator-outcome confounders
               + knn_nhs1_agebld + knn_nhs1_agebld_sqr 
               + knn_nhs1_percent_pov_cens_76 + knn_nhs1_percent_unem_cens_76 + knn_nhs1_percent_college_cens_76 + knn_nhs1_percent_white_cens_76
               + knn_nhs1_bmi18,         
                 data=curr.data, family=binomial(link="logit"))
  
   # Mediation analysis
   my.med.model.iorw <- cmest(data = curr.data, model = "iorw",
                         outcome = "met", exposure = "race_bw_f",
                         mediator = c("knn_nhs1_smk_per_day",
                                      "knn_nhs1_alc",
                                      "knn_nhs1_ahei2010_noALC",
                                      "knn_nhs1_act",
                                      "knn_nhs1_ln_income_cens",
                                      "knn_nhs1_percent_college_cens",
                                      "knn_nhs1_percent_white_cens",
                                      "knn_nhs1_hysthx","knn_nhs1_oophhx",
                                      "knn_nhs1_menopmh_88_der",
                                      "knn_nhs1_wtchg",
                                      "knn_nhs1_hypthx","knn_nhs1_dbhx","knn_nhs1_cholhx"),
                         basec = c("knn_nhs1_agebld","knn_nhs1_agebld_sqr",
                                   "knn_nhs1_percent_pov_cens_76","knn_nhs1_percent_unem_cens_76","knn_nhs1_percent_college_cens_76","knn_nhs1_percent_white_cens_76",
                                   "knn_nhs1_bmi18"),
                         yreg = yreg1,
                         ereg = ereg1,
                         astar = 0, a = 1,
                         estimation = "imputation",
                         inference = "bootstrap",
                         nboot = 1000)
   my.model.model.iorw <- summary(my.med.model.iorw)  

  # Export results for IORW approach: all mediators
  nde.model.iorw[i,] <- c(my.model.model.iorw$effect.pe["pnde"],my.model.model.iorw$effect.ci.low["pnde"],my.model.model.iorw$effect.ci.high["pnde"],my.model.model.iorw$effect.pval["pnde"])
  nie.model.iorw[i,] <- c(my.model.model.iorw$effect.pe["tnie"],my.model.model.iorw$effect.ci.low["tnie"],my.model.model.iorw$effect.ci.high["tnie"],my.model.model.iorw$effect.pval["tnie"])
  te.model.iorw[i,] <- c(my.model.model.iorw$effect.pe["te"],my.model.model.iorw$effect.ci.low["te"],my.model.model.iorw$effect.ci.high["te"],my.model.model.iorw$effect.pval["te"])
  pm.model.iorw[i,] <- c(my.model.model.iorw$effect.pe["pm"],my.model.model.iorw$effect.ci.low["pm"],my.model.model.iorw$effect.ci.high["pm"],my.model.model.iorw$effect.pval["pm"])
}



#### 3 - SAVE ALL RESULTS ####

# Extract metabolite names
fdata <- fData(known.eset.final)
myvars <- c("class_broad","hmdb_id","metabolite_name")
metab.name <- fdata[myvars]
# Recode PC plasmalogen 
metab.name$class_broad <- ifelse(metab.name$class_broad=="PC plasmalogens", "Phosphatidylcholine plasmalogens", metab.name$class_broad)


# Format results for printing
list.results <- list(nde.model.iorw, nie.model.iorw, te.model.iorw, pm.model.iorw)
list.results <- lapply(list.results, as.data.frame)
list.results.format <- lapply(list.results, function(x) cbind(format(round(x[c("effect.estimate","L_CI","U_CI")],2),nsmall=2),
                                                              "p.value"=ifelse(round(x$p.value,2)>0.99, ">0.99",
                                                                               ifelse(x$p.value>=0.01, format(round(x$p.value,2),nsmall=2),
                                                                                      ifelse((x$p.value>=0.001 & x$p.value<0.01), format(round(x$p.value,3),nsmall=3),
                                                                                             ifelse(x$p.value==0, "<0.001",
                                                                                                  format(x$p.value,digits=3)))))))

list.results.ci <- lapply(list.results.format, function(x) cbind(x, "conf.int"=paste("(", x$L_CI, ", ", x$U_CI, ")",sep = "")))
final.results <- lapply(list.results.ci, function(x) x[,c("effect.estimate","conf.int","p.value")])

# Merge IORW-based results: all mediators
mediation_results_iorw <- (cbind(final.results[[1]],final.results[[2]],final.results[[3]],final.results[[4]]))
colnames(mediation_results_iorw)[1:3] <- paste(colnames(mediation_results_iorw)[1:3], "nde", sep = ".")
colnames(mediation_results_iorw)[4:6] <- paste(colnames(mediation_results_iorw)[4:6], "nie", sep = ".")
colnames(mediation_results_iorw)[7:9] <- paste(colnames(mediation_results_iorw)[7:9], "te", sep = ".")
colnames(mediation_results_iorw)[10:12] <- paste(colnames(mediation_results_iorw)[10:12], "pm", sep = ".")
row.names(mediation_results_iorw) <- rownames(met.vals)

# Save all results
write.xlsx(mediation_results_iorw, file = file.path(main.dir,"data/mediation_results_iorw_based_all_mediators.xlsx"),
           overwrite=TRUE)

#### 4 - CREATE SUPPLEMENTAL TABLE ####

# Export TE and NDE, the residual difference measure, using the IORW-based approach with all mediators
results_nde <- cbind(final.results[[3]],final.results[1])
colnames(results_nde)[1:3] <- paste(colnames(results_nde)[1:3], "te", sep = ".")
colnames(results_nde)[4:6] <- paste(colnames(results_nde)[4:6], "nde", sep = ".")

# Merge metabolite names with results
row.names(results_nde) <- rownames(met.vals)
results.final <- merge(metab.name, results_nde, by="row.names")
results.final <- results.final[ , !(names(results.final) %in% c("Row.names","Row.names.y"))]

# Arrange by metabolite class and hmdb_id
results.final <- arrange(results.final, class_broad, hmdb_id)

# Format missing HMDB ID and class information
results.final$hmdb_id <- ifelse(is.na(results.final$hmdb_id) == TRUE, "Unavailable", results.final$hmdb_id)
results.final$class_broad <- ifelse(is.na(results.final$class_broad) == TRUE, "Unclassified", results.final$class_broad)

# Save table
write.xlsx(results.final, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_3.xlsx"),
           overwrite=TRUE)
