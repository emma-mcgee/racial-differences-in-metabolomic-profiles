##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_10.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: observed differences in individual known metabolites with >=10% missingness, logistic regression for metabolite presence vs. absence (Supplemental Table 10)

# Statistical Analyses:
    # Logistic regression with metabolite presence vs. absence as the outcome variable

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####
  
rm(list = ls())
main.dir <- ""
results.dir <- file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(tidyverse)
library(Biobase)
library(openxlsx)

# Load eset of known metabolites
load(file.path(main.dir, "data/known.eset.final.RData"))

# Prepare data
pdata <- pData(known.eset.final)
met.vals <- exprs(known.eset.final)

# Restrict to metabolites with >=10% missing data
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
met.vals <- met.vals[which(missing >= 10),]

# Create present vs. absent indicator
met.miss <- ifelse(is.na(met.vals)==TRUE,0,1)

# Prepare regression output files
lin.res.model0 <- matrix(0,ncol = 5,nrow = dim(met.miss)[1])
colnames(lin.res.model0) <- c("effect.estimate","standard.error","p.value","L_CI","U_CI")
rownames(lin.res.model0) <- rownames(met.miss)

#### 2 - QUANTIFY OBSERVED DIFFERENCES USING LINEAR REGRESSION ####

# Run all linear regressions
for(i in 1:dim(met.miss)[1]){
  curr.data <- cbind(pdata[,c("race_bw",
                             "knn_nhs1_agebld","knn_nhs1_agebld_sqr")],met = met.miss[i,])
  
  # NHS: adjusted for age
  my.lm.model0 <- glm(met~race_bw+
                      knn_nhs1_agebld+knn_nhs1_agebld_sqr,
                    data = curr.data,family = binomial(link = "logit"))
  my.model.model0 <- summary(my.lm.model0)
  coefficient <- my.model.model0$coefficients
  Estimate <- exp(coefficient["race_bw1",c("Estimate")])
  my.ci.model0 <- exp(confint(my.lm.model0))
  
  lin.res.model0[i,] <- c(Estimate,my.model.model0$coefficients["race_bw1",c("Std. Error","Pr(>|z|)")],my.ci.model0["race_bw1",])
}

#### 3 - CREATE TABLE ####

# Prepare data for table

# Extract metabolite names
fdata <- fData(known.eset.final)
myvars <- c("class_broad","hmdb_id","metabolite_name")
metab.name <- fdata[myvars]

# Format results for printing
list.results <- list(lin.res.model0)
list.results <- lapply(list.results, as.data.frame)
list.results.format <- lapply(list.results, function(x) cbind(format(round(x[c("effect.estimate","L_CI","U_CI")],2),nsmall=2),
                                                              "p.value"=ifelse(round(x$p.value,2)>0.99, ">0.99",
                                                                               ifelse(x$p.value>=0.01, format(round(x$p.value,2),nsmall=2),
                                                                                      ifelse((x$p.value>=0.001 & x$p.value<0.01), format(round(x$p.value,3),nsmall=3),
                                                                                             format(x$p.value,digits=3))))))
#list.results.format <- lapply(list.results, function(x) cbind(format(round(x[c("effect.estimate","L_CI","U_CI")],2),nsmall=2),
#                                                              "p.value"=ifelse(round(x$p.value,3)>0.000, format(round(x$p.value,3),nsmall=3),
#                                                                               format(x$p.value,digits=3))))
list.results.ci <- lapply(list.results.format, function(x) cbind(x, "conf.int"=paste("(", x$L_CI, ", ", x$U_CI, ")",sep = "")))
final.results <- lapply(list.results.ci, function(x) x[,c("effect.estimate","conf.int","p.value")])

# Merge metabolite names with results
results.final <- merge(metab.name, final.results[[1]], by="row.names")
row.names(results.final) <- results.final$Row.names
results.final <- results.final[ , !(names(results.final) %in% c("Row.names"))]

# Add number of metabolites present and absent
met.miss.dat <- as.data.frame(t(met.miss))
n <- matrix(0,ncol = 2,nrow = dim(met.vals)[1])
colnames(n) <- c("n_present","n_absent")
rownames(n) <- rownames(met.vals)
for(i in 1:dim(met.vals)[1]){
  present <- table(met.miss.dat[,i])[2]
  absent <- table(met.miss.dat[,i])[1]
  n[i,] = c(present, absent)
}

# Merge results
results.final.n <- merge(n, results.final, by="row.names")

# Order results
results.final.n.sub <- results.final.n[,c("class_broad","hmdb_id","metabolite_name","n_present","n_absent","effect.estimate","conf.int","p.value")]

# Arrange by metabolite class and hmdb_id
results.final.n.sub <- arrange(results.final.n.sub, class_broad, hmdb_id)

# Format missing HMDB ID and class information
results.final.n.sub$hmdb_id <- ifelse(is.na(results.final.n.sub$hmdb_id) == TRUE, "Unavailable", results.final.n.sub$hmdb_id)
results.final.n.sub$class_broad <- ifelse(is.na(results.final.n.sub$class_broad) == TRUE, "Unclassified", results.final.n.sub$class_broad)

# Save table
write.xlsx(results.final.n.sub, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_10.xlsx"), overwrite=TRUE)
