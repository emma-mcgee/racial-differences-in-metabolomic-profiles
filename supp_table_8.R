##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_8.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: quantify observed differences in individual known metabolites, stratified by birth cohort (Supplemental Table 8)

# Statistical Analyses:
    # Linear regression with each ln Z-score transformed metabolite as the outcome variable
    # Stratified by birth cohort in decades (20s, 30s, 40s)

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

# Check birth cohort data
table(known.eset.final$nhs1_birth_cohort, useNA="always")
table(known.eset.final$nhs1_birth_cohort, known.eset.final$race_bw)

# Prepare data
pdata <- pData(known.eset.final)
met.vals <- exprs(known.eset.final)

# Restrict to metabolites with 0% missing data after imputation (<10% missing in raw data)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
met.vals <- met.vals[which(missing == 0),]
dim(met.vals)

# Prepare regression output files
lin.res.model0 <- matrix(0,ncol = 3,nrow = dim(met.vals)[1])
colnames(lin.res.model0) <- c("effect.estimate","L_CI","U_CI")
rownames(lin.res.model0) <- rownames(met.vals)
lin.res.model.20 <- lin.res.model0
lin.res.model.30 <- lin.res.model0
lin.res.model.40 <- lin.res.model0
lin.res.model.int.30 <- lin.res.model0
lin.res.model.int.40 <- lin.res.model0
lin.res.model.30w <- lin.res.model0
lin.res.model.40w <- lin.res.model0
lin.res.model.20b <- lin.res.model0
lin.res.model.30b <- lin.res.model0
lin.res.model.40b <- lin.res.model0

p.het <- matrix(0,ncol = 1,nrow = dim(met.vals)[1])
colnames(p.het) <- c("p.het")
rownames(p.het) <- rownames(met.vals)

#### 2 - QUANTIFY OBSERVED DIFFERENCES USING LINEAR REGRESSION ####

# Run linear regressions to estimate observed differences in metabolite levels
for(i in 1:dim(met.vals)[1]){
  curr.data <- cbind(pdata[,c("race_bw",
                             "nhs1_birth_cohort",
                             "knn_nhs1_agebld","knn_nhs1_agebld_sqr")],met = met.vals[i,])
  
  # Model with no product term, all women
  my.lm.model0 <- lm(met ~ race_bw + nhs1_birth_cohort +
                      knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                    data = curr.data)
  my.model.model0 <- summary(my.lm.model0)      
  my.ci.model0 <- confint(my.lm.model0)
  
  # Model with product term, all women
  my.lm.model.int <- lm(met ~ race_bw + nhs1_birth_cohort + nhs1_birth_cohort*race_bw +
                      knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                    data = curr.data)
  my.model.model.int <- summary(my.lm.model.int)      
  my.ci.model.int <- confint(my.lm.model.int)
  vcov.model.int <- vcov(my.lm.model.int)
  
  # Calculating estimates and CIs not provided directly by the model
  # Estimate - Black woman born in the 30s vs. White woman born in the 20s
  estimate30b <- (my.model.model.int$coefficients["race_bw1",c("Estimate")] + my.model.model.int$coefficients["nhs1_birth_cohort30s",c("Estimate")] + my.model.model.int$coefficients["race_bw1:nhs1_birth_cohort30s",c("Estimate")])
  # CI
  var30b <- (vcov.model.int["race_bw1","race_bw1"]+ vcov.model.int["nhs1_birth_cohort30s","nhs1_birth_cohort30s"]+ vcov.model.int["race_bw1:nhs1_birth_cohort30s","race_bw1:nhs1_birth_cohort30s"]+ 2*vcov.model.int["race_bw1","nhs1_birth_cohort30s"]+ 2*vcov.model.int["race_bw1","race_bw1:nhs1_birth_cohort30s"]+ 2*vcov.model.int["nhs1_birth_cohort30s","race_bw1:nhs1_birth_cohort30s"])
  l_ci30b <- estimate30b - qnorm(0.975)*sqrt(var30b)
  u_ci30b <- estimate30b + qnorm(0.975)*sqrt(var30b)
  
  # Estimate - Black woman born in the 40s vs. White woman born in the 20s
  estimate40b <- (my.model.model.int$coefficients["race_bw1",c("Estimate")] + my.model.model.int$coefficients["nhs1_birth_cohort40s",c("Estimate")]+ my.model.model.int$coefficients["race_bw1:nhs1_birth_cohort40s",c("Estimate")])
  # CI
  var40b <- (vcov.model.int["race_bw1","race_bw1"]+ vcov.model.int["nhs1_birth_cohort40s","nhs1_birth_cohort40s"]+ vcov.model.int["race_bw1:nhs1_birth_cohort40s","race_bw1:nhs1_birth_cohort40s"]+ 2*vcov.model.int["race_bw1","nhs1_birth_cohort40s"]+ 2*vcov.model.int["race_bw1","race_bw1:nhs1_birth_cohort40s"]+ 2*vcov.model.int["nhs1_birth_cohort40s","race_bw1:nhs1_birth_cohort40s"])
  l_ci40b <- estimate40b - qnorm(0.975)*sqrt(var40b)
  u_ci40b <- estimate40b + qnorm(0.975)*sqrt(var40b)
  
  # ANOVA F-test
  anova.results <- anova(my.lm.model0, my.lm.model.int)
  
  # Model for Black-White mean difference among those born in the 20s
  my.lm.model20 <- lm(met ~ race_bw +
                       knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                     data = curr.data, subset = nhs1_birth_cohort=="20s")
  my.model.model20 <- summary(my.lm.model20)      
  my.ci.model20 <- confint(my.lm.model20)
  
  # Model for Black-White mean difference among those born in the 30s
  my.lm.model30 <- lm(met ~ race_bw +
                        knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                      data = curr.data, subset = nhs1_birth_cohort=="30s")
  my.model.model30 <- summary(my.lm.model30)      
  my.ci.model30 <- confint(my.lm.model30)
  
  # Model for Black-White mean difference among those born in the 40s
  my.lm.model40 <- lm(met ~ race_bw +
                        knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                      data = curr.data, subset = nhs1_birth_cohort=="40s")
  my.model.model40 <- summary(my.lm.model40)      
  my.ci.model40 <- confint(my.lm.model40)
  
  lin.res.model.30w[i,] <- c(my.model.model.int$coefficients["nhs1_birth_cohort30s",c("Estimate")],my.ci.model.int["nhs1_birth_cohort30s",])
  lin.res.model.40w[i,] <- c(my.model.model.int$coefficients["nhs1_birth_cohort40s",c("Estimate")],my.ci.model.int["nhs1_birth_cohort40s",])
  lin.res.model.20b[i,] <- c(my.model.model.int$coefficients["race_bw1",c("Estimate")],my.ci.model.int["race_bw1",])
  lin.res.model.30b[i,] <- c(estimate30b, l_ci30b, u_ci30b)
  lin.res.model.40b[i,] <- c(estimate40b, l_ci40b, u_ci40b)
  lin.res.model.20[i,] <- c(my.model.model20$coefficients["race_bw1",c("Estimate")],my.ci.model20["race_bw1",])
  lin.res.model.30[i,] <- c(my.model.model30$coefficients["race_bw1",c("Estimate")],my.ci.model30["race_bw1",])
  lin.res.model.40[i,] <- c(my.model.model40$coefficients["race_bw1",c("Estimate")],my.ci.model40["race_bw1",])
  lin.res.model.int.30[i,] <- c(my.model.model.int$coefficients["race_bw1:nhs1_birth_cohort30s",c("Estimate")],my.ci.model.int["race_bw1:nhs1_birth_cohort30s",])
  lin.res.model.int.40[i,] <- c(my.model.model.int$coefficients["race_bw1:nhs1_birth_cohort40s",c("Estimate")],my.ci.model.int["race_bw1:nhs1_birth_cohort40s",])
  p.het[i,] <- c(anova.results["Pr(>F)"][2,])
}


#### 3 - CREATE TABLE ####

# Extract metabolite names
fdata <- fData(known.eset.final)
myvars <- c("class_broad","hmdb_id","metabolite_name")
fdata <- fdata[myvars]
# Recode PC plasmalogen 
fdata$class_broad <- ifelse(fdata$class_broad=="PC plasmalogens", "Phosphatidylcholine plasmalogens", fdata$class_broad)


# Format results for printing
list.results <- list(lin.res.model.30w, lin.res.model.40w, lin.res.model.20b, lin.res.model.30b, lin.res.model.40b,
                     lin.res.model.20, lin.res.model.30, lin.res.model.40, lin.res.model.int.30, lin.res.model.int.40)
list.results <- lapply(list.results, as.data.frame)
list.results.format <- lapply(list.results, function(x) cbind(format(round(x[c("effect.estimate","L_CI","U_CI")],2),nsmall=2)))
final.results <- lapply(list.results.format, function(x) cbind(x, "col"=paste(x$effect.estimate," (", x$L_CI, ", ", x$U_CI, ")",sep = "")))
p.val <- as.data.frame(p.het)
p.val$p.het <- ifelse(round(p.val$p.het,2)>0.99, ">0.99",
                 ifelse(p.val$p.het>=0.01, format(round(p.val$p.het,2),nsmall=2),
                        ifelse((p.val$p.het>=0.001 & p.val$p.het<0.01), format(round(p.val$p.het,3),nsmall=3),
                               format(p.val$p.het,digits=3))))

# Combine results - first row for each metabolite
results <- merge(final.results[[3]][c("col")], final.results[[6]][c("col")], by="row.names", suffixes = c(".2",".3"))
col.1 <- "Ref."
col.4 <- "Ref."
col.5 <- ""
results.all.cols <- cbind(col.1,results,col.4,col.5)
results.all.cols <- results.all.cols[c("Row.names","col.1","col.2","col.3","col.4","col.5")]

# Merge metabolite info with results
row.names(results.all.cols) <- results.all.cols$Row.names
results.final.cols <- merge(fdata, results.all.cols, by="row.names")
results.final.cols <- results.final.cols[ , !(names(results.final.cols) %in% c("Row.names","Row.names.y"))]

# Add numeric identifiers for metabolites with missing HMDB IDs
for(i in 1:dim(results.all.cols)[1]){
  results.final.cols[i,] <- c(results.final.cols[i,1],
                                ifelse(is.na(results.final.cols$hmdb_id[i])==TRUE, paste("Unavailable", i, sep=""), results.final.cols$hmdb_id[i]),
                                results.final.cols[i,3:8])
}
results.final.cols.2 <- results.final.cols

# Add second row for each metabolite
for(i in 1:dim(results.all.cols)[1]){
  results.final.cols.2 <- results.final.cols.2 %>% add_row(class_broad = results.final.cols$class_broad[i],
                                                   hmdb_id = paste(results.final.cols$hmdb_id[i],"_2"),
                                                   metabolite_name = "",
                                                   col.1 = final.results[[1]][i,c("col")],
                                                   col.2 = final.results[[4]][i,c("col")],
                                                   col.3 = final.results[[7]][i,c("col")],
                                                   col.4 = final.results[[9]][i,c("col")],
                                                   col.5 = "",
                                                   .before = i)
}
results.final.cols.3 <- results.final.cols.2

# Add third row for each metabolite
for(i in 1:dim(results.all.cols)[1]){
  results.final.cols.3 <- results.final.cols.3 %>% add_row(class_broad = results.final.cols$class_broad[i],
                                                     hmdb_id = paste(results.final.cols$hmdb_id[i],"_3"),
                                                     metabolite_name = "",
                                                     col.1 = final.results[[2]][i,c("col")],
                                                     col.2 = final.results[[5]][i,c("col")],
                                                     col.3 = final.results[[8]][i,c("col")],
                                                     col.4 = final.results[[10]][i,c("col")],
                                                     col.5 = p.val$p.het[i],
                                                     .before = i)
}

# Arrange by metabolite class and hmdb_id
results.final.cols.3 <- arrange(results.final.cols.3, class_broad, hmdb_id)

# Save table
write.xlsx(results.final.cols.3, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_8.xlsx"), overwrite=TRUE)

# Check number of p-values for heterogeneity that are nominally significant at alpha = 0.05
p_cutoff <- 0.05
full.adjust.results <- as.data.frame(p.het)
length(full.adjust.results[ which(full.adjust.results$p.het<p_cutoff), ])
