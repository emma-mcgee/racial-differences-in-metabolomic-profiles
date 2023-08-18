##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_9.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: quantify observed differences in individual known metabolites, stratified by geographic region (Supplemental Table 9)

# Statistical Analyses:
    # Linear regression with each ln Z-score transformed metabolite as the outcome variable
    # Stratified by geographic region of residence (20s, 30s, 40s)

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

# Check geographic region data
table(known.eset.final$nhs1_region, useNA="always")
table(known.eset.final$nhs1_region, known.eset.final$race_bw)

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
lin.res.model.northeast <- lin.res.model0
lin.res.model.midwest <- lin.res.model0
lin.res.model.south <- lin.res.model0
lin.res.model.west <- lin.res.model0
lin.res.model.int.midwest <- lin.res.model0
lin.res.model.int.south <- lin.res.model0
lin.res.model.int.west <- lin.res.model0
lin.res.model.midwest.w <- lin.res.model0
lin.res.model.south.w <- lin.res.model0
lin.res.model.west.w <- lin.res.model0
lin.res.model.northeast.b <- lin.res.model0
lin.res.model.midwest.b <- lin.res.model0
lin.res.model.south.b <- lin.res.model0
lin.res.model.west.b <- lin.res.model0

p.het <- matrix(0,ncol = 1,nrow = dim(met.vals)[1])
colnames(p.het) <- c("p.het")
rownames(p.het) <- rownames(met.vals)

#### 2 - QUANTIFY OBSERVED DIFFERENCES USING LINEAR REGRESSION ####

# Run linear regressions to estimate observed differences in metabolite levels
for(i in 1:dim(met.vals)[1]){
  curr.data <- cbind(pdata[,c("race_bw",
                             "nhs1_region",
                             "knn_nhs1_agebld","knn_nhs1_agebld_sqr")],met = met.vals[i,])
  
  # Model with no product term, all women
  my.lm.model0 <- lm(met ~ race_bw + nhs1_region +
                      knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                    data = curr.data)
  my.model.model0 <- summary(my.lm.model0)      
  my.ci.model0 <- confint(my.lm.model0)
  
  # Model with product term, all women
  my.lm.model.int <- lm(met ~ race_bw + nhs1_region + nhs1_region*race_bw +
                      knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                    data = curr.data)
  my.model.model.int <- summary(my.lm.model.int)      
  my.ci.model.int <- confint(my.lm.model.int)
  vcov.model.int <- vcov(my.lm.model.int)
  
  # Calculating estimates and CIs not provided directly by the model
  # Estimate - Black woman living in Midwest vs. White woman living in Northeast
  estimate.midwest.b <- (my.model.model.int$coefficients["race_bw1",c("Estimate")] + my.model.model.int$coefficients["nhs1_region2=Midwest",c("Estimate")] + my.model.model.int$coefficients["race_bw1:nhs1_region2=Midwest",c("Estimate")])
  # CI
  var.midwest.b <- (vcov.model.int["race_bw1","race_bw1"]+ vcov.model.int["nhs1_region2=Midwest","nhs1_region2=Midwest"]+ vcov.model.int["race_bw1:nhs1_region2=Midwest","race_bw1:nhs1_region2=Midwest"]+ 2*vcov.model.int["race_bw1","nhs1_region2=Midwest"]+ 2*vcov.model.int["race_bw1","race_bw1:nhs1_region2=Midwest"]+ 2*vcov.model.int["nhs1_region2=Midwest","race_bw1:nhs1_region2=Midwest"])
  l_ci.midwest.b <- estimate.midwest.b - qnorm(0.975)*sqrt(var.midwest.b)
  u_ci.midwest.b <- estimate.midwest.b + qnorm(0.975)*sqrt(var.midwest.b)
  
  # Estimate - Black woman living in South vs. White woman living in Northeast
  estimate.south.b <- (my.model.model.int$coefficients["race_bw1",c("Estimate")] + my.model.model.int$coefficients["nhs1_region3=South",c("Estimate")]+ my.model.model.int$coefficients["race_bw1:nhs1_region3=South",c("Estimate")])
  # CI
  var.south.b <- (vcov.model.int["race_bw1","race_bw1"]+ vcov.model.int["nhs1_region3=South","nhs1_region3=South"]+ vcov.model.int["race_bw1:nhs1_region3=South","race_bw1:nhs1_region3=South"]+ 2*vcov.model.int["race_bw1","nhs1_region3=South"]+ 2*vcov.model.int["race_bw1","race_bw1:nhs1_region3=South"]+ 2*vcov.model.int["nhs1_region3=South","race_bw1:nhs1_region3=South"])
  l_ci.south.b <- estimate.south.b - qnorm(0.975)*sqrt(var.south.b)
  u_ci.south.b <- estimate.south.b + qnorm(0.975)*sqrt(var.south.b)
  
  # Estimate - Black woman living in West vs. White woman living in Northeast
  estimate.west.b <- (my.model.model.int$coefficients["race_bw1",c("Estimate")] + my.model.model.int$coefficients["nhs1_region3=South",c("Estimate")]+ my.model.model.int$coefficients["race_bw1:nhs1_region3=South",c("Estimate")])
  # CI
  var.west.b <- (vcov.model.int["race_bw1","race_bw1"]+ vcov.model.int["nhs1_region4=West","nhs1_region4=West"]+ vcov.model.int["race_bw1:nhs1_region4=West","race_bw1:nhs1_region4=West"]+ 2*vcov.model.int["race_bw1","nhs1_region4=West"]+ 2*vcov.model.int["race_bw1","race_bw1:nhs1_region4=West"]+ 2*vcov.model.int["nhs1_region4=West","race_bw1:nhs1_region4=West"])
  l_ci.west.b <- estimate.west.b - qnorm(0.975)*sqrt(var.west.b)
  u_ci.west.b <- estimate.west.b + qnorm(0.975)*sqrt(var.west.b)
  
  # ANOVA F-test
  anova.results <- anova(my.lm.model0, my.lm.model.int)
  
  # Model for Black-White mean difference among those living in the Northeast
  my.lm.model.northeast <- lm(met ~ race_bw +
                              knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                            data = curr.data, subset = nhs1_region=="1=Northeast")
  my.model.model.northeast <- summary(my.lm.model.northeast)      
  my.ci.model.northeast <- confint(my.lm.model.northeast)
  
  # Model for Black-White mean difference among those living in the Midwest
  my.lm.model.midwest <- lm(met ~ race_bw +
                       knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                     data = curr.data, subset = nhs1_region=="2=Midwest")
  my.model.model.midwest <- summary(my.lm.model.midwest)      
  my.ci.model.midwest <- confint(my.lm.model.midwest)
  
  # Model for Black-White mean difference among those living in the South
  my.lm.model.south <- lm(met ~ race_bw +
                        knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                      data = curr.data, subset = nhs1_region=="3=South")
  my.model.model.south <- summary(my.lm.model.south)      
  my.ci.model.south <- confint(my.lm.model.south)
  
  # Model for Black-White mean difference among those living in the West
  my.lm.model.west <- lm(met ~ race_bw +
                        knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                      data = curr.data, subset = nhs1_region=="4=West")
  my.model.model.west <- summary(my.lm.model.west)      
  my.ci.model.west <- confint(my.lm.model.west)
  
  lin.res.model.midwest.w[i,] <- c(my.model.model.int$coefficients["nhs1_region2=Midwest",c("Estimate")],my.ci.model.int["nhs1_region2=Midwest",])
  lin.res.model.south.w[i,] <- c(my.model.model.int$coefficients["nhs1_region3=South",c("Estimate")],my.ci.model.int["nhs1_region3=South",])
  lin.res.model.west.w[i,] <- c(my.model.model.int$coefficients["nhs1_region4=West",c("Estimate")],my.ci.model.int["nhs1_region4=West",])
  lin.res.model.northeast.b[i,] <- c(my.model.model.int$coefficients["race_bw1",c("Estimate")],my.ci.model.int["race_bw1",])
  lin.res.model.midwest.b[i,] <- c(estimate.midwest.b, l_ci.midwest.b, u_ci.midwest.b)
  lin.res.model.south.b[i,] <- c(estimate.south.b, l_ci.south.b, u_ci.south.b)
  lin.res.model.west.b[i,] <- c(estimate.west.b, l_ci.west.b, u_ci.west.b)
  lin.res.model.northeast[i,] <- c(my.model.model.northeast$coefficients["race_bw1",c("Estimate")],my.ci.model.northeast["race_bw1",])
  lin.res.model.midwest[i,] <- c(my.model.model.midwest$coefficients["race_bw1",c("Estimate")],my.ci.model.midwest["race_bw1",])
  lin.res.model.south[i,] <- c(my.model.model.south$coefficients["race_bw1",c("Estimate")],my.ci.model.south["race_bw1",])
  lin.res.model.west[i,] <- c(my.model.model.west$coefficients["race_bw1",c("Estimate")],my.ci.model.west["race_bw1",])
  lin.res.model.int.midwest[i,] <- c(my.model.model.int$coefficients["race_bw1:nhs1_region2=Midwest",c("Estimate")],my.ci.model.int["race_bw1:nhs1_region2=Midwest",])
  lin.res.model.int.south[i,] <- c(my.model.model.int$coefficients["race_bw1:nhs1_region3=South",c("Estimate")],my.ci.model.int["race_bw1:nhs1_region3=South",])
  lin.res.model.int.west[i,] <- c(my.model.model.int$coefficients["race_bw1:nhs1_region4=West",c("Estimate")],my.ci.model.int["race_bw1:nhs1_region4=West",])
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
list.results <- list(lin.res.model.midwest.w, lin.res.model.south.w, lin.res.model.west.w,
                     lin.res.model.northeast.b, lin.res.model.midwest.b, lin.res.model.south.b, lin.res.model.west.b,
                     lin.res.model.northeast, lin.res.model.midwest, lin.res.model.south, lin.res.model.west,
                     lin.res.model.int.midwest, lin.res.model.int.south, lin.res.model.int.west)
list.results <- lapply(list.results, as.data.frame)
list.results.format <- lapply(list.results, function(x) cbind(format(round(x[c("effect.estimate","L_CI","U_CI")],2),nsmall=2)))
final.results <- lapply(list.results.format, function(x) cbind(x, "col"=paste(x$effect.estimate," (", x$L_CI, ", ", x$U_CI, ")",sep = "")))
p.val <- as.data.frame(p.het)
p.val$p.het <- ifelse(round(p.val$p.het,2)>0.99, ">0.99",
                      ifelse(p.val$p.het>=0.01, format(round(p.val$p.het,2),nsmall=2),
                             ifelse((p.val$p.het>=0.001 & p.val$p.het<0.01), format(round(p.val$p.het,3),nsmall=3),
                                    format(p.val$p.het,digits=3))))

# Combine results - first row for each metabolite
results <- merge(final.results[[4]][c("col")], final.results[[8]][c("col")], by="row.names", suffixes = c(".2",".3"))
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
                                                   col.2 = final.results[[5]][i,c("col")],
                                                   col.3 = final.results[[9]][i,c("col")],
                                                   col.4 = final.results[[12]][i,c("col")],
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
                                                     col.2 = final.results[[6]][i,c("col")],
                                                     col.3 = final.results[[10]][i,c("col")],
                                                     col.4 = final.results[[13]][i,c("col")],
                                                     col.5 = "",
                                                     .before = i)
}
results.final.cols.4 <- results.final.cols.3

# Add fourth row for each metabolite
for(i in 1:dim(results.all.cols)[1]){
  results.final.cols.4 <- results.final.cols.4 %>% add_row(class_broad = results.final.cols$class_broad[i],
                                                           hmdb_id = paste(results.final.cols$hmdb_id[i],"_4"),
                                                           metabolite_name = "",
                                                           col.1 = final.results[[3]][i,c("col")],
                                                           col.2 = final.results[[7]][i,c("col")],
                                                           col.3 = final.results[[11]][i,c("col")],
                                                           col.4 = final.results[[14]][i,c("col")],
                                                           col.5 = p.val$p.het[i],
                                                           .before = i)
}

# Arrange by metabolite class and hmdb_id
results.final.cols.4 <- arrange(results.final.cols.4, class_broad, hmdb_id)

# Save table
write.xlsx(results.final.cols.4, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_9.xlsx"), overwrite=TRUE)


# Check number of p-values for heterogeneity that are nominally significant at alpha = 0.05
p_cutoff <- 0.05
full.adjust.results <- as.data.frame(p.het)
length(full.adjust.results[ which(full.adjust.results$p.het<p_cutoff), ])
