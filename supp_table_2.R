##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_2.R

# Programmer: Emma McGee 

# Date Last Updated: February 8, 2023

# Purpose of Program: quantify observed differences in individual known metabolites between Black and White women (Supplemental Table 2)

# Statistical Analyses:
    # Linear regression with each ln Z-score transformed metabolite as the outcome variable
    # Inverse probability weighting used to account for the original case-control design
    # Valid confidence intervals and p-values are obtained from the robust sandwich variance
    # Multiple testing correction using the number of effective tests with 99.5% of the variance explained

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####

main.dir <- ""
results.dir <- file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(tidyverse)
library(Biobase)
library(geepack)
library(openxlsx)

#### NHS ####
# Load eset of known metabolites
load(file.path(main.dir, "data/known.eset.final.RData"))

# Prepare data
pdata <- pData(known.eset.final)
met.vals <- exprs(known.eset.final)

# Restrict to metabolites with 0% missing data after imputation (<10% missing in raw data)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
met.vals <- met.vals[which(missing == 0),]
dim(met.vals)

# Prepare regression output files - NHS
lin.res.model0 <- matrix(0,ncol = 5,nrow = dim(met.vals)[1])
colnames(lin.res.model0) <- c("effect.estimate","standard.error","p.value","L_CI","U_CI")
rownames(lin.res.model0) <- rownames(met.vals)

# Prepare regression output files - NHS IPW
lin.res.model1 <- matrix(0,ncol = 5,nrow = dim(met.vals)[1])
colnames(lin.res.model1) <- c("effect.estimate","standard.error","p.value","L_CI","U_CI")
rownames(lin.res.model1) <- rownames(met.vals)


#### 2 - QUANTIFY OBSERVED DIFFERENCES USING LINEAR REGRESSION ####

#### NHS ####
# Run linear regressions to estimate observed differences in metabolite levels
for(i in 1:dim(met.vals)[1]){
  curr.data <- cbind(pdata[,c("race_bw",
                             "ipw_c",
                             "id",
                             "knn_nhs1_agebld","knn_nhs1_agebld_sqr","knn_nhs1_menopmh_bld","knn_nhs1_fast","knn_nhs1_bldmonth",
                                              "knn_nhs1_btimelr")], met = met.vals[i,])
  
  # NHS: adjusted for age
  my.lm.model0 <- lm(met~race_bw
                     + knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                     # + knn_nhs1_menopmh_bld + knn_nhs1_fast + knn_nhs1_bldmonth + knn_nhs1_btimelr,
                      data = curr.data)
  my.model.model0 <- summary(my.lm.model0)      
  my.ci.model0 <- confint(my.lm.model0)
  
  # NHS: adjusted for age + IPW to account for sampling
  my.lm.model1 <- geeglm(met ~ race_bw + knn_nhs1_agebld + knn_nhs1_agebld_sqr, data=curr.data, weights=ipw_c, id=id, corstr="independence")
  my.model.model1 <- summary(my.lm.model1) 
  beta <- coef(my.lm.model1)
  SE <- coef(summary(my.lm.model1))[,"Std.err"]
  lcl <- beta-qnorm(0.975)*SE 
  ucl <- beta+qnorm(0.975)*SE
  my.ci.model1 <- cbind(lcl, ucl)
  
  lin.res.model0[i,] <- c(my.model.model0$coefficients["race_bw1",c("Estimate","Std. Error","Pr(>|t|)")],my.ci.model0["race_bw1",])
  lin.res.model1[i,] <- c(data.matrix(my.model.model1$coefficients["race_bw1",c("Estimate","Std.err","Pr(>|W|)")]),my.ci.model1["race_bw1",])
}

#### 3 - CORRECT FOR MULTIPLE TESTING ####

# Calculate the number of effective tests with 99.5% of the variance explained, as described by Gao et al.
# Divide nominal p-value of 0.05 by the number of effective tests
# Use the resulting p-value as the threshold for multiple testing correction

# Restrict to metabolites with 0% missing data after imputation (<10% missing in raw data)
vals <- exprs(known.eset.final)
missing <- apply(vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
vals.pca <- vals[which(missing == 0),]
dim(vals.pca)
 
# Perform PCA
pca.res <- prcomp(t(vals.pca))
 
# Calculate cumulative sum of % variance explained
cum.exp.var <- cumsum(pca.res$sdev^2/sum(pca.res$sdev^2))
 
# Return number of PCAs that explain 99.5% of the variance
nef <- which(cum.exp.var>=0.995)[1]
nef

# Define threshold based on number of effective tests
p_cutoff <- 0.05/nef
p_cutoff

# Check # of metabolites that pass NEF
length(as.data.frame(lin.res.model0)$p.value[which(as.data.frame(lin.res.model0)$p.value<p_cutoff)])
length(as.data.frame(lin.res.model1)$p.value[which(as.data.frame(lin.res.model1)$p.value<p_cutoff)])

# Check correlation
all.res <- merge(lin.res.model1[,c("effect.estimate")], lin.res.model0[,c("effect.estimate")], by="row.names")
cor(all.res$x, all.res$y, method = "spearman")


#### 4 - CREATE TABLE ####

# Prepare data for table

fdata <- fData(known.eset.final)
myvars <- c("class_broad","hmdb_id","metabolite_name")
metab.fdata <- fdata[myvars]
# Recode PC plasmalogen 
metab.fdata$class_broad <- ifelse(metab.fdata$class_broad=="PC plasmalogens", "Phosphatidylcholine plasmalogens", metab.fdata$class_broad)

# Format results for printing
list.results <- list(lin.res.model0, lin.res.model1)
list.results <- lapply(list.results, as.data.frame)
list.results.format <- lapply(list.results, function(x) cbind(format(round(x[c("effect.estimate","L_CI","U_CI")],2),nsmall=2),
                                                              "p.value"=ifelse(round(x$p.value,2)>0.99, ">0.99",
                                                                               ifelse(x$p.value>=0.01, format(round(x$p.value,2),nsmall=2),
                                                                                      ifelse((x$p.value>=0.001 & x$p.value<0.01), format(round(x$p.value,3),nsmall=3),
                                                                                             format(x$p.value,digits=3))))))
                            
list.results.ci <- lapply(list.results.format, function(x) cbind(x, "conf.int"=paste("(", x$L_CI, ", ", x$U_CI, ")",sep = "")))
final.results <- lapply(list.results.ci, function(x) x[,c("effect.estimate","conf.int","p.value")])


# Combine results
results <- merge(final.results[[1]], final.results[[2]], by="row.names", all.x = TRUE, suffixes = c(".main",".ipw"))

# Merge metabolite names with results
row.names(results) <- results$Row.names
results.final <- merge(metab.fdata, results, by="row.names")
results.final <- results.final[ , !(names(results.final) %in% c("Row.names","Row.names.y"))]

# Arrange by metabolite class and hmdb_id
results.final <- arrange(results.final, class_broad, hmdb_id)

# Format missing HMDB ID and class information
results.final$hmdb_id <- ifelse(is.na(results.final$hmdb_id) == TRUE, "Unavailable", results.final$hmdb_id)
results.final$class_broad <- ifelse(is.na(results.final$class_broad) == TRUE, "Unclassified", results.final$class_broad)

# Save table
write.xlsx(results.final, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_2.xlsx"),
           overwrite=TRUE)
