##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_6.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: quantify observed differences in individual known metabolites between Black and White women in WHI (Supplemental Table 6)

# Statistical Analyses:
    # Linear regression with each ln Z-score transformed metabolite as the outcome variable
    # 95% prediction intervals to evaluate replication

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

#### WHI ####
source(file = file.path(main.dir,"whi_data.R"))
whi.met.vals.orig <- whi.met.vals
whi.met.vals <- t(whi.met.vals.orig)

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

# Prepare regression output files - WHI
lin.res.model1 <- matrix(0,ncol = 5,nrow = dim(whi.met.vals)[1])
colnames(lin.res.model1) <- c("effect.estimate","standard.error","p.value","L_CI","U_CI")
rownames(lin.res.model1) <- colnames(whi.met.vals.orig)


#### 2 - QUANTIFY OBSERVED DIFFERENCES USING LINEAR REGRESSION ####

#### NHS ####
# Run linear regressions to estimate observed differences and prediction intervals
for(i in 1:dim(met.vals)[1]){
  curr.data <- cbind(pdata[,c("race_bw",
                             "ipw_c",
                             "id",
                             "knn_nhs1_agebld","knn_nhs1_agebld_sqr","knn_nhs1_menopmh_bld","knn_nhs1_fast","knn_nhs1_bldmonth",
                                              "knn_nhs1_btimelr")], met = met.vals[i,])
  
  # NHS: adjusted for age
  my.lm.model0 <- lm(met~race_bw
                     + knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                      data = curr.data)
  my.model.model0 <- summary(my.lm.model0)      
  
  # Prediction interval
  # Using the formula for a prediction interval for a mean
  # w/ replication sample of different size than the original study
  # Based on methods described by:
  # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0162874#:~:text=The%20prediction%20interval%20provides%20an,due%20to%20sample%20error%20alone. 
  
  n1 <- dim(pdata)[1]
  n2 <- dim(whi.met.vals)[2]
  t <- qt(p=0.05, df=summary(my.lm.model0)$df[2], lower.tail=FALSE)
  var <- sqrt((my.model.model0$coefficients["race_bw1","Std. Error"]*sqrt(n1))^2/n1 + (my.model.model0$coefficients["race_bw1","Std. Error"]*sqrt(n1))^2/n2)
  U_CI <- my.model.model0$coefficients["race_bw1","Estimate"] + (t*var)
  L_CI <- my.model.model0$coefficients["race_bw1","Estimate"] - (t*var)
  
  lin.res.model0[i,] <- c(my.model.model0$coefficients["race_bw1",c("Estimate","Std. Error","Pr(>|t|)")],L_CI,U_CI)
}


#### WHI - Replication for overlapping metabolites ####
# Run linear regressions to estimate observed differences in metabolite levels
for(i in 1:dim(whi.met.vals)[1]){
  curr.data <- cbind(whi.pdata[,c("race_bw",
                              "age","age_sqr")],met = whi.met.vals[i,])
  
  # WHI: adjusted for age
  my.lm.model1 <- lm(met ~ race_bw + age + age_sqr,
                     data = curr.data)
  my.model.model1 <- summary(my.lm.model1)      
  my.ci.model1 <- confint(my.lm.model1)

  
  lin.res.model1[i,] <- c(my.model.model1$coefficients["race_bw1",c("Estimate","Std. Error","Pr(>|t|)")],my.ci.model1["race_bw1",])
}

### CHECK REPLICATION ###

# Check how many WHI values fall within NHS prediction interval
rep.check <- merge(lin.res.model0, lin.res.model1, by="row.names", all.x = TRUE, suffixes = c(".nhs",".whi"))
results.final.rep <- rep.check[which((rep.check$effect.estimate.whi > rep.check$L_CI.nhs) & 
                                       (rep.check$effect.estimate.whi < rep.check$U_CI.nhs)),]
# Number within prediction interval
dim(results.final.rep)[1]
# Proportion within prediction interval (of all replicated metabs in WHI)
dim(results.final.rep)[1]/dim(lin.res.model1)[1]
# Note that some are NOT within prediction interval b/c the association was stronger in the replication

# Number within prediction interval or more extreme
results.final.rep <- rep.check[which((rep.check$effect.estimate.whi > rep.check$L_CI.nhs) & 
                                       (rep.check$effect.estimate.whi < rep.check$U_CI.nhs) | #falls within prediction interval
                                  (rep.check$effect.estimate.whi > rep.check$U_CI.nhs) & # or is >0 in both and more extreme in WHI
                                      (sign(rep.check$effect.estimate.whi)==1) &
                                      (sign(rep.check$effect.estimate.whi)==sign(rep.check$effect.estimate.nhs)) |
                                 (rep.check$effect.estimate.whi < rep.check$L_CI.nhs) & # or is <0 in both and more extreme in WHI
                                 (sign(rep.check$effect.estimate.whi)==-1) &
                                 (sign(rep.check$effect.estimate.whi)==sign(rep.check$effect.estimate.nhs))),]
# Number within prediction interval OR more extreme
dim(results.final.rep)[1]
# Proportion within prediction interval OR more extreme (of all replicated metabs in WHI)
dim(results.final.rep)[1]/dim(lin.res.model1)[1]

# Check correlation between WHI and NHS estimates
all.res <- merge(lin.res.model1[,c("effect.estimate")], lin.res.model0[,c("effect.estimate")], by="row.names")
cor(all.res$x, all.res$y, method = "spearman")


# Check how many of those with absolute differences >=0.5 in NHS are replicated
results.final.rep <- rep.check[which((rep.check$effect.estimate.whi > rep.check$L_CI.nhs) & 
                                       (rep.check$effect.estimate.whi < rep.check$U_CI.nhs) &
                                       (abs(rep.check$effect.estimate.nhs) >= 0.5)),]
# Number within prediction interval
dim(results.final.rep)[1]
# Proportion within prediction interval (of all replicated metabs in WHI)
dim(results.final.rep)[1]/dim(rep.check[which(abs(rep.check$effect.estimate.nhs) >= 0.5 & is.na(rep.check$effect.estimate.whi) ==FALSE),])[1]

# Number within prediction interval or more extreme (of those with effect estimate >=0.5 SD)
results.final.rep <- rep.check[which((rep.check$effect.estimate.whi > rep.check$L_CI.nhs) &
                                       (abs(rep.check$effect.estimate.nhs) >= 0.5) & 
                                       (rep.check$effect.estimate.whi < rep.check$U_CI.nhs) | #falls within prediction interval
                                       (rep.check$effect.estimate.whi > rep.check$U_CI.nhs) & # or is >0 in both and more extreme in WHI
                                       (abs(rep.check$effect.estimate.nhs) >= 0.5) & 
                                       (sign(rep.check$effect.estimate.whi)==1) &
                                       (sign(rep.check$effect.estimate.whi)==sign(rep.check$effect.estimate.nhs)) |
                                       (rep.check$effect.estimate.whi < rep.check$L_CI.nhs) & # or is <0 in both and more extreme in WHI
                                       (abs(rep.check$effect.estimate.nhs) >= 0.5) & 
                                       (sign(rep.check$effect.estimate.whi)==-1) &
                                       (sign(rep.check$effect.estimate.whi)==sign(rep.check$effect.estimate.nhs))),]
dim(results.final.rep)[1]
# Proportion within prediction interval (of all replicated metabs in WHI)
dim(results.final.rep)[1]/dim(rep.check[which(abs(rep.check$effect.estimate.nhs) >= 0.5 & is.na(rep.check$effect.estimate.whi) ==FALSE),])[1]
dim(rep.check[which(abs(rep.check$effect.estimate.nhs) >= 0.5 & is.na(rep.check$effect.estimate.whi) ==FALSE),])[1]

#### 3 - CREATE TABLE ####

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

# Combine results from NHS and WHI, restricted to metabolites included in WHI
results <- merge(final.results[[1]], final.results[[2]], by="row.names", all.y = TRUE, suffixes = c(".nhs",".whi"))

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
write.xlsx(results.final, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_6.xlsx"),
           overwrite=TRUE)


