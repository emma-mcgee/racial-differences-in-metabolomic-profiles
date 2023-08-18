##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: table_2.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Create table of observed and residual differences in individual metabolites (Table 2)
                   # Display results for metabolites with absolute differences >=|0.5| z-scores

# Statistical Analyses:
    # None, subset of data analyses run in other files only

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####

rm(list = ls())
main.dir = ""
results.dir = file.path(main.dir,"final_results")
setwd(main.dir)

# Source observed difference estimates
source(file = file.path(main.dir,"supp_table_2.R"))

# Load known eset
load(file.path(main.dir,"data/known.eset.final.RData"))

# Remove metabolites with >=10% missing values
met.vals <- exprs(known.eset.final)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
miss.names <- row.names(met.vals[which(missing >= 10),])
source(file.path(main.dir,"function_removeMetabolitesFromESet.R"))
known.eset.final = removeMetabolitesFromESet(eset = known.eset.final, metabolites.to.remove = miss.names)

# Extract fdata
fdata <- fData(known.eset.final)

# Source residual differences -  bootstrapping takes a long time to run, so results are sourced
library("readxl")
results_residual <- read_excel(file.path(main.dir,"data/mediation_results_iorw_based_all_mediators.xlsx"))
results_residual <- as.data.frame(results_residual)


# Format results for printing
list.results <- list(lin.res.model0)
list.results <- lapply(list.results, as.data.frame)
list.results.format <- lapply(list.results, function(x) cbind(format(round(x[c("effect.estimate","L_CI","U_CI")],2),nsmall=2),
                                                              "p.value"=ifelse(round(x$p.value,3)>0.000, format(round(x$p.value,3),nsmall=3),
                                                                               format(x$p.value,digits=3))))
list.results.ci <- lapply(list.results.format, function(x) cbind(x, "conf.int"=paste("(", x$L_CI, ", ", x$U_CI, ")",sep = "")))
final.results <- lapply(list.results.ci, function(x) x[,c("effect.estimate","conf.int")])

# Format residual results
results.residual <- results_residual[,c("effect.estimate.nde","conf.int.nde")]
row.names(results.residual) <- row.names(fdata)

# Combine results
results <- merge(final.results[[1]], results.residual, by="row.names", all.x = TRUE)

# Merge metabolite names with results
row.names(results) <- results$Row.names
results.final <- merge(metab.fdata, results, by="row.names")
results.final <- results.final[ , !(names(results.final) %in% c("Row.names","Row.names.y"))]

# Arrange by metabolite class and hmdb_id
results.final <- arrange(results.final, class_broad, hmdb_id)

# Format missing HMDB ID and class information
results.final$hmdb_id <- ifelse(is.na(results.final$hmdb_id) == TRUE, "Unavailable", results.final$hmdb_id)
results.final$class_broad <- ifelse(is.na(results.final$class_broad) == TRUE, "Unclassified", results.final$class_broad)

# Subset to data for table 2 - subset of metabolites with absolute differences >= |0.5| z-scores (observed)
# row.names(results.final) <- row.names(fdata)
results.final$abs.effect.estimate <- abs(as.numeric(results.final$effect.estimate))
results.observed.0.5 <- results.final %>%
  arrange(class_broad, hmdb_id) %>%
  filter(abs.effect.estimate >= 0.5) %>%
  select(-abs.effect.estimate) 

# Save table
write.xlsx(results.observed.0.5, file = file.path(main.dir,"final_results/tables/table_2.xlsx"), overwrite=TRUE)

