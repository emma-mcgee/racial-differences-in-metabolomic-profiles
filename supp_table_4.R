##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_4.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Metabolite set enrichment analysis for NHS, observed and residual (Supplemental Table 4)

# Statistical Analyses:
    # Metabolite set enrichment analysis based an metabolite classes obtained from the Broad lab
    # Regression coefficients are based on linear regression models for observed differences and counterfacutla inverse odds ratio weighting-based models for residual differences

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####

main.dir <- ""
results.dir <- file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(fgsea)
library(reshape)
library(RColorBrewer)
library(openxlsx)
library(stringr)

# Source observed and residual difference estimates - NHS

# WHI observed
source(file = file.path(main.dir,"supp_table_6.R"))
results.whi <- as.data.frame(lin.res.model1)

# NHS observed
source(file = file.path(main.dir,"supp_table_2.R"))
results.nhs <- as.data.frame(lin.res.model0)

# NHS residual
library("readxl")
results_residual <- read_excel(file.path(main.dir,"data/mediation_results_iorw_based_all_mediators.xlsx"))
results.nhs.resid <- as.data.frame(results_residual)
row.names(results.nhs.resid) <- row.names(results.nhs)

# Prepare fdata to be added to difference data
load(file.path(main.dir,"data/known.eset.final.RData"))
fdata <- fData(known.eset.final)
myvars <- c("hmdb_id","metabolite_name","class_broad")
fdata <- fdata[myvars]

# Merge fdata with results
results.nhs.f <- merge(fdata, results.nhs, by="row.names")
results.whi.f <- merge(fdata, results.whi, by="row.names")
results.nhs.resid.f <- merge(fdata, results.nhs.resid, by="row.names")

list.results <- list(results.nhs.f, results.whi.f, results.nhs.resid.f)

# Recode PC plasmalogen 
list.results <- lapply(list.results, function(x) {
  x[,"class_broad"] <-  ifelse(x$class_broad=="PC plasmalogens", "Phosphatidylcholine plasmalogens", x$class_broad)
  x
})

# Extract # of double bonds for each TAG and
# Recode TAGs into two separate groups (<3 double bonds vs. >=3 double bonds)
list.results <- lapply(list.results, function(x) {
   x$db <-  ifelse(x$class_broad=="Triglycerides", 
                                str_extract(x$metabolite_name, "(?<=:).*(?= TAG)"), 
                                NA)
  x$db <- ifelse(x$class_broad=="Triglycerides" & str_sub(x$metabolite_name, -19,-1) == "DAG_or_TAG_fragment", 
                 str_extract(x$metabolite_name, "(?<=:).*(?= DAG_or_TAG_fragment)"), 
                 x$db)
  x$class_broad_new <- ifelse((x$class_broad=="Triglycerides" & (x$db==0 | x$db==1 | x$db==2)), "Triglycerides with <3 double bonds", 
                              ifelse(x$class_broad=="Triglycerides", "Triglycerides with ≥3 double bonds", x$class_broad))
  x
})

results.nhs.f <- list.results[[1]]
results.whi.f <- list.results[[2]]
results.nhs.resid.f <- list.results[[3]]

# Format data for FGSEA package - annotations stored in met.classes
met.classes <- vector("list",length = length(unique(results.nhs.f$class_broad_new))) 
names(met.classes) <- unique(results.nhs.f$class_broad_new)
met.classes[["NA"]] = NULL

for(i in 1:length(met.classes)){
  met.classes[[i]] = results.nhs.f$metabolite_name[which(results.nhs.f$class_broad_new == names(met.classes)[i])]
}

# Create wrapper function for FGSEA analysis
runFGSEA <- function(pathways = met.classes, stats = betas, minSize=3, maxSize=500, nperm=1000, file.result){
  # pathways = met.classes (where annotations are stored)
  # stats = betas (use beta coefficients from linear regression to order the metabolites)
          # prefer to use Betas to order metabolites rather than p-values since p-values don't convey directionality
  # minSize = 3 (only include classes that have at least 3 metabolites mapped to it)
  # maxSize = 500 (only include classes that have less than 500 metabolites mapped to it)
  # nperm = 1000 (number of permutations run to calculate p-values)
  # file.result = file.path(XX.dir,"XXX.csv")
  
  set.seed(08347) # Set seed
  
  fgsea.res <- fgsea(pathways = pathways, stats = stats, minSize=minSize, maxSize=maxSize, nperm=nperm) # Run MSEA
  fgsea.res <- fgsea.res %>%   # Sort by abs(NES) and number of metabolites (in case of NES ties)
    arrange(desc(abs(round(NES,2))), desc(size))
  to.print <- fgsea.res[,-8]  # Drop unnecessary columns
  
  write.table(x = to.print, file = file.result ,col.names = T, row.names = F, sep = "\t", quote = F) # write and print table
  return(to.print)
  
}


#### 2 - RUN MSEA ANALYSES ####

# NHS - extract data
betas <- as.numeric(results.nhs.f$effect.estimate) 
names(betas) <- results.nhs.f$metabolite_name
# Run MSEA
fgsea.nhs <- runFGSEA(pathways = met.classes, stats = betas, file.result = file.path(results.dir,"figures/figure_2_data/msea_results_nhs.csv"))

# WHI - extract data
betas <- as.numeric(results.whi.f$effect.estimate)
names(betas) <- results.whi.f$metabolite_name
# Run MSEA
fgsea.whi <- runFGSEA(pathways = met.classes, stats = betas, file.result = file.path(results.dir,"figures/figure_2_data/msea_results_whi.csv"))

# NHS residual - extract data
betas <- as.numeric(results.nhs.resid.f$effect.estimate.nde)
names(betas) <- results.nhs.resid.f$metabolite_name
# Run MSEA
fgsea.nhs.resid <- runFGSEA(pathways = met.classes, stats = betas, file.result = file.path(results.dir,"figures/figure_2_data/msea_results_nhs_residual.csv"))



#### 2 - CREATE SUPPLEMENTAL TABLE 4 ####

# Format numbers
list.results <- list(fgsea.nhs, fgsea.whi, fgsea.nhs.resid)
list.df <- lapply(list.results, as.data.frame)
list.results.format <- lapply(list.df, function(x) cbind(x[c("pathway")], x[c("size")], format(round(x[c("NES")],2),nsmall=2),
                                                              "pval"=ifelse(round(x$pval,3)>0.000, format(round(x$pval,3),nsmall=3),
                                                                               format(x$pval,digits=3)),
                                                              "padj"=ifelse(round(x$padj,3)>0.000, format(round(x$padj,3),nsmall=3),
                                                                                format(x$padj,digits=3))))

# Merge data
all_msea <- merge(list.results.format[[1]][,c("pathway","size","NES","pval")],
                  list.results.format[[3]][,c("pathway","NES","pval")], by="pathway", suffixes = c(".obs",".resid"))


# Sort by abs(NES) and number of metabolites (in case of NES ties), based on the observed differences in NHS
all_msea <- all_msea %>% 
  arrange(desc(abs(as.numeric(NES.obs))), desc(size))

# Format p-values
all_msea$pval.resid <- as.numeric(all_msea$pval.resid)
all_msea$pval.resid <- ifelse(round(all_msea$pval.resid,2)>0.99, ">0.99",
                            ifelse(all_msea$pval.resid>=0.01, format(round(all_msea$pval.resid,2),nsmall=2),
                                   ifelse((all_msea$pval.resid>=0.001 & all_msea$pval.resid<0.01), format(round(all_msea$pval.resid,3),nsmall=3),
                                                                                            format(all_msea$pval.resid,digits=3))))
all_msea$pval.obs <- as.numeric(all_msea$pval.obs)
all_msea$pval.obs <- ifelse(round(all_msea$pval.obs,2)>0.99, ">0.99",
                              ifelse(all_msea$pval.obs>=0.01, format(round(all_msea$pval.obs,2),nsmall=2),
                                     ifelse((all_msea$pval.obs>=0.001 & all_msea$pval.obs<0.01), format(round(all_msea$pval.obs,3),nsmall=3),
                                            format(all_msea$pval.obs,digits=3))))

# Save table
write.xlsx(all_msea, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_4.xlsx"),
           rowNames=FALSE, overwrite = TRUE)
