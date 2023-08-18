##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_table_7.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Metabolite set enrichment analysis for WHI replication (Supplemental Table 7)

# Statistical Analyses:
    # Metabolite set enrichment analysis based an metabolite classes obtained from the Broad lab
    # Regression coefficients are based on linear regression models for observed differences and counterfactual inverse odds ratio weighting-based analyses for residual differences

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####

# Source MSEA
source(file = file.path(main.dir,"supp_table_4.R"))

# Prep data
all_msea <- list.results.format[[2]][,c("pathway","size","NES","pval")]

# Sort by abs(NES) (and number of metabolites, in case of NES ties)
all_msea <- all_msea %>% 
  arrange(desc(abs(as.numeric(NES))), desc(size))

# Format p-values
all_msea$pval <- as.numeric(all_msea$pval)
all_msea$pval <- ifelse(round(all_msea$pval,2)>0.99, ">0.99",
                            ifelse(all_msea$pval>=0.01, format(round(all_msea$pval,2),nsmall=2),
                                   ifelse((all_msea$pval>=0.001 & all_msea$pval<0.01), format(round(all_msea$pval,3),nsmall=3),
                                          format(all_msea$pval,digits=3))))

# Save table
write.xlsx(all_msea, file = file.path(main.dir,"final_results/supplemental_tables/supp_table_7.xlsx"),
           rowNames=FALSE, overwrite = TRUE)
