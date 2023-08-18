##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_figure_6.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: quantify observed differneces in individual known metabolites between Black and White women, restricted to fasting samples (Supplemental Figure 6)

# Statistical Analyses:
    # Linear regression with each ln Z-score transformed metabolite as the outcome variable

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
library(ggplot2)
library(RColorBrewer)

# Load eset of known metabolites
load(file.path(main.dir, "data/known.eset.final.RData"))

# Restrict to women who are fasting at blood collection (n=1442)
table(known.eset.final$knn_nhs1_fast, useNA="always")
table(known.eset.final$knn_nhs1_fast, known.eset.final$race_bw)

dim(known.eset.final)
known.eset.final.fast <- known.eset.final[ , (known.eset.final$knn_nhs1_fast == "1")]
dim(known.eset.final.fast)

# Prepare data
pdata <- pData(known.eset.final.fast)
met.vals <- exprs(known.eset.final.fast)

# Restrict to metabolites with 0% missing data after imputation (<10% missing in raw data)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
met.vals <- met.vals[which(missing == 0),]
dim(met.vals)

# Prepare regression output files
fast.model0 <- matrix(0,ncol = 5,nrow = dim(met.vals)[1])
colnames(fast.model0) <- c("effect.estimate","standard.error","p.value","L_CI","U_CI")
rownames(fast.model0) <- rownames(met.vals)


#### 2 - QUANTIFY OBSERVED DIFFERENCES USING LINEAR REGRESSION ####

# Run linear regressions to estimate observed differences in metabolite levels (accounting for matching factors using IPW + adjusted for age)
for(i in 1:dim(met.vals)[1]){
  curr.data <- cbind(pdata[,c("race_bw",
                              "ipw_c",
                              "id",
                              "knn_nhs1_agebld","knn_nhs1_agebld_sqr")],met = met.vals[i,])
  
  # NHS: adjusted for age
  my.lm.model0 <- lm(met~race_bw
                     + knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                     data = curr.data)
  my.model.model0 <- summary(my.lm.model0)      
  my.ci.model0 <- confint(my.lm.model0)
  
  fast.model0[i,] <- c(my.model.model0$coefficients["race_bw1",c("Estimate","Std. Error","Pr(>|t|)")],my.ci.model0["race_bw1",])
}


# Calculate Pearson correlation coefficients for effect estimates restricted to fasting vs. in main analysis
source(file = file.path(main.dir,"supp_table_2.R"))
all.res <- merge(fast.model0[,c("effect.estimate")], lin.res.model0[,c("effect.estimate")], by="row.names")
cor(all.res$x, all.res$y, method = "spearman")


#### 3 - PLOT RESULTS ####

# Extract data for plotting
fdata <- fData(known.eset.final)

# Recode metabolite classes
fdata <- mutate_at(fdata, c("class_broad") , ~replace(., is.na(.), "Missing"))
fdata <- mutate(fdata, class= factor(class_broad),
                class = fct_recode(class,
                                   "Carboxylic acids and derivatives" = "Carboxylic acids and derivatives",
                                   "Carnitines" = "Carnitines",
                                   "Cholesteryl esters" = "Cholesteryl esters",
                                   "Diglycerides" = "Diglycerides",
                                   "Organic acids and derivatives" = "Organic acids and derivatives",
                                   "Organoheterocyclic compounds" = "Organoheterocyclic compounds",
                                   "Phosphatidylcholine plasmalogens" = "Phosphatidylcholine plasmalogens",
                                   "Phosphatidylcholine plasmalogens" = "PC plasmalogens",
                                   "Phosphatidylcholines" = "Phosphatidylcholines",
                                   "Triglycerides" = "Triglycerides",
                                   "Other" = "Alkaloids and derivatives",
                                   "Other" = "Azoles",
                                   "Other" = "Benzene and substituted derivatives",
                                   "Other" = "Benzenoids",
                                   "Other" = "Ceramides",
                                   "Other" = "Diazines",
                                   "Other" = "Fatty Acyls",
                                   "Other" = "Imidazopyrimidines",
                                   "Other" = "Lysophosphatidylcholines",
                                   "Other" = "Lysophosphatidylethanolamines",
                                   "Other" = "Nucleosides, nucleotides, and analogues",
                                   "Other" = "Organic oxygen compounds",
                                   "Other" = "Organonitrogen compounds",
                                   "Other" = "Phenols",
                                   "Other" = "Phosphatidylethanolamine plasmalogens",
                                   "Other" = "Phosphatidylethanolamines",
                                   "Other" = "Phosphatidylinositols",
                                   "Other" = "Phosphatidylserine plasmalogens",
                                   "Other" = "Phosphatidylserines",
                                   "Other" = "Pyridines and derivatives",
                                   "Other" = "Sphingomyelins",
                                   "Other" = "Steroids and steroid derivatives",
                                   "Unclassified" = "Missing"),
                class = fct_relevel(class,
                                    "Carboxylic acids and derivatives",
                                    "Carnitines",
                                    "Cholesteryl esters",
                                    "Diglycerides",
                                    "Organic acids and derivatives",
                                    "Organoheterocyclic compounds",
                                    "Phosphatidylcholine plasmalogens",
                                    "Phosphatidylcholines", 
                                    "Triglycerides",
                                    "Other",
                                    "Unclassified"))

# Subset fdata to relevant variables
myvars <- c("hmdb_id","metabolite_name","class")
fdata  <- fdata[myvars]

# Merge fdata with results
results.all <- merge(fast.model0, fdata, by="row.names")

# Order results by metabolite class, and HMDB ID
results.all <- results.all %>% 
  arrange(class, hmdb_id) 

# Create sequence numbers for plotting
results.all$seq <- seq(1, dim(results.all)[1], by=1)


#### 2 - CREATE FIGURE ####

# Define color scheme
mycolors <- c(brewer.pal(7,"Dark2"),"lightskyblue","steelblue","#FFD700","#696969")

# Create plot
barplot<-ggplot(data=results.all,
                aes(x=seq, y=effect.estimate, ymin=L_CI, ymax=U_CI,
                    fill=class)) +
  geom_bar(stat="identity") +
  geom_errorbar(width=0.3, size=0.3, color="black") +
  scale_y_continuous(limits=c(-1.13,1.13), breaks = c(-1,-0.5,0,0.5,1),
                     name="Mean difference") +
  geom_hline(yintercept=0, color="black", linetype="solid") +
  labs(x = "Metabolites", y = "Mean difference", fill="") +
  scale_fill_manual(values = mycolors)+
  theme_minimal() +
  theme( 
    legend.position = "bottom",
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20, face="bold"),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 26),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgrey", size=1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(color="black", size=20),
    axis.title.y = element_text(color="black", size=25, face="bold"),
    axis.title.x = element_text(color="black", size=25, face="bold"),
    axis.text.x = element_blank())
barplot
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_6.png"), device = "png", width = 20, height = 10, units = "in", plot = barplot)  
