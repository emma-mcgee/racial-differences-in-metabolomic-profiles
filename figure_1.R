#########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: figure_1.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Plot of observed and residual differences in individual metabolites (Figure 1)

# Statistical Analyses:
    # None, visualization of data

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####

rm(list = ls())
main.dir = ""
results.dir = file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(tidyverse)
library(Biobase)
library(ggplot2)
library(RColorBrewer)
library(grid)

# Source observed difference estimates
source(file = file.path(main.dir,"supp_table_2.R"))

# Source results for residual- bootstrapping takes a long time to run, so results are sourced
library("readxl")
results_residual <- read_excel(file.path(main.dir,"data/mediation_results_iorw_based_all_mediators.xlsx"))
results_residual <- as.data.frame(results_residual)
row.names(results_residual) <- row.names(results)
results_residual$effect.estimate.nde <- as.numeric(results_residual$effect.estimate.nde)

# Load known eset
load(file.path(main.dir,"data/known.eset.final.RData"))

# Remove metabolites with >=10% missing values
met.vals <- exprs(known.eset.final)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
miss.names <- row.names(met.vals[which(missing >= 10),])
source(file.path(main.dir,"function_removeMetabolitesFromESet.R"))
known.eset.final = removeMetabolitesFromESet(eset = known.eset.final, metabolites.to.remove = miss.names)

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
table(fdata$class)
# Subset fdata to relevant variables
myvars <- c("hmdb_id","metabolite_name","class")
fdata  <- fdata[myvars]

# Merge fdata, NHS observed, and NHS residual results

# Observed data
results.nhs <- as.data.frame(lin.res.model0)
results.nhs$study <- "Observed"
results.observed <- merge(results.nhs, fdata, by="row.names")
results.observed <- select(results.observed, -c("standard.error"))

# Residual data
results.resid <- as.data.frame(results_residual)
results.resid$study <- "Residual"
results.resid.f <- merge(results.resid, fdata, by="row.names")

# Format residual differences for merging
results.resid.f$L_CI <- as.numeric(substr(results.resid.f$conf.int.nde,2,6))
results.resid.f$U_CI <- as.numeric(substr(results.resid.f$conf.int.nde,8,13))
results.resid.f  <- results.resid.f[,c("Row.names","effect.estimate.nde","p.value.nde",
                                       "L_CI","U_CI","study","hmdb_id","metabolite_name","class")]
colnames(results.resid.f) <- colnames(results.observed)

# Combine
results.all <- rbind(results.observed, results.resid.f)

# Re-level obs/residual indicator
results.all <- mutate(results.all, study=factor(study),
                      study = fct_relevel(study,
                                          "Observed",
                                          "Residual"))

# Order results by obs/residual, metabolite class, and HMDB ID
results.all <- results.all %>% 
  arrange(study, class, hmdb_id) 

# Create sequence numbers for plotting
results.all$seq <- rep(seq(1, dim(results.nhs)[1], by=1), 2)

#### 2 - CREATE FIGURE ####

#### Panel a) Observed ####
results.obs <- results.all[which(results.all$study=="Observed"),]
dim(results.obs)
table(results.obs$study)

# Define color scheme
mycolors <- c(brewer.pal(7,"Dark2"),"lightskyblue","steelblue","#FFD700","#696969")

# Create plot
barplot<-ggplot(data=results.obs,
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
    legend.position = "none",
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20, face="bold"),
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgrey", size=1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(color="black", size=18),
    axis.title.y = element_text(color="black", size=22, face="bold"),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.text.x = element_blank())
barplot
ggsave(filename=file.path(results.dir,"figures/figure_1_a.png"), device = "png", width = 18, height = 4.5, units = "in", plot = barplot)  


#### Panel b) Residual ####
results.res <- results.all[which(results.all$study=="Residual"),]
dim(results.res)
table(results.res$study)

# Create plot
barplot<-ggplot(data=results.res,
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
    legend.position = "none",
    legend.text = element_text(color = "black", size = 20),
    legend.title = element_text(color = "black", size = 20, face="bold"),
    axis.line = element_line(colour = "black"),
    strip.text.x = element_text(size = 26),
    panel.border = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgrey", size=1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(color="black", size=18),
    axis.title.y = element_text(color="black", size=22, face="bold"),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.text.x = element_blank())
barplot
ggsave(filename=file.path(results.dir,"figures/figure_1_b.png"), device = "png", width = 18, height = 4.5, units = "in", plot = barplot)  


# Plot legend
legendplot<-ggplot(data=results.obs,
                aes(x=seq, y=effect.estimate, ymin=L_CI, ymax=U_CI,
                    fill=class)) +
  geom_bar(stat="identity") +
  geom_errorbar(width=0.3, size=0.3, color="black") +
  scale_y_continuous(limits=c(-1.13,1.13), breaks = c(-1,-0.5,0,0.5,1),
                     name="Mean difference") +
  geom_hline(yintercept=0, color="black", linetype="solid") +
  scale_fill_manual(values=mycolors, labels=c("Carboxylic acids and derivatives                 ",
                                               "Carnitines",
                                               "Cholesteryl esters",
                                               "Diglycerides",
                                               "Organic acids and derivatives    ",
                                               "Organoheterocyclic compounds             ",
                                               "Phosphatidylcholine plasmalogens          ",
                                               "Phosphatidylcholines", 
                                               "Triglycerides",
                                               "Other",
                                               "Unclassified"))+
  theme_minimal() +
  theme( 
    legend.position = "bottom",
    legend.text = element_text(color = "black", size = 16),
    legend.title = element_blank(),
    legend.direction="horizontal",
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold"),
    axis.text.y = element_text(size=14, color="black"))
legend <- cowplot::get_legend(legendplot)
grid.newpage()
grid.draw(legend)
ggsave(filename=file.path(results.dir,"figures/figure_1_legend.png"), device = "png", width = 18, height = 3, units = "in", plot = legend)  

