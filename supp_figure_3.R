##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_figure_3.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Plot frequencies of known metabolites by class, annotated at the Broad Institute (Supplemental Figure 3)

# Statistical Analyses:
    # Descriptive statistics

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

# Load known eset
load(file.path(main.dir,"data/known.eset.final.RData"))

# Remove metabolites with >=10% missing values
met.vals <- exprs(known.eset.final)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
miss.names <- row.names(met.vals[which(missing >= 10),])
source(file.path(main.dir,"function_removeMetabolitesFromESet.R"))
known.eset.final = removeMetabolitesFromESet(eset = known.eset.final, metabolites.to.remove = miss.names)

# Extract metabolite fdata
fdata <- fData(known.eset.final)

# Recode metabolite classes
table(fdata$class_broad, useNA="always")
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
                                    "Unclassified",
                                    "Other",
                                    "Triglycerides",
                                    "Phosphatidylcholines",
                                    "Phosphatidylcholine plasmalogens",
                                    "Organoheterocyclic compounds",
                                    "Organic acids and derivatives",
                                    "Diglycerides",
                                    "Cholesteryl esters",
                                    "Carnitines",
                                    "Carboxylic acids and derivatives"))
table(fdata$class, useNA = "always")


##### 2 - CREATE FIGURE ####

# Define color scheme
mycolors <- c(brewer.pal(7,"Dark2"),"lightskyblue","steelblue","#FFD700","#696969")

# Plot figure
hist <- ggplot(fdata, aes(y=class, fill=class)) + 
  geom_bar(stat="count") +
  scale_fill_manual(values = rev(mycolors))+
  labs(x = "Count", y = "Class", fill="class_broad") +
  scale_x_continuous(breaks=c(0,25,50,75,100,125,150,175,200)) +
  stat_count(geom = "text", colour = "black", size = 8,
             aes(label = ..count..),position=position_dodge(width=0.9),hjust=-0.25)+
  theme_classic() +
  theme( 
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(color="black", size=18),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(color="black", size=22))
plot(hist)

# Save figure
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_3.png"), device = "png", width = 20, height = 9, units = "in", plot = hist)  
