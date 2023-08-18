##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_figure_4.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Plot CVs, ICCs, and % missing for known and unknown metabolite, stratified by metabolite class (Supplemental Figure 4)

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
library(ggplot2)
library(tidyverse)
library(Biobase)
library(RColorBrewer)
library(grid)

#### KNOWN METABOLITES ####
# Load known eset
load(file.path(main.dir,"data/known.eset.final.RData"))

# Remove metabolites with >=10% missing values
met.vals <- exprs(known.eset.final)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
miss.names <- row.names(met.vals[which(missing >= 10),])
source(file.path(main.dir,"function_removeMetabolitesFromESet.R"))
known.eset.final = removeMetabolitesFromESet(eset = known.eset.final, metabolites.to.remove = miss.names)

# Add CV and ICC data and class
known.fdata <- fData(known.eset.final)
myvars <- c("mean_cv","mean_icc","class_broad")
known.data <- known.fdata[myvars]

# Recode metabolite classes
known.data <- mutate_at(known.data, c("class_broad") , ~replace(., is.na(.), "Missing"))
known.data <- mutate(known.data, class= factor(class_broad),
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
table(known.data$class)
known.data <- select(known.data, -c(class_broad))

# Check CV, ICC
summary(known.data$mean_cv)
summary(known.data$mean_icc)

#### 2 - CREATE FIGURE ####

# Define color scheme
mycolors <- c(brewer.pal(7,"Dark2"),"lightskyblue","steelblue","#FFD700","#696969","Black")

# Plot CVs by class
cvplot <- ggplot(known.data, aes(x = class, y = mean_cv, fill = as.factor(class))) +
  geom_boxplot(na.rm = TRUE)+
  scale_fill_manual(values=mycolors)+
  labs(x="", y = "CV") +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),labels = c("0%", "25%", "50%", "75%", "100%")) +
  theme_minimal() +
  theme( 
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.title.y = element_text(color="black", size=22, face="bold"),
    axis.text.y = element_text(size=18, color="black"))
print(cvplot)
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_4_CV.png"), device = "png", width = 20, height = 5.7, units = "in", plot = cvplot)  


# Plot ICCs
iccplot <- ggplot(known.data, aes(x = class, y = mean_icc, fill = as.factor(class))) +
  geom_boxplot(na.rm = TRUE)+
  scale_fill_manual(values=mycolors)+
  labs(x = "", y = "ICC") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  theme( 
    legend.position = "none",
    axis.line = element_line(colour = "black"),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(color="black", size=22, face="bold"),
    axis.title.y = element_text(color="black", size=22, face="bold"),
    axis.text.y = element_text(size=18, color="black"))
print(iccplot)
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_4_ICC.png"), device = "png", width = 20, height = 5.7, units = "in", plot = iccplot)  



# Plot legend
legendplot <- ggplot(known.data, aes(x = class, y = mean_icc, color = as.factor(class))) +
  geom_point(shape=15)+
  scale_color_manual(values=mycolors, labels=c("Carboxylic acids and derivatives (n=16)            ",
                                              "Carnitines (n=20)",
                                              "Cholesteryl esters (n=12)",
                                              "Diglycerides (n=11)",
                                              "Organic acids and derivatives (n=28)    ",
                                              "Organoheterocyclic compounds (n=8)             ",
                                              "Phosphatidylcholine plasmalogens (n=13)          ",
                                              "Phosphatidylcholines (n=21)", 
                                              "Triglycerides (n=82)",
                                              "Other (n=84)",
                                              "Unclassified (n=39)"))+
  labs(x = "", y = "% missing", color="") +
  scale_y_continuous(limits = c(0, 100)) +
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme_minimal() +
  theme( 
    legend.position = "bottom",
    legend.text = element_text(color = "black", size = 16),
    legend.title = element_text(color = "black", size = 12),
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
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_4_legend.png"), device = "png", width = 18, height = 3, units = "in", plot = legend)  
