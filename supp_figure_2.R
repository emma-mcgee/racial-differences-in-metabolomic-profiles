##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_figure_2.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Flow diagram for selection of metabolites (Supplemental Figure 2)

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
library(openxlsx)

# Set up empty grid plot for flow diagram
data <- tibble(x= 1:120, y= 1:120)

data %>% 
  ggplot(aes(x, y)) +
  scale_x_continuous(minor_breaks = seq(10, 120, 10)) +
  scale_y_continuous(minor_breaks = seq(10, 120, 10)) +
  theme_linedraw() ->
  grid

# Load eset data for calculations
load(file.path(main.dir, "data/merged.eset.raw.RData"))
load(file.path(main.dir, "data/known.eset.final.RData"))
load(file.path(main.dir, "data/unknown.eset.final.RData"))

# Set N for main box #1 - all metabolites/peaks
main_1 <- dim(merged.eset.raw)[1]
main_1

# Set N for exclusion box #1 - unknown metabolites
exclude_1 <- dim(unknown.eset.final)[1]
exclude_1

# Set N for main box #2 - all known metabolites
main_2 <- main_1 - exclude_1
main_2

# Set N for exclusion box #2 - failed processing delay pilot
fdata <- fData(merged.eset.raw)
exclude_2 <- dim(fdata[ which(fdata$pm_pilot_fail_heparin=='1'), ])[1]
exclude_2

# Set N for final included metabolites - all known metabolites included in the analysis
group_1 <- main_2 - exclude_2
group_1

# Set N for included known metabolites <10% missing
met.vals <- exprs(known.eset.final)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
group_1_lt10 <- dim(met.vals[which(missing < 10),])[1]
group_1_lt10

# Set N for excluded known metabolites >=10% missing
group_1_ge10 <- dim(met.vals[which(missing >= 10),])[1]
group_1_ge10


#### 2 - CREATE FIGURE ####

# Create flow diagram
flow <- grid +
  # Main box 1: All metabolites
  geom_rect(xmin = 13, xmax=57, ymin=90, ymax=100, color='black',
            fill='white', size=0.25) +
  annotate('text', x=35, y=95,label= paste("All metabolites/features\n",
                                            "(n=", format(main_1, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
  # Exclude box 1: unknown metabolites/features
  geom_rect(xmin =60, xmax=118, ymin=80, ymax=90, color='black',
            fill='white', size=0.25) +
  annotate('text', x=89, y=85,label= paste('Excluded: chemically unidentified\n',
                                             "metabolomic features (n=", format(exclude_1, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
  # Main box 2: All known metabolites
  geom_rect(xmin = 13, xmax=57, ymin=70, ymax=80, color='black',
            fill='white', size=0.25) +
  annotate('text', x=35, y=75,label= paste('Chemically identified metabolites\n',
                                            "(n=", format(main_2, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
  # Exclude box 2: Known metabolites that fail processing delay
  geom_rect(xmin =60, xmax=118, ymin=60, ymax=70, color='black',
            fill='white', size=0.25) +
  annotate('text', x=89, y=65,label= paste("Excluded: metabolites which failed the\n",
                                            "blood processing delay pilot study* (n=", format(exclude_2, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
  # Group box 1: Final known metabolites
  geom_rect(xmin = 13, xmax=57, ymin=50, ymax=60, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 35, y=55,label= paste('Chemically identified metabolites\n',
                                            "which passed the pilot study* (n=", format(group_1, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
  # Final: <10% missing known metabolites
  geom_rect(xmin = 13, xmax=57, ymin=30, ymax=40, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 35, y=35,label= paste('Included metabolites\n',
                                            "(n=", format(group_1_lt10, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
  # Exclude box 3 >=10% missing known metabolites
  geom_rect(xmin =60, xmax=118, ymin=40, ymax=50, color='black',
            fill='white', size=0.25) +
  annotate('text', x= 89, y=45,label= paste('Excluded: metabolites with \u226510% missing values\n',
                                            "(n=", format(group_1_ge10, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
  # Arrow main 1 --> main 2
  geom_segment(x=35, xend=35, y=90, yend=80.5, 
               size=0.15, linejoin = "mitre", lineend = "butt",
               arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  # Arrow exclude 1
  geom_segment(x=35, xend=59.5, y=85, yend=85, 
               size=0.15, linejoin = "mitre", lineend = "butt",
               arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  # Arrow exclude 2
  geom_segment(x=35, xend=59.5, y=65, yend=65, 
               size=0.15, linejoin = "mitre", lineend = "butt",
               arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  # Arrow exclude 3
  geom_segment(x=35, xend=59.5, y=45, yend=45, 
               size=0.15, linejoin = "mitre", lineend = "butt",
               arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  # Arrow main 2 --> group 1
  geom_segment(x=35, xend=35, y=70, yend=60.5, 
               size=0.15, linejoin = "mitre", lineend = "butt",
               arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  # Arrow main 3 --> Final
  geom_segment(x=35, xend=35, y=50, yend=40.5, 
               size=0.15, linejoin = "mitre", lineend = "butt",
               arrow = arrow(length = unit(1, "mm"), type= "closed")) +
  # Remove background
  theme_void()
print(flow)
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_2.png"), device = "png", width = 8, height = 7, units = "in", plot = flow)  

# Footnote - HMDB IDs of excluded metabolites
excluded  <- as.data.frame(rownames(fdata[ which(fdata$pm_pilot_fail_heparin=='1'), ]))
write.xlsx(excluded, file = file.path(main.dir,"final_results/supplemental_figures/supp_figure_2_footnote.xlsx"), overwrite = TRUE)

