##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_figure_1.R

# Programmer: Emma McGee 

# Date Last Updated: February 8, 2023

# Purpose of Program: Flow diagram illustrating selection of participants (Supplemental Figure 1)

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

# Set up empty grid to plot flow diagram
data <- tibble(x= 1:110, y= 1:110)

data %>% 
  ggplot(aes(x, y)) +
  scale_x_continuous(minor_breaks = seq(10, 110, 10)) +
  scale_y_continuous(minor_breaks = seq(10, 110, 10)) +
  theme_linedraw() ->
  grid

# Load covariate data from blood questionnaire for calculations
nhs1.cov <- read.csv(file = file.path(main.dir, "data/NHS1_covariate_data_blood_cohort.csv"), sep = ",", header = T, na.strings = "")

# Set N for main box #1 - provided a blood sample
main_1 <- dim(nhs1.cov)[1]
main_1

# Set N for exclusion box #1 - prior history of cancer at blood collection
table(nhs1.cov$canhx)
exclude_1 <- dim(nhs1.cov[(nhs1.cov$canhx == "1"),])[1]
exclude_1

# Set N for main box #2 - cancer free at blood collection
main_2 <- main_1 - exclude_1
main_2

# Set N for exclusion box #2 - not selected into case-control study or a potentially eligible Black woman
load(file.path(main.dir, "data/merged.eset.raw.RData"))
caco_n <- length(merged.eset.raw$id)
exclude_2 <- main_2 - caco_n
exclude_2

# Set N for main box #3 - selected into case-control study or potentially eligible Black woman
main_3 <- main_2 - exclude_2
main_3

# Set N for exclusion box #3 - race other than Black or White
merged.eset.raw <- merged.eset.raw[ , !(merged.eset.raw$nhs1_race == "1" | merged.eset.raw$nhs1_race == "2")]
exclude_3 <- length(merged.eset.raw$id)
exclude_3

# Set N for main box #3 - all eligible and included study participants
main_4 <- main_3 - exclude_3
main_4

# Set N for group #1 - Black women
load(file.path(main.dir, "data/known.eset.final.RData"))
group_1 <- length(known.eset.final$race_bw[(known.eset.final$race_bw == "1")])
group_1

# Set N for group #2 - White women
group_2 <- length(known.eset.final$race_bw[(known.eset.final$race_bw == "0")])
group_2

#### 2 - CREATE FIGURE ####

# Create flow diagram
flow <- grid +
              # Main box 1: NHS blood cohort
              geom_rect(xmin = 16, xmax=60, ymin=90, ymax=100, color='black',
                        fill='white', size=0.25) +
              annotate('text', x= 38, y=95,label= paste("Provided blood sample\n",
                        "(n=", format(main_1, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
              # Main box 2: breast cancer free at blood collection
              geom_rect(xmin = 16, xmax=60, ymin=70, ymax=80, color='black',
                        fill='white', size=0.25) +
              annotate('text', x= 38, y=75,label= paste('Cancer-free at blood collection\n',
                        "(n=", format(main_2, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
              # Main box 3: selected into a case-control study or eligible to be selected
              geom_rect(xmin = 16, xmax=60, ymin=50, ymax=60, color='black',
                        fill='white', size=0.25) +
              annotate('text', x= 38, y=55,label= paste('Underwent metabolomic profiling\n',
                        "for a prior study* (n=", format(main_3, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
              # Main box 4: eligible and included in the analysis
              geom_rect(xmin = 16, xmax=60, ymin=30, ymax=40, color='black',
                        fill='white', size=0.25) +
              annotate('text', x= 38, y=35,label= paste('Self-identified as Black or White\n',
                        "(n=", format(main_4, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
              # Exclude box 1: prior cancer
              geom_rect(xmin = 66, xmax=110, ymin=80, ymax=90, color='black',
                      fill='white', size=0.25) +
              annotate('text', x= 88, y=85,label= paste("Excluded: prior history of cancer\n",
                        "(n=", format(exclude_1, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +            
              # Exclude box 2: not selected into case-control study of Black woman
              geom_rect(xmin = 66, xmax=110, ymin=60, ymax=70, color='black',
                        fill='white', size=0.25) +
              annotate('text', x= 88, y=65,label= paste("Excluded: did not undergo metabolomic\n",
                        "profiling for a prior study* (n=", format(exclude_2, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
              # Exclude box 3: not Black or White
              geom_rect(xmin = 66, xmax=110, ymin=40, ymax=50, color='black',
                        fill='white', size=0.25) +
              annotate('text', x= 88, y=45,label= paste("Excluded: did not self-identify as Black\n",
                        "or White (n=", format(exclude_3, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
              # Group box 1: Black women
              geom_rect(xmin = 1, xmax=35, ymin=10, ymax=20, color='black',
                        fill='white', size=0.25) +
              annotate('text', x= 18, y=15,label= paste('Black women\n',
                        "(n=", format(group_1, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
              # Group box 2: White women
              geom_rect(xmin = 41, xmax=75, ymin=10, ymax=20, color='black',
                        fill='white', size=0.25) +
              annotate('text', x= 58, y=15,label= paste('White women\n',
                        "(n=", format(group_2, big.mark=",", scientific=FALSE), ")", sep=""), size=3.4) +
              # Arrow main 1 --> main 2
              geom_segment(x=38, xend=38, y=90, yend=80.5, 
                          size=0.15, linejoin = "mitre", lineend = "butt",
                          arrow = arrow(length = unit(1, "mm"), type= "closed")) +
              # Arrow main 2 --> main 3
              geom_segment(x=38, xend=38, y=70, yend=60.5, 
                          size=0.15, linejoin = "mitre", lineend = "butt",
                          arrow = arrow(length = unit(1, "mm"), type= "closed")) +
              # Arrow main 3 --> main 4
              geom_segment(x=38, xend=38, y=50, yend=40.5, 
                          size=0.15, linejoin = "mitre", lineend = "butt",
                          arrow = arrow(length = unit(1, "mm"), type= "closed")) +
              # Arrow exclude 1
              geom_segment(x=38, xend=65.5, y=85, yend=85, 
                         size=0.15, linejoin = "mitre", lineend = "butt",
                         arrow = arrow(length = unit(1, "mm"), type= "closed")) +
              # Arrow exclude 2
              geom_segment(x=38, xend=65.5, y=65, yend=65, 
                          size=0.15, linejoin = "mitre", lineend = "butt",
                          arrow = arrow(length = unit(1, "mm"), type= "closed")) +
              # Arrow exclude 3
              geom_segment(x=38, xend=65.5, y=45, yend=45, 
                          size=0.15, linejoin = "mitre", lineend = "butt",
                          arrow = arrow(length = unit(1, "mm"), type= "closed")) +
              # Arrow groups box 1 & 2
              geom_segment(x=38, xend=38, y=30, yend=25, 
                          size=0.15, linejoin = "mitre", lineend = "butt") +
              geom_segment(x=18, xend=38, y=25, yend=25, 
                          size=0.15, linejoin = "mitre", lineend = "butt") +
              geom_segment(x=38, xend=58, y=25, yend=25, 
                          size=0.15, linejoin = "mitre", lineend = "butt") +
              geom_segment(x=18, xend=18, y=25, yend=20.5, 
                          size=0.15, linejoin = "mitre", lineend = "butt",
                          arrow = arrow(length = unit(1, "mm"), type= "closed")) +
              geom_segment(x=58, xend=58, y=25, yend=20.5, 
                          size=0.15, linejoin = "mitre", lineend = "butt",
                          arrow = arrow(length = unit(1, "mm"), type= "closed")) +
              # Remove background
              theme_void()
print(flow)
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_1.png"), device = "png", width = 8, height = 7, units = "in", plot = flow)  
