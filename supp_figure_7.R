##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_figure_7.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: quantify observed differences in individual unknown metabolites/peaks between Black and White women (Supplemental Figure 7)

# Statistical Analyses:
    # Linear regression with each ln Z-score transformed metabolite as the outcome variable
    # Restricted to unknown metabolites/peaks with <10% missing data

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####
main.dir <- ""
results.dir <- file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(tidyverse)
library(Biobase)
library(ggplot2)

# Load eset of unknown metabolites/peaks
load(file.path(main.dir, "data/unknown.eset.final.RData"))

# Prepare data
pdata <- pData(unknown.eset.final)
met.vals <- exprs(unknown.eset.final)
dim(met.vals)

# Restrict to metabolites with <10% missing data
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
met.vals <- met.vals[which(missing < 10),]
dim(met.vals)

# Prepare regression output files
lin.res.unknown <- matrix(0,ncol = 4,nrow = dim(met.vals)[1])
colnames(lin.res.unknown) <- c("p.value","effect.estimate","L_CI","U_CI")
rownames(lin.res.unknown) <- rownames(met.vals)


#### 2 - QUANTIFY OBSERVED DIFFERENCES USING LINEAR REGRESSION ####

# Run linear regressions to estimate observed differences in metabolite levels (accounting for matching factors using IPW + adjusted for age)
for(i in 1:dim(met.vals)[1]){
  curr.data <- cbind(pdata[,c("race_bw",
                              "id",
                              "knn_nhs1_agebld","knn_nhs1_agebld_sqr","knn_nhs1_menopmh_bld","knn_nhs1_fast","knn_nhs1_bldmonth",
                              "knn_nhs1_btimelr")],met = met.vals[i,])
  
  # NHS: adjusted for age
  my.lm.model0 <- lm(met~race_bw+
                      knn_nhs1_agebld + knn_nhs1_agebld_sqr,
                      data = curr.data)
  my.model.model0 <- summary(my.lm.model0)      
  my.ci.model0 <- confint(my.lm.model0)
  
  lin.res.unknown[i,] <- c(my.model.model0$coefficients["race_bw1",c("Pr(>|t|)","Estimate")],my.ci.model0["race_bw1",])
}


#### 3 - CREATE FIGURE ####

results.unknown <- as.data.frame(lin.res.unknown)
summary(results.unknown$L_CI)

# Order by effect size
results.unknown <- arrange(results.unknown, effect.estimate)

# Create sequence numbers for plotting
results.unknown$seq <- seq(1, dim(results.unknown)[1], by=1)

# Create plot
barplot<-ggplot(data=results.unknown,
                aes(x=seq, y=effect.estimate, ymin=L_CI, ymax=U_CI)) +
  geom_bar(stat="identity", fill="dodgerblue4") +
  geom_errorbar(width=0.05, size=0.05, color="black") +
  scale_y_continuous(limits=c(-1.35,1.35), breaks = c(-1,-0.5,0,0.5,1),
                     name="Mean difference") +
  geom_hline(yintercept=0, color="black", linetype="solid") +
  labs(x = "Metabolomic features", y = "Mean difference") +
  theme_minimal() +
  theme( 
    legend.position = "right",
    legend.text = element_text(color = "black", size = 24),
    legend.title = element_text(color = "black", size = 30, face="bold"),
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
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_7.png"), device = "png", width = 20, height = 10, units = "in", plot = barplot)  
