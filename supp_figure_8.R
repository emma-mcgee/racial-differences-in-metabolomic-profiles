##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_figure_8.R

# Programmer: Emma McGee 

# Date Last Updated: February 8, 2023

# Purpose of Program: create matrix of Spearman correlation coefficients between known and unknown metabolites with observed differences >= |0.5| (Supplemental Figure 8)

# Statistical Analyses:
    # Linear regression with each ln Z-score transformed metabolite as the outcome variable
    # Restricted to unknown metabolites/peaks with <10% missing data

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####
rm(list = ls())
main.dir <- ""
results.dir <- file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(gplots)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(reshape)
library(Biobase)
library(tidyverse)

# Source observed difference estimates for known metabolites
source(file = file.path(main.dir,"supp_table_2.R"))

# Source observed difference estimates for unknown metabolites/peaks
source(file = file.path(main.dir,"supp_figure_7.R"))

# Load known metabolite data and subset to observed differences >= |0.5|
source(file.path(main.dir,"function_removeMetabolitesFromESet.R"))
load(file.path(main.dir, "data/known.eset.final.RData"))
met.vals <- exprs(known.eset.final)
missing <- apply(met.vals, 1, FUN = function(x) return(length(which(is.na(x)))/length(x)*100))
missing.names <- rownames(met.vals)[which(missing >= 10)]
known.eset.miss <- removeMetabolitesFromESet(eset = known.eset.final, metabolites.to.remove = missing.names)
known.names <- rownames(results)[which(abs(as.numeric(results$effect.estimate.main)) < 0.5)]
known.eset.0.5 <- removeMetabolitesFromESet(eset = known.eset.miss, metabolites.to.remove = known.names)

# Load unknown metabolite data and subset to observed differences >= |0.5|
load(file.path(main.dir, "data/unknown.eset.final.RData"))
met.vals <- exprs(unknown.eset.final)
missing <- apply(met.vals, 1, FUN = function(x) return(length(which(is.na(x)))/length(x)*100))
missing.names <- rownames(met.vals)[which(missing >= 10)]
unknown.eset.miss <- removeMetabolitesFromESet(eset = unknown.eset.final, metabolites.to.remove = missing.names)
unknown.names <- rownames(results.unknown)[which(abs(results.unknown$effect.estimate) < 0.5)]
unknown.eset.0.5 <- removeMetabolitesFromESet(eset = unknown.eset.miss, metabolites.to.remove = unknown.names)

# Extract expression data and merge
vals.0.5.known <- exprs(known.eset.0.5)
vals.0.5.unknown <- exprs(unknown.eset.0.5)
vals.0.5 <- as.data.frame(rbind(vals.0.5.known, vals.0.5.unknown))


#### 2 - CREATE CORRELATION MATRIX ####

# Add metabolite class
fdata.known <- fData(known.eset.0.5)
fdata.known <- mutate_at(fdata.known, c("class_broad") , ~replace(., is.na(.), "Unclassified"))
fdata.unknown <- fData(unknown.eset.0.5)
fdata.unknown <- mutate_at(fdata.unknown, c("class_broad") , ~replace(., is.na(.), "Unknown"))
fdata <- rbind(fdata.known, fdata.unknown)
vals.0.5$class_broad <- fdata$class_broad


# Recode metabolite classes
table(vals.0.5$class_broad, useNA="always")
vals.0.5 <- mutate(vals.0.5, class= factor(class_broad),
                    class = fct_recode(class,
                                       "Carboxylic acids and derivatives" = "Carboxylic acids and derivatives",
                                       "Cholesteryl esters" = "Cholesteryl esters",
                                       "Diglycerides" = "Diglycerides",
                                       "Organic acids and derivatives" = "Organic acids and derivatives",
                                       "Organoheterocyclic compounds" = "Organoheterocyclic compounds",
                                       "Phosphatidylcholine plasmalogens" = "Phosphatidylcholine plasmalogens",
                                       "Phosphatidylcholine plasmalogens" = "PC plasmalogens",
                                       "Phosphatidylcholines" = "Phosphatidylcholines",
                                       "Triglycerides" = "Triglycerides",
                                       "Other" = "Alkaloids and derivatives",
                                       "Other" = "Lysophosphatidylcholines",
                                       "Other" = "Lysophosphatidylethanolamines",
                                       "Other" = "Phosphatidylethanolamine plasmalogens",
                                       "Other" = "Phosphatidylethanolamines",
                                       "Unclassified" = "Unclassified",
                                       "Unknown" = "Unknown"),
                    class = fct_relevel(class,
                                        "Carboxylic acids and derivatives",
                                        "Cholesteryl esters",
                                        "Diglycerides",
                                        "Organic acids and derivatives",
                                        "Organoheterocyclic compounds",
                                        "Phosphatidylcholine plasmalogens",
                                        "Phosphatidylcholines", 
                                        "Triglycerides",
                                        "Other",
                                        "Unclassified",
                                        "Unknown"))
table(vals.0.5$class, useNA="always")
vals.0.5 <- select(vals.0.5, -c(class_broad))

# Order metabolites by class for correlation matrix
vals.0.5.ordered <- vals.0.5[order(vals.0.5$class),]

# Calculate Spearman correlations
cor.mat <- cor(t(vals.0.5.ordered[,1:(dim(vals.0.5.ordered)[2]-1)]), method ="spearman", use="p")


## Plot correlation matrix ##

# Format data for plotting
to.plot <- melt(cor.mat)
to.plot$X1 <- factor(to.plot[,1], levels = rownames(vals.0.5.ordered)[1:dim(vals.0.5.ordered)[1]])
to.plot$X2 <- factor(to.plot[,2], levels = rownames(vals.0.5.ordered)[dim(vals.0.5.ordered)[1]:1])

# Create color palette
hm.palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")), bias = 1)

# Create plot
p <- ggplot(data = to.plot, aes(x = X1, y = X2, fill = value)) + 
  geom_tile(color = "white", size = 0.1) + 
  scale_x_discrete(position = "top") +
  scale_fill_gradientn(colours = hm.palette(11), limits = c(-1, 1), 
                       breaks = c(-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1), 
                       labels = c("–1.0","–0.8","–0.6","–0.4","–0.2", "0","0.2","0.4","0.6","0.8","1.0"))+
  theme_minimal()+ 
  labs(x = "", y= "", fill = "")+
  guides(fill = guide_colorbar(label=T, barheight=29, barwidth=0.5, ticks=TRUE, ticks.colour="Black",
                               title.position="right"))+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, 
                               size = 2, hjust = 0, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 2),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12, angle = 90),
    legend.text = element_text(face = "bold", size = 12))

# Plot figure
plot(p)

# Export figure
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_8.png"), device = "png", width = 7.5, height = 6.56, units = "in", plot = p)  

