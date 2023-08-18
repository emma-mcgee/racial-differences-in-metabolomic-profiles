##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: supp_figure_5.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Create matrices of Spearman correlation coefficients for known metabolites, with separate matrices by race (Supplemental Figure 5)

# Statistical Analyses:
    # Spearman correlation coefficients

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####

rm(list = ls())
main.dir = ""
results.dir = file.path(main.dir,"final_results")
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

# Load known metabolite data
load(file.path(main.dir, "data/known.eset.final.RData"))

# Remove metabolites with >=10% missing values
met.vals <- exprs(known.eset.final)
missing <- apply(met.vals, 1, function(x) round(length(which(is.na(x)))/length(x)*100,2))
miss.names <- row.names(met.vals[which(missing >= 10),])
source(file.path(main.dir,"function_removeMetabolitesFromESet.R"))
known.eset.final = removeMetabolitesFromESet(eset = known.eset.final, metabolites.to.remove = miss.names)

# Extract expression, phenotype, and feature data
vals <- exprs(known.eset.final)
pdata <- pData(known.eset.final)
fdata <- fData(known.eset.final)


#### 2 - CREATE CORRELATION MATRIX FOR BLACK WOMEN ####

# Subset metabolite values to Black women
black.val <- as.data.frame(vals[,which(pdata$race_bw == "1")])

# Add metabolite class
black.val$class_broad <- fdata$class_broad

# Recode metabolite classes
black.val <- mutate_at(black.val, c("class_broad") , ~replace(., is.na(.), "Missing"))
black.val <- mutate(black.val, class= factor(class_broad),
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
black.val <- select(black.val, -c(class_broad))
table(black.val$class)

# Order metabolites by class for correlation matrix
black.val.ordered <- black.val[order(black.val$class),]

# Calculate Spearman correlations
cor.mat <- cor(t(black.val.ordered[,1:(dim(black.val.ordered)[2]-1)]), method ="spearman", use="p")

## Plot correlation matrix ##

# Format data for plotting
to.plot <- melt(cor.mat)
to.plot$X1 <- factor(to.plot[,1], levels = rownames(black.val.ordered)[1:dim(black.val.ordered)[1]])
to.plot$X2 <- factor(to.plot[,2], levels = rownames(black.val.ordered)[dim(black.val.ordered)[1]:1])

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

# Export figure - black women
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_5a_black.png"), device = "png", width = 7.5, height = 6.56, units = "in", plot = p)  


#### 3 - CREATE CORRELATION MATRIX FOR WHITE WOMEN ####

# Subset metabolite values to White women
white.val <- as.data.frame(vals[,which(pdata$race_bw == "0")])

# Add metabolite class
white.val$class_broad <- fdata$class_broad

# Recode metabolite classes
white.val <- mutate_at(white.val, c("class_broad") , ~replace(., is.na(.), "Missing"))
white.val <- mutate(white.val, class= factor(class_broad),
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
white.val <- select(white.val, -c(class_broad))
table(white.val$class)

# Order metabolites by class for correlation matrix
white.val.ordered <- white.val[order(white.val$class),]

# Calculate Spearman correlations
cor.mat <- cor(t(white.val.ordered[,1:(dim(white.val.ordered)[2]-1)]), method ="spearman", use="p")

## Plot correlation matrix ##

# Format data for plotting
to.plot <- melt(cor.mat)
to.plot$X1 <- factor(to.plot[,1], levels = rownames(white.val.ordered)[1:dim(white.val.ordered)[1]])
to.plot$X2 <- factor(to.plot[,2], levels = rownames(white.val.ordered)[dim(white.val.ordered)[1]:1])

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

# Export figure - White women
ggsave(filename=file.path(results.dir,"supplemental_figures/supp_figure_5b_white.png"), device = "png", width = 7.5, height = 6.56, units = "in", plot = p)  


