##########################################################################

# Title: Differences in Metabolomic Profiles Between Black and White Women

##########################################################################

# Program: figure_2.R

# Programmer: Emma McGee

# Date Last Updated: February 8, 2023

# Purpose of Program: Metabolite set enrichment analysis for known metabolites (Figure 2)

# Statistical Analyses:
    # Metabolite set enrichment analysis based an metabolite classes obtained from the Broad lab
    # Regression coefficients are from linear regression models adjusted for matching factors
    # Multiple testing correction using the Benjamini-Hochberg procedure with a False Discovery Rate of 20%

# Additional Study Information: See all_data.R

##########################################################################

#### 1 - PREPARE DATA ####

rm(list = ls())
main.dir <- ""
results.dir <- file.path(main.dir,"final_results")
setwd(main.dir)

# Load packages
library(grid)

# Source observed and residual disparity estimates
source(file = file.path(main.dir,"supp_table_4.R"))


#### 2 - PLOT RESULTS ####

# Format results for plotting
fgsea.list <- list(fgsea.nhs, fgsea.whi, fgsea.nhs.resid)
data.nes <- lapply(fgsea.list, function(x) x[,c(1,5)])
nes <- lapply(data.nes, function(x) melt(x))

# Create results for each study
# NHS
nes.nhs <- nes[[1]]
nes.nhs$Pathway <- factor(nes.nhs$pathway, levels = data.nes[[1]]$pathway)
nes.nhs$study <- "NHS"
range(nes.nhs$value)

# WHI
nes.whi <- nes[[2]]
nes.whi$Pathway <- factor(nes.whi$pathway, levels = data.nes[[2]]$pathway)
nes.whi$study <- "WHI"
range(nes.whi$value)

# NHS - residual
nes.nhs.resid <- nes[[3]]
nes.nhs.resid$Pathway <- factor(nes.nhs.resid$pathway, levels = data.nes[[3]]$pathway)
nes.nhs.resid$study <- "NHS "
range(nes.nhs.resid$value)

# Combine all
nes.all <- rbind(nes.nhs, nes.whi, nes.nhs.resid)

# Re-level study indicator
nes.all <- mutate(nes.all, study=factor(study),
                      study = fct_relevel(study,
                                          "NHS",
                                          "WHI",
                                          "NHS "))

# Define color palette
hm.palette <- rev(brewer.pal(11, "RdBu"))

# Create plot
p <- ggplot(data = nes.all, aes(x = variable, y = Pathway, fill = value)) +
            geom_tile(color = "white", size = 0.1) + 
            scale_fill_gradientn(colours = hm.palette, limits = c(-2.91,2.91))+
            labs(x = "Normalized Enrichment Score (NES)", y = "", fill="NES") +
            facet_grid(cols=vars(study), scales= "free", space="free") +
            guides(fill=guide_colorbar(title = "NES",title.position = "top",nrow = 1, barheight=10, 
                                       barwidth=1, label.position = "right"))+
            theme_minimal()+ 
                theme(axis.text = element_text(angle = 0, vjust = 0.5, 
                                 size = 11, face = "bold", colour = "black"),
                axis.title.x = element_text(face="bold", size=11, color="black"),
                axis.text.x = element_blank(),
                legend.position = "right",
                legend.title = element_text(face = "bold", size = 11, color = "black"),
                legend.text = element_text(face = "bold", size = 11, color = "black"),
                strip.text.x = element_text(size = 11))

# Plot the figure
grid.draw(p)

# Save the figure
ggsave(filename=file.path(results.dir,"figures/figure_2.png"), 
       plot = p, device = "png", dpi = 300, width = 10, height = 5.4)  
