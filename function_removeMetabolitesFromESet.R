##########################################################################

# Title: Remove Metabolites From a Metabolomics ESet

##########################################################################

# Program: function_removeMetabolitesFromESet.R

# Programmer: Emma McGee 
              # Adapted from code written by Oana Zeleznik

# Date Last Updated: July 27, 2021

# Purpose of Program: Create a function to remove metabolites from a metabolomics ESet (removeMetabolitesFromESet)

# Statistical Analyses: None

# Additional Study Information: See all_data.R

##########################################################################

removeMetabolitesFromESet = function(eset, metabolites.to.remove){
  
  ##### eset:
    # eset from which metabolites are to be removed
  
  ##### metabolites.to.remove:
    # list of rownames of metabolites to be removed
  
  # Extract data
  eset.fdata = fData(eset)
  eset.exprs = exprs(eset)
  
  # Remove metabolites
  met.to.remove = intersect(rownames(eset.fdata), metabolites.to.remove)
  help.mets = setdiff(metabolites.to.remove, rownames(eset.fdata))
  
  cat("Starting to remove ", length(metabolites.to.remove), " metabolites.\n")
  if(length(help.mets)>0)
    cat("Cannot remove the following metabolites as these have not been found in the data set:",help.mets,"\n")
  
  eset.fdata = subset(eset.fdata, !rownames(eset.fdata) %in% met.to.remove)
  eset.exprs = subset(eset.exprs, !rownames(eset.exprs) %in% met.to.remove)
  
  # Save new ESet
  new.eset = ExpressionSet(assayData = eset.exprs,
                           phenoData = AnnotatedDataFrame(pData(eset)), 
                           featureData = AnnotatedDataFrame(eset.fdata))
  
  cat("Dimension of new eset:", dim(new.eset),"\n")
  return(new.eset)
}