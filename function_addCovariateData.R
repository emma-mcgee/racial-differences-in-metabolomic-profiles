##########################################################################

# Title: Add Covariate Data to Metabolomics ESet

##########################################################################

# Program: function_addCovariateData.R

# Programmer: Emma McGee 
              # Adapted from code written by Oana Zeleznik

# Date Last Updated: July 27, 2021

# Purpose of Program: Create a function to add covariate data to metabolomics ESet (addCovariateData)

# Statistical Analyses: None

# Additional Study Information: See all_data.R

##########################################################################

addCovariateData = function(covariate.data = NULL, eset.data = NULL, prefix = NULL, verbose = F){
  
  ##### eset.data:
    # eset to be annotated with the covariate data
  
  ##### covariate data:
    # data to be merged to the sample annotation data in eset.data
    # it has to contain a column named "id" where the sample ids are stored as integers
  
  ##### prefix
    # prefix to add to covariate names to know from which file these came from (to find data dictionary easier)
  
  ##### verbose
    # if true will print out information during the execution of the function
  
  # Print warnings if non-ESet object, no id variable, or missing parameters
  if(class(eset.data) != "ExpressionSet"){
    cat("ERROR: the parameter eset.data has to have the type Expression Set.\n")
    return(0)
  }
  
  if(! "id" %in% colnames(covariate.data)){
    cat("ERROR: the parameter covariate.data has to include a column with the name id.\nThese ids will be used to match the covariate data to the Expression Set data.\n")
    return(0)
  }
  
  if(is.null(covariate.data) | is.null(eset.data) | is.null(prefix)){
    cat("ERROR: Plese provide ALL the following parameters: covariate.data, eset.data, prefix.\n")
    return(0)
  }
  
  if(verbose) cat("Adding", dim(covariate.data)[2]-1, "covariates with prefix \"",prefix,"\". \n")
  
  # Extract data from ESet
  metabolites.values = exprs(eset.data)
  samples.annot = pData(eset.data)
  samples.annot$id = as.integer(samples.annot$id)
  
  # Check that IDs are overlapping between ESet and covariates
  if(length(intersect(samples.annot$id, covariate.data$id)) == 0){
    print("ERROR: The IDs in your Expression Set are not annotated in the covariate data you provided.\n")
    return(0)
  }
  
  # Set index for ID variable
  id.index = which(colnames(covariate.data) == "id")
  
  # Add prefix to the new covariates (excluding ID variable)
  colnames(covariate.data)[-1*id.index] = paste(prefix, colnames(covariate.data)[-1*id.index], sep = "")
  
  # Merge Eset annotations with covariates - maintain order and original samples from samples.annot pData usng left_join instead of merge
  samples.annot.new = left_join(x = samples.annot, y = covariate.data, by = "id", copy=FALSE)
  
  # Set rownames (NOTE: original row order is maintained in samples.annot.new)
  rownames(samples.annot.new) = rownames(samples.annot)
  
  # Check for IDs that are not annotated in covariate data
  ids.not.annotated = setdiff(samples.annot$id, covariate.data$id)
  
  if(length(ids.not.annotated) >=1){
    if(verbose) cat("The following IDs are not annotated in the covariate data:\n", 
                    paste(ids.not.annotated, sep = " ", collapse = ", "),".\n",
                    "The not annotated samples will be set to NA in the corresponding variables.\n")
  }
  
  identical(rownames(metabolites.values), rownames(fData(eset.data)))
  identical(colnames(metabolites.values), rownames(samples.annot.new))
  
  
  # Create new ESet
  new.eset = ExpressionSet(assayData = metabolites.values,
                           phenoData = AnnotatedDataFrame(samples.annot.new), 
                           featureData = AnnotatedDataFrame(fData(eset.data)))
  return(new.eset)
  
}