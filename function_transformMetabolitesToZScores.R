##########################################################################

# Title: Z-score Transform Metabolites in an ESet

##########################################################################

# Program: function_transformMetabolitesToZScores.R

# Programmer: Emma McGee
              # Adapted from code written by Oana Zeleznik

# Date Last Updated: July 27, 2021

# Purpose of Program: Create a function to z-score transform metabolites in an ESet (transformMetabolitesToZScores)

# Statistical Analyses: None

# Additional Study Information: See all_data.R

##########################################################################

transformMetabolitesToZScores = function(eset){
  
  my.met.exp = exprs(eset)
  
  z.scores.met = t(apply(my.met.exp,1,FUN = function(x) return((x-mean(x, na.rm = T))/sd(x, na.rm = T))))
  
  result.eset = ExpressionSet(assayData = z.scores.met,
                              phenoData = AnnotatedDataFrame(pData(eset)),
                              featureData = AnnotatedDataFrame(fData(eset)))
}
