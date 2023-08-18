##########################################################################

# Title: Natural Log and Z-score Transform Metabolites in an ESet

##########################################################################

# Program: function_transformMetabolitesToLNZScores.R

# Programmer: Emma McGee
              # Adapted from code written by Oana Zeleznik

# Date Last Updated: July 27, 2021

# Purpose of Program: Create a function to natural log and z-score transform metabolites in an ESet (transformMetabolitesToLnZScores)

# Statistical Analyses: None

# Additional Study Information: See all_data.R

##########################################################################


transformMetabolitesToLnZScores = function(eset){

  my.met.exp = exprs(eset)
  
  ln.vals = t(apply(my.met.exp,1,FUN = function(x) return(log(x, base = exp(1)))))
  ln.z.scores.met = t(apply(ln.vals,1,FUN = function(x) return((x-mean(x, na.rm = T))/sd(x, na.rm = T))))
  
  result.eset = ExpressionSet(assayData = ln.z.scores.met,
                              phenoData = AnnotatedDataFrame(pData(eset)),
                              featureData = AnnotatedDataFrame(fData(eset)))
}