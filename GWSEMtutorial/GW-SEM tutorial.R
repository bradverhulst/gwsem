# The goal of the tutorial is to teach you how to use the GW-SEM software with very little guidance.
library(gwsem)

# The line below will simulate phenotypic data using 2000 SNPs for 6000 people
source("GW-SEMtutorialDataSim.R")

# Read the phenotype data into R and look at the data
phenoData <- read.table("oneFacphenoData.txt", header=TRUE)

# The first thing you should do is build the GWAS model. 

                                                  
FacMod <- buildOneFac(                # what the data dataframe called?
                                      # what are the items called?
                                      # what covariates will you want to include?
                      )               # use WLS (it is much faster than ML)!

# Use mxRun() to make sure that the model runs. 
# What are the parameter estimates?
# What does the factor model look like? Are there any issues?




# Provided that the model looks reasonable, plug the model that you have built into the GWAS function
GWAS(  ,              # what is the model object?
	   ,              # what is the path to the snpData file?
	   ,              # Where would you like to save the results?
	   )              # Which snps would you like to fit?

# Read the results into R 
# replace the XXX with an appropriate name

XXX <- loadResults()
XXX <- signif()

# Plot the results
plot(-log10()


# What are the Most significant SNPs?

head(succinct[order(succinct$Z, decreasing = T),], 10)

#     MxComputeLoop1 CHR   BP     SNP A1 A2 statusCode catch1   snp_to_F Vsnp_to_F:snp_to_F        Z            P
#  1:           1104   1 1103 snp1103  A  B         OK     NA 0.16651323       0.0005721841 6.961148 3.375113e-12
#  2:           1703   1 1702 snp1702  B  A         OK     NA 0.15197069       0.0005636670 6.401011 1.543519e-10
#  3:           1390   1 1389 snp1389  A  B         OK     NA 0.15442192       0.0005830522 6.395216 1.603206e-10
#  4:           1149   1 1148 snp1148  A  B         OK     NA 0.14511672       0.0005640739 6.110116 9.955879e-10
#  5:            450   1  449  snp449  B  A         OK     NA 0.14306095       0.0005704739 5.989674 2.102616e-09
#  6:           1347   1 1346 snp1346  A  B         OK     NA 0.14142158       0.0005689580 5.928920 3.049333e-09
#  7:            884   1  883  snp883  A  B         OK     NA 0.11231702       0.0005716700 4.697566 2.632805e-06
#  8:           1231   1 1230 snp1230  A  B         OK     NA 0.08927399       0.0005648222 3.756376 1.723916e-04
#  9:           1883   1 1882 snp1882  A  B         OK     NA 0.06359136       0.0005636696 2.678464 7.396072e-03
# 10:           1118   1 1117 snp1117  A  B         OK     NA 0.06008892       0.0005556720 2.549090 1.080045e-02




### Run the Residuals Model

                                                                       
FacRes <- buildOneFacRes(         ,       # what is data dataframe called?  
                                  ,       # what are the items called?  
                                  ,       # what covariates will you want to include?  
                     fitfun = "WLS")      # use WLS (it is much faster than ML)!  


# Make sure the model is reasonable by fitting it using mxRun() 
XXX <- mxRun(XXX)
summary(xxx)


# Run the gwas model
GWAS()                 

# It is instructive to look at the results again, but this might be easier to do in linux rather than R
# This is how you load it into R
XXX <- read.delim()                            
head(XXX)

# Load the various parameters into R

XXX <- loadResults()
XXX <- signif()
XXX <- loadResults()
XXX <- signif()
XXX <- loadResults()
XXX <- signif()
# Now we can plot all the residual manhattan plots


par(mfrow = c(3, 1))
plot(-log10(succinctTob$P), main = "Tobacco")
plot(-log10(succinctCan$P), main = "Cannabis")
plot(-log10(succinctAlc$P), main = "Alcohol")



# Look at the data from the residuals model
head(succinctTob[order(succinctTob$Z, decreasing = T),], 10)   
head(succinctCan[order(succinctCan$Z, decreasing = T),], 10)
head(succinctAlc[order(succinctAlc$Z, decreasing = T),], 10)



## Don't look down here. This is where the answers are.

# hits            # [1] 1242  417
# hits_tob        # [1] 292 451
# hits_can        # [1] 1860 1301
# hits_alc        # [1] 1486 1397