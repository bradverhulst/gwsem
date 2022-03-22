library(gwsem)


# We cannot share real genetic data (but there are lots of places where you can get some).
# Therefore, for the tutorial, we are going to use simulated data.

# There are 2000 SNPs and 6000 people

source("GW-SEMtutorialDataSim.R")

#We have simulated 2 "hits" for each trait

hits            # [1] 1703  450
hits_tob        # [1] 1104 1149
hits_can        # [1] 1347 1231
hits_alc        # [1] 1390  884


# Read the phenotype data into R and look at the data
phenoData <- read.table("oneFacphenoData.txt", header=TRUE)
head(phenoData)

# This would be a great time to recode the data, transform it, and tell R if we have ordinal or binary indicators using mxFactor()
# The data for this example are simulated and continuous, and therefore, we will not be doing anything now, but if you would like 
# to chop the indicators up into binary or ordinal variables, this would be a great time to do it.


# The first thing that we would want to do is build a general model. 
# This model takes a series of arguments: 

                                                                                     # You must tell GW-SEM:
addFac <- buildOneFac(phenoData = phenoData,                                         # what the data object is (which you read in above)
                     itemNames = c("tobacco", "cannabis", "alcohol"),                # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                    # what covariates that you want to include in the analysis
                     fitfun = "WLS")                                                 # and the fit function that you would like to use (WLS is much faster than ML)

# You can take this object (addFac) which is technically an OpenMx model, and simply fit it using mxRun() 
addFacFit <- mxRun(addFac)
summary(addFacFit)

# This is strongly advised, as it is a great time to learn if your model is being specified in a way that you didn't expect or that 
# the model is not giving you the answers that you are expecting

# Now, provided that the model looks reasonable, you can plug the model that you have built into the GWAS function
                                                                                     # To fit this function, you must tell GW-SEM :
GWAS(model = addFac,                                                                 # what model object you would like to fit
	snpData = 'example.pgen',                                                        # that path to the snpData file.
	out="latFac.log",                                                                # the file that you would like to save the full results into
	SNP=1:2000
	)                                                                       # the index of the snps (how many) you would like to fit

# Note about snpData: The path to your snpData will likely include switching to a different directory (as you will likely do your analysis in a different
# folder than your SNP data). All you need to do is point to the data using relative paths. Further, it is able to take plink bed/bim/fam or pgen/psam/pvar 
# data or bgen data (oxford format)

# Note about SNP: This can be used to run a limited number of SNP (i.e. not the whole snp file). This is particularly useful if you would like to run chop up a 
# chr into several parts without cutting you actual genotype data into seperate files.

# This will take a while and frequently be done on a cluster.  If you chop up this script and remove the comments and insert you data, you can use it on your cluster


# The next step is to read the data in. While you are unlikely to want to read all of the output into R, it is useful to do this on a small test GWAS analysis
# to ensure that you are getting reasonable results.

FullResult <- read.delim("latFac.log")                            # This is didactically useful, but it contains much more information than most people want
head(FullResult) 

# Then, we can read the results into R for a specific parameter. This function takes takes two arguments: the path to the data and the column in the results 
# file for the parameter that you want to examine.

succinct <- loadResults(path = "latFac.log", focus =  "snp_to_F")
succinct <- signif(succinct, focus =  "snp_to_F")

# Then, you can plot the results
plot(-log10(succinct$P))

# If this was real data, we would plot a Manhattan Plot, but with simulated data it doesn't really work like that.
# here is some useful data in case you want to do that
#png("my_Manhattan_plot.png", width=1500, height=750)
#par(cex.lab = 2, mai = c(1, 1, .1, .1) + 0.1, bg="white")
#manhattan(myResults, p = "P")
#dev.off()



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

hits            # [2] 1703 (.15)  [5]   450 (.12)
hits_tob        # [1] 1104 (.13)  [4]  1149 (.10)
hits_can        # [6] 1347 (.11)  [8]  1231 (.09)
hits_alc        # [3] 1390 (.14)  [7]   884 (.12)



### The Residuals Model

# The process is very similar for the Latent factor model except that we are regressing the individual items on the SNP

# Again, the first thing that we would want to do is build a general model. 
# The function tells GW-SEM that you want to run the residuals model. This model takes a series of arguments: 
                                                                                       # You must tell GW-SEM:
addFacRes <- buildOneFacRes(phenoData,                                                 # what the data object is (which you read in above)
                     c("tobacco", "cannabis", "alcohol"),                              # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                      # what covariates that you want to include in the analysis
                     fitfun = "WLS")                                                   # and the fit function that you would like to use (WLS is much faster than ML)


# Again, You can take this object (addFac) which is technically an OpenMx model, and simply fit it using mxRun() 
addFacResFit <- mxRun(addFacRes)
summary(addFacResFit)


# Then you will use the same function that we used for the latent variable model, with a different model object

GWAS(model = addFacRes, snpData = 'example.pgen', out="facRes.log")                 

# It is instructive to look at the results again, but this might be easier to do in linux rather than R
# This is how you load it into R
ResResult <- read.delim("facRes.log")                            
head(ResResult)

# Now that we have run the residual's model, we have multiple parameters that we want to load into R (Three in this case)
# The function to do this is the same as for one parameter, but you need to do it several times (once for each paramter). This gives several R objects.

succinctTob <- loadResults(path = "facRes.log", focus =  "snp_to_tobacco")
succinctTob <- signif(succinctTob, , focus =  "snp_to_tobacco")
succinctCan <- loadResults(path = "facRes.log", focus =  "snp_to_cannabis")
succinctCan <- signif(succinctCan, focus =  "snp_to_cannabis")
succinctAlc <- loadResults(path = "facRes.log", focus =  "snp_to_alcohol")
succinctAlc <- signif(succinctAlc, "snp_to_alcohol")
# Now we can plot all the residual manhattan plots


par(mfrow = c(3, 1))
plot(-log10(succinctTob$P), main = "Tobacco")
plot(-log10(succinctCan$P), main = "Cannabis")
plot(-log10(succinctAlc$P), main = "Alcohol")


# hits            # [2] 1703 (.15)  [5]   450 (.12)
# hits_tob        # [1] 1104 (.13)  [4]  1149 (.10)
# hits_can        # [6] 1347 (.11)  [8]  1231 (.09)
# hits_alc        # [3] 1390 (.14)  [7]   884 (.12)

head(succinctTob[order(succinctTob$Z, decreasing = T),], 10)   

#      MxComputeLoop1 CHR   BP     SNP A1 A2 statusCode catch1 snp_to_tobacco Vsnp_to_tobacco:snp_to_tobacco        Z            P
#T  1:           1104   1 1103 snp1103  A  B         OK     NA     0.14900067                   0.0005031415 6.642677 3.080357e-11
#F  2:           1703   1 1702 snp1702  B  A         OK     NA     0.12801211                   0.0004915599 5.773815 7.749655e-09
#A  3:           1390   1 1389 snp1389  A  B         OK     NA     0.12872182                   0.0005046529 5.730015 1.004216e-08
#T  4:           1149   1 1148 snp1148  A  B         OK     NA     0.12203242                   0.0004850003 5.541205 3.003975e-08
#F  5:            450   1  449  snp449  B  A         OK     NA     0.11785668                   0.0004916381 5.315345 1.064554e-07
#C  6:           1347   1 1346 snp1346  A  B         OK     NA     0.11413135                   0.0004936741 5.136707 2.795941e-07
#A  7:            884   1  883  snp883  A  B         OK     NA     0.09218042                   0.0004982755 4.129561 3.634567e-05
#C  8:           1231   1 1230 snp1230  A  B         OK     NA     0.07483854                   0.0005013840 3.342259 8.309954e-04
#   9:            240   1  239  snp239  B  A         OK     NA     0.06124436                   0.0005023504 2.732516 6.285262e-03
#  10:            416   1  415  snp415  B  A         OK     NA     0.06060229                   0.0004945189 2.725195 6.426355e-03

head(succinctCan[order(succinctCan$Z, decreasing = T),], 10)

#      MxComputeLoop1 CHR   BP     SNP A1 A2 statusCode catch1 snp_to_cannabis Vsnp_to_cannabis:snp_to_cannabis        Z            P
#F   1:           1703   1 1702 snp1702  B  A         OK     NA      0.13426467                     0.0004993773 6.008241 1.875465e-09
#T   2:           1104   1 1103 snp1103  A  B         OK     NA      0.12940679                     0.0004986048 5.795339 6.818328e-09
#C   3:           1347   1 1346 snp1346  A  B         OK     NA      0.12604411                     0.0004941459 5.670155 1.426681e-08
#A   4:           1390   1 1389 snp1389  A  B         OK     NA      0.12868773                     0.0005172070 5.658547 1.526598e-08
#T   5:           1149   1 1148 snp1148  A  B         OK     NA      0.11607221                     0.0004938467 5.223146 1.759088e-07
#F   6:            450   1  449  snp449  B  A         OK     NA      0.11476962                     0.0004938006 5.164772 2.407319e-07
#A   7:            884   1  883  snp883  A  B         OK     NA      0.09491221                     0.0005078154 4.211813 2.533288e-05
#C   8:           1231   1 1230 snp1230  A  B         OK     NA      0.08519646                     0.0004954528 3.827546 1.294274e-04
#    9:           1262   1 1261 snp1261  B  A         OK     NA      0.06574975                     0.0004907133 2.968111 2.996358e-03
#   10:            319   1  318  snp318  B  A         OK     NA      0.06413396                     0.0005018366 2.862905 4.197770e-03

head(succinctAlc[order(succinctAlc$Z, decreasing = T),], 10)

#      MxComputeLoop1 CHR   BP     SNP A1 A2 statusCode catch1 snp_to_alcohol Vsnp_to_alcohol:snp_to_alcohol        Z            P
#F   1:            450   1  449  snp449  B  A         OK     NA     0.13780440                   0.0004910780 6.218531 5.018295e-10
#T   2:           1149   1 1148 snp1148  A  B         OK     NA     0.13296277                   0.0004830044 6.049988 1.448567e-09
#A   3:           1390   1 1389 snp1389  A  B         OK     NA     0.13267422                   0.0005146149 5.848512 4.959911e-09
#C   4:           1347   1 1346 snp1346  A  B         OK     NA     0.12392506                   0.0004974231 5.556434 2.753417e-08
#F   5:           1703   1 1702 snp1702  B  A         OK     NA     0.10720084                   0.0004901968 4.841868 1.286243e-06
#A   6:            884   1  883  snp883  A  B         OK     NA     0.10409697                   0.0005041460 4.636176 3.549139e-06
#T   7:           1104   1 1103 snp1103  A  B         OK     NA     0.09829096                   0.0004906294 4.437484 9.101664e-06
#    8:            904   1  903  snp903  B  A         OK     NA     0.06736399                   0.0005032610 3.002833 2.674794e-03
#    9:           1896   1 1895 snp1895  A  B         OK     NA     0.06456991                   0.0004907543 2.914729 3.559980e-03
#   10:           1370   1 1369 snp1369  A  B         OK     NA     0.06346269                   0.0004919247 2.861338 4.218574e-03


# Additional considerations: 

# Because we are working with numerous chromosomes, which we have chopped up into different chunks, we have several log files for each chromosome. For example, 
# we may have a set of log files such as: fac1a.log, fac1b.log, fac1c.log, fac1d.log, fac1e.log, fac1f.log, fac1g.log, fac1h.log.
# To simplify things, we can construct a list of file names for each chr with the following code.

c1  <- paste0(paste0("fac1",  letters)[1:8], ".log")
c2  <- paste0(paste0("fac2",  letters)[1:9], ".log")
c3  <- paste0(paste0("fac3",  letters)[1:7], ".log")
c4  <- paste0(paste0("fac4",  letters)[1:8], ".log")
c5  <- paste0(paste0("fac5",  letters)[1:7], ".log")
c6  <- paste0(paste0("fac6",  letters)[1:7], ".log")
c7  <- paste0(paste0("fac7",  letters)[1:6], ".log")
c8  <- paste0(paste0("fac8",  letters)[1:6], ".log")
c9  <- paste0(paste0("fac9",  letters)[1:5], ".log")
c10 <- paste0(paste0("fac10", letters)[1:5], ".log") 
c11 <- paste0(paste0("fac11", letters)[1:5], ".log")
c12 <- paste0(paste0("fac12", letters)[1:5], ".log")
c13 <- paste0(paste0("fac13", letters)[1:4], ".log")
c14 <- paste0(paste0("fac14", letters)[1:4], ".log")
c15 <- paste0(paste0("fac15", letters)[1:3], ".log")
c16 <- paste0(paste0("fac16", letters)[1:3], ".log")
c17 <- paste0(paste0("fac17", letters)[1:3], ".log")
c18 <- paste0(paste0("fac18", letters)[1:3], ".log")
c19 <- paste0(paste0("fac19", letters)[1:3], ".log")
c20 <- paste0(paste0("fac20", letters)[1:2], ".log")
c21 <- paste0(paste0("fac21", letters)[1:2], ".log")
c22 <- paste0(paste0("fac22", letters)[1:2], ".log")

# Then we can plug these files into the loadResults function so that all the data is loaded into a single object.

res <- loadResults(c(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,
	                 c11, c12, c13, c14, c15, c16, c17, c18, c19, c20,
					 c21, c22), "snp_to_F")



