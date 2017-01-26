###############################################################################
###  For license and how to run this script, see 
###  Jags-MetaAnalysis-TwoProportions-README.txt
###############################################################################

# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

# The following file must be in R's current working directory:
source("Jags-MetaAnalysis-TwoProportions-Utilities.R") 

fileNameRoot="Jags-MetaAnalysis-TwoProportions" # For output file names.

#------------------------------------------------------------------------------
# Specify the data. 

dataSource = c( "TowelReuse" , "BetaBlocker" )[2]

if ( dataSource=="TowelReuse" ) {
  fileNameRoot=paste0(fileNameRoot,"-TowelReuse")
  # Towel re-use data are from Benjamin Scheibehenne, Tahira Jamil, and Eric-Jan
  # Wagenmakers (2016). Bayesian evidence synthesis can reconcile seemingly 
  # inconsistent results: The case of hotel towel reuse. Psychological Science, 
  # vol. 27 no. 7, pp. 1043-1046. doi 10.1177/0956797616644081
  theData = matrix( c(
    28,   2, 101,  31,
    123,  24, 472, 104,
    21,   4,  21,   3,
    82, 105, 278, 277,
    77,  58, 406, 249,
    103, 174, 587, 731,
    74, 137,  98, 124 ) , ncol=4, byrow=TRUE )
  colnames( theData ) = c("controlReuse", "controlThrow", "normReuse", "normThrow" )
  theData = cbind( StudyID = 1:nrow(theData) , theData )
  theData = cbind( theData , 
                   controlTotal=theData[,"controlReuse"]+theData[,"controlThrow"] )
  theData = cbind( theData , 
                   normTotal=theData[,"normReuse"]+theData[,"normThrow"] )
  theData = data.frame( theData )
  # Package the data for JAGS:
  dataList = list(
    nS = length(theData$StudyID) ,
    zC = theData$controlReuse ,
    nC = theData$controlTotal , 
    zT = theData$normReuse ,
    nT = theData$normTotal 
  )
}

if ( dataSource=="BetaBlocker" ) {
  fileNameRoot=paste0(fileNameRoot,"-BetaBlocker")
  # Data from Yusuf et al. (1985), as reported in Table 5.4, p. 124, of Gelman, 
  # Carlin, Stern, Dunson, Vehtari, and Rubin (2014). Bayesian Data Analysis,
  # Third Edition. Boca Raton: CRC Press.
  # http://www.stat.columbia.edu/~gelman/book/data/meta.asc
  theData = matrix( c(
     1,     3, 39,        3, 38,
     2,     14, 116,      7, 114, 
     3,     11, 93,         5, 69, 
     4,     127, 1520,  102, 1533, 
     5,     27, 365,      28, 355, 
     6,     6, 52,          4, 59, 
     7,     152, 939,     98, 945, 
     8,     48, 471,      60, 632, 
     9,     37, 282,      25, 278, 
     10,    188, 1921,  138, 1916, 
     11,    52, 583,      64, 873, 
     12,    47, 266,      45, 263, 
     13,    16, 293,       9, 291, 
     14,    45, 883,      57, 858, 
     15,    31, 147,      25, 154, 
     16,    38, 213,      33, 207, 
     17,    12, 122,      28, 251, 
     18,    6, 154,        8, 151, 
     19,    3, 134,        6, 174, 
     20,    40, 218,      32, 209, 
     21,    43, 364,      27, 391, 
     22,    39, 674,      22, 680 
  ) , ncol=5, byrow=TRUE )
  colnames(theData) = c( "StudyID" , "ControlDeaths" , "ControlTotal" , 
                       "TreatedDeaths" , "TreatedTotal" )
  theData = data.frame(theData)
  # Package the data for JAGS:
  dataList = list(
    nS = length(theData$StudyID) ,
    zC = theData$ControlDeaths ,
    nC = theData$ControlTotal , 
    zT = theData$TreatedDeaths ,
    nT = theData$TreatedTotal 
  )
}

#------------------------------------------------------------------------------

# Define the JAGS model:
#
# thetaC[s] is success probability in control group, study s.
# thetaT[s] is success probability in treatment group, study s.
# rho[s] is the difference of log-odds between groups:
#   rho[s] = logit(thetaT[s]) - logit(thetaC[s])
#   Re-arranged, the equation expresses the relation of thetaT[s] to thetaC[s]:
#   thetaT[s] = logistic( rho[s] + logit(thetaC[s]) )
#   This relation is a "natural" way to represent the dependency of the 
#   probabilities between groups because the relation is (i) symmetric with
#   respect to what outcome is defined as success because logit(.5+p) =
#   -logit(.5-p), and (ii) symetric with respect to which group is defined as
#   the treatment by reversing the sign of rho. 
#   Note that rho[s] is the so-called log odds ratio across groups: 
#   rho[s] = log( (thetaT[s]/(1-thetaT[s])) / (thetaC[s]/(1-thetaC[s])) )
# thetaComega is modal thetaC[s] across control studies.
# thetaCkappa is concentration (consistency) of thetaC[s] across studies.
# rhoMu is the mean rho[s] across studies.
# rhoSD is the standard deviation rho[s] across studies.
#
modelString = "
  model {
    for ( s in 1:nS ) {
      zC[s] ~ dbin( thetaC[s] , nC[s] )
      zT[s] ~ dbin( thetaT[s] , nT[s] )
      thetaC[s] ~ dbeta( thetaComega*(thetaCkappa-2)+1 ,
                         (1-thetaComega)*(thetaCkappa-2)+1 )
      thetaT[s] <- ilogit( rho[s] + logit( thetaC[s] ) ) # ilogit is logistic
      rho[s] ~ dnorm( rhoMu , 1/rhoSD^2 )
    }
    # Prior rhoMeta, rhoSD:
    rhoMu ~ dnorm( 0 , 1/10^2 )
    rhoSD ~ dgamma(1.64,0.64) # mode=1,sd=2
    # Prior on thetaComega and thetaCkappa:
    thetaComega ~ dbeta(1.01,1.01) 
    thetaCkappa <- thetaCkappaMinusTwo + 2
    thetaCkappaMinusTwo ~ dgamma(2.618,0.162) # mode=10 , sd=10
    # Derived variables:
    for ( s in 1:nS ) { 
      mu[s] <- thetaT[s] / thetaC[s]  # risk ratio
    }
    thetaTomega <- ilogit( rhoMu + logit( thetaComega ) )
    muMeta <- thetaTomega / thetaComega  # risk ratio
  }
" # close quote for modelString
writeLines( modelString , con="TEMPmodel.txt" )

# Variables to monitor:
parameters = c( "thetaC" , "thetaComega" , "thetaCkappa" , 
                "thetaT" , "thetaTomega" ,
                "rho" , "rhoMu" , "rhoSD" , 
                "mu" , "muMeta" )

# Initial values for MCMC chains:
initsList = list( 
  thetaC=rep( sum(dataList$zC)/sum(dataList$nC) , dataList$nS ) ,
  thetaComega=sum(dataList$zC)/sum(dataList$nC) , 
  rho = rep( 0 , dataList$nS ) ,
  rhoMu = 0 ,
  rhoSD = 1 ,
  thetaCkappaMinusTwo = 5
)

# Run the chains:
adaptSteps = 1000 
burnInSteps = 1000
numSavedSteps = 30000
thinSteps=10
nChains = 3

runJagsOut <- run.jags( method=c("rjags","parallel")[2] ,
                        model="TEMPmodel.txt" , 
                        monitor=parameters , 
                        data=dataList ,  
                        inits=initsList , 
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps , 
                        sample=ceiling(numSavedSteps/nChains) ,
                        thin=thinSteps ,
                        summarise=FALSE ,
                        plots=FALSE )
codaSamples = as.mcmc.list( runJagsOut )

save( codaSamples , file=paste0(fileNameRoot,"-Mcmc.Rdata") )
mcmcMat = as.matrix(codaSamples)

# Convergence diagnostics:

# Graphical view of selected variables:
for ( parName in c( "thetaC[1]" , "thetaComega" , "thetaCkappa" , 
                    "thetaT[1]" , "thetaTomega" ,
                    "rho[1]" , "rhoMu" , "rhoSD" ,
                    "mu[1]" , "muMeta"  ) ) {
  diagMCMC( codaObject=codaSamples , parName=parName )
}

#------------------------------------------------------------------------------
# Numerical view of all variables:
# Using runjags package summary method, which includes SSeff and psrf:
( summaryMCMC = summary( runJagsOut )) 
# Any effective samples sizes less than 10,000?
which( summaryMCMC[,"SSeff"] < 10000 )
summaryMCMC[ which( summaryMCMC[,"SSeff"] < 10000 ) , "SSeff" ]
# Any R-hats greater than 1.05?
which( summaryMCMC[,"psrf"] > 1.05 )
summaryMCMC[ which( summaryMCMC[,"psrf"] > 1.05 ) , "psrf" ]

#--------------------------------------------------------------------------
# # Posterior graphs:

# Plot correlation of thetaC and thetaT:
openGraph(height=3.75,width=7)
layout(matrix(1:2,nrow=1,byrow=TRUE))
pltIdx = round(seq(1,nrow(mcmcMat),length=1000))
par( mar=c(3.0,4.0,0.5,0.5) , mgp=c(2.0,0.7,0) , pty="s" )
#
plot( mcmcMat[pltIdx,"thetaC[1]"] , mcmcMat[pltIdx,"thetaT[1]"] , asp=1 ,
      xlab="thetaC[1]" , ylab="thetaT[1]" , main=""  , col="skyblue" )
abline(0,1,lty="dashed")
#
plot( mcmcMat[pltIdx,"thetaComega"] , mcmcMat[pltIdx,"thetaTomega"] , asp=1 ,
      xlab="thetaComega" , ylab="thetaTomega" , 
      main="" , col="skyblue" )
abline(0,1,lty="dashed")
#
graphName = paste0(fileNameRoot,"-Corr")
saveGraph( file=graphName , type="pdf" )
saveGraph( file=graphName , type="eps" )
saveGraph( file=graphName , type="png" )

#--------------------------------------------------------------------------

# Function for plotting a single horizontal violin within a plot:
plotViolin = function( mcmcSamp , plotHt=0 , dataVal=NULL , dataN=500 ,  
                       cenTend=c("mode","median")[1] ,
                       cexMin=0.75 , cexMax=2.0 , threshN=700 , slopeN=6/1000 ,
                       topOnly=FALSE , violinWd=c(0.90,1.50)[topOnly+1] ) {
  # Plot density curve:
  densCurve = density( mcmcSamp )
  polygon( densCurve$x , plotHt+densCurve$y/max(densCurve$y)*(violinWd/2) , 
           col="skyblue" , border="skyblue" )
  if ( !topOnly ) {
    polygon( densCurve$x , plotHt-densCurve$y/max(densCurve$y)*(violinWd/2) , 
             col="skyblue" , border="skyblue" )
  }
  # Plot HDI:
  hdiLim = HDIofMCMC( mcmcSamp )
  lines( hdiLim , rep(plotHt,2) , lend=1 , lwd=4 )
  # Plot central tendency of study sIdx:
  mcmcDensity = density(mcmcSamp)
  mcmcMode = mcmcDensity$x[which.max(mcmcDensity$y)]
  if ( cenTend=="mode" ) { mcmcCenTend = mcmcMode }
  if ( cenTend=="median" ) { mcmcCenTend = median( mcmcSamp ) }
  lines( rep(mcmcCenTend,2) , c(plotHt,plotHt+(violinWd/2)) , 
         lwd=2 , lend=1 )
  if ( !topOnly ) {
    lines( rep(mcmcCenTend,2) , c(plotHt-(violinWd/2),plotHt) , 
           lwd=2 , lend=1 )
  }
  # Display numerical details of mode and HDI:
  text( par('usr')[2] , plotHt , 
        labels=paste0( format(round(mcmcCenTend,2),nsmall=2)
                       ," (", format(round((hdiLim[1]),2),nsmall=2)
                       ,",", format(round((hdiLim[2]),2),nsmall=2) ,")" ) , 
        adj=c(1.1,-0.5) , cex=0.75 )
  # Plot data point of study sIdx:
  if ( !is.null(dataVal) ) {
    if ( !topOnly ) { pChar=23 } else { pChar=24 }
    points(  dataVal , plotHt , 
             pch=pChar , col="black" , lwd=1 , bg="grey" ,
             cex = cexMin+(cexMax-cexMin)/(1+exp(-slopeN*(dataN-threshN))) )
    #text( dataVal , plotHt , bquote(N==.(dataN)) , 
    #      cex=0.5 , adj=c(0.5,2.0) )
    text( dataVal , plotHt , bquote("N="*.(dataN)) , 
          cex=0.5 , adj=c(0.5,1.85) )
  }
}

#--------------------------------------------------------------------------
# Forest plot of mu:

openGraph( height=10 , width=5 )
par( mar=c(3.0,4.0,0.5,0.5) , mgp=c(2.0,0.7,0) )
muLim = quantile( mcmcMat[ , grep( "mu\\[" , colnames(mcmcMat) ) ] , 
                 probs=c(0.001,0.999) )
xLim = c( muLim[1] , muLim[2] + 0.5*(muLim[2]-muLim[1]) )
plot( -1,-1, 
      xlim=xLim ,
      xlab="Multiplier for Treatment" ,
      ylim=c(1,dataList$nS+2.5) , ylab="Study ID Number" ,
      yaxt="n" )
axis( side=2 , at=c(1:dataList$nS,dataList$nS+2) ,
      labels=c(as.character(1:dataList$nS),"Overall") , las=1 )
abline( v=1.0 , lty="dotted" )
for ( sIdx in 1:dataList$nS ) {
  plotViolin( mcmcSamp=mcmcMat[,paste0("mu[",sIdx,"]")] , plotHt=sIdx ,
              dataVal=((dataList$zT[sIdx]/dataList$nT[sIdx])
                        /(dataList$zC[sIdx]/dataList$nC[sIdx])) , 
              dataN=(dataList$nT[sIdx]+dataList$nC[sIdx]) ,
              topOnly=TRUE )
}
# Overall:
plotViolin( mcmcSamp=mcmcMat[,"muMeta"] , plotHt=dataList$nS+2 ,
            topOnly=TRUE ) 
mcmcDensity = density(mcmcMat[,"muMeta"])
mcmcMode = mcmcDensity$x[which.max(mcmcDensity$y)]
#abline( v=mcmcMode , lty="dashed" )
abline(h=sIdx+1)

graphName = paste0(fileNameRoot,"-Forest-Mu")
saveGraph( file=graphName , type="pdf" )
saveGraph( file=graphName , type="eps" )
saveGraph( file=graphName , type="png" )

#-----------------------------------------------------------------------------
# Forest plot of log-odds:

openGraph( height=10 , width=5 )
par( mar=c(3.0,4.0,0.5,0.5) , mgp=c(2.0,0.7,0) )
muLim = quantile( mcmcMat[ , grep( "rho\\[" , colnames(mcmcMat) ) ] , 
                  probs=c(0.001,0.999) )
xLim = c( muLim[1] , muLim[2] + 0.5*(muLim[2]-muLim[1]) )
plot( -1,-1, 
      xlim=c(-1.0,1.5) , xlab="Log Odds Ratio" ,
      ylim=c(1,length(theData$StudyID)+2.5) , ylab="Study ID Number" ,
      yaxt="n" )
axis( side=2 , at=c(1:dataList$nS,dataList$nS+2) ,
      labels=c(as.character(1:dataList$nS),"Overall") , las=1 )
abline( v=0.0 , lty="dotted" )
for ( sIdx in 1:dataList$nS ) {
  plotViolin( mcmcSamp=mcmcMat[,paste0("rho[",sIdx,"]")] , plotHt=sIdx ,
              dataVal=log( (dataList$zT[sIdx]/(dataList$nT[sIdx]-dataList$zT[sIdx]))
                       / (dataList$zC[sIdx]/(dataList$nC[sIdx]-dataList$zC[sIdx])) ) , 
              dataN=(dataList$nT[sIdx]+dataList$nC[sIdx]) ,
              topOnly=TRUE )
}
# Overall:
plotViolin( mcmcSamp=mcmcMat[,"rhoMu"] , plotHt=dataList$nS+2 ,
            topOnly=TRUE ) 
mcmcDensity = density(mcmcMat[,"rhoMu"])
mcmcMode = mcmcDensity$x[which.max(mcmcDensity$y)]
#abline( v=mcmcMode , lty="dashed" )
abline(h=sIdx+1)

graphName = paste0(fileNameRoot,"-Forest-LogOdds")
saveGraph( file=graphName , type="pdf" )
saveGraph( file=graphName , type="eps" )
saveGraph( file=graphName , type="png" )


################################################################################
# Plot thetaC - thetaT, like Figure 3 of Brophy, Joseph, & Rouleau (2001).

openGraph(width=7,height=3.5)
layout(matrix(1:2,nrow=1,byrow=TRUE))
par( mar=c(4.0,4.0,1.5,0.5) , mgp=c(2.0,0.7,0) )
#
smoothScatter( mcmcMat[,"thetaComega"] , mcmcMat[,"thetaTomega"] , 
      xlab=bquote(omega[C]) , ylab=bquote(omega[T]) ,
      asp=1 , col="skyblue" , cex.lab=1.5 )
abline(0,1,lty="dashed")
#abline(-0.02,1)
#
plotPost( mcmcMat[,"thetaTomega"]-mcmcMat[,"thetaComega"] , 
          xlab=bquote(omega[T]-omega[C]) , showCurve=TRUE )
# #
# plot( mcmcMat[,"thetaC[10]"] , mcmcMat[,"thetaT[10]"] , 
#       xlab=bquote(thetaC[10]) , ylab=bquote(thetaT[10]) ,
#       asp=1 , col="skyblue" )
# abline(0,1)
# #abline(-0.02,1)
# #
# plotPost( mcmcMat[,"thetaC[10]"]-mcmcMat[,"thetaT[10]"] , 
#           xlab=bquote(thetaC[10]-thetaT[10]) ,
#           compVal=0.0 #, 
#           #ROPE=c(-0.02,0.02) 
#           )
graphName = paste0(fileNameRoot,"-ControlMinusTreatment")
saveGraph( file=graphName , type="pdf" )


