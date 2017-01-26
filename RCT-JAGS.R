###############################################################################
###  For license and how to run this script, see 
###  RCT-README.txt
###############################################################################


#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

# Specify name and format for saved files:
fileNameRoot = "RCT-JAGS-" 
graphFileType = "pdf" 

#------------------------------------------------------------------------------- 
#Load The data file 
myDataFrame = read.csv("RCT-DataGen-Data.csv")
myDataFrame$Time = factor( 
  myDataFrame$Time , 
  levels=c("PreTreat","PostTreat","FollowUp1","FollowUp2") )
myDataFrame$Treatment = factor( 
  myDataFrame$Treatment , 
  levels=c("Control","Medication","Counseling") )

# Specify some contrasts of interest:
yName = "Anxiety"
xBetweenName="Treatment" 
xWithinName="Time" 
xSubjectName="Subj"
xBetweenContrasts = list(
  list( c("Medication") , c("Counseling") , 
        #compVal=0.0 , 
        ROPE=c(-2,2) ) ,
  list( c("Control") , c("Medication","Counseling") , 
        #compVal=0.0 , 
        ROPE=c(-2,2) ) 
)
xWithinContrasts = list( 
  list( c("PreTreat") , c("PostTreat","FollowUp1","FollowUp2") , 
        #compVal=0.0 , 
        ROPE=c(-2,2) ) ,
  list( c("PostTreat") , c("FollowUp1","FollowUp2") , 
        #compVal=0.0 , 
        ROPE=c(-2,2) ) 
)
xBetweenWithinContrasts = list(
  list( list( c("Medication") , c("Control") ) ,
        list(  c("PreTreat") , c("PostTreat") ) ,
        #compVal=0.0 , 
        ROPE=c(-2,2) ) ,
  list( list(  c("Counseling") , c("Medication") ) ,
        list(  c("PreTreat") , c("PostTreat") ) ,
        #compVal=0.0 , 
        ROPE=c(-2,2) ) ,
  list( list(  c("Counseling") , c("Medication") ) ,
        list(  c("PostTreat") , c("FollowUp1","FollowUp2") ) ,
        #compVal=0.0 , 
        ROPE=c(-2,2) ) 
)
# Simple contrasts are specified at end of script!

#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
# A variant of Jags-Ymet-XnomSplitPlot-MnormalHom.R from DBDA2E.
source("RCT-JAGS-functions.R") 

#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
mcmcCoda = genMCMC( datFrm=myDataFrame , 
                    yName=yName , xBetweenName=xBetweenName , 
                    xWithinName=xWithinName , xSubjectName=xSubjectName ,
                    numSavedSteps=20000 , thinSteps=20 , 
                    saveName=fileNameRoot )

#-------------------------------------------------------------------------------
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda)
show( parameterNames ) # show all parameter names, for reference
for ( parName in c("mSxW[1,1]","b0","bB[1]","bW[1]","bS[1]","bBxW[1,1]",
                   "sigma","sigmaB","sigmaW","sigmaS","sigmaBxW") ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName ,
            saveName=fileNameRoot , saveType=graphFileType )
}

#-------------------------------------------------------------------------------
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda ,
                        datFrm=myDataFrame , xBetweenName=xBetweenName ,
                        xWithinName=xWithinName , xSubjectName=xSubjectName ,
                        xBetweenContrasts=xBetweenContrasts ,
                        xWithinContrasts=xWithinContrasts ,
                        xBetweenWithinContrasts=xBetweenWithinContrasts ,
                        saveName=fileNameRoot )
show(summaryInfo)
# Show contrast info:
show( tail( summaryInfo , 7 )[,c("Mode","HDIlow","HDIhigh","ESS")] )

#------------------------------------------------------------------------------- 
# Display posterior information:
plotMCMC( mcmcCoda , 
          datFrm=myDataFrame , yName=yName , xBetweenName=xBetweenName , 
          xWithinName=xWithinName , xSubjectName=xSubjectName ,
          xBetweenContrasts=xBetweenContrasts , 
          xWithinContrasts=xWithinContrasts , 
          xBetweenWithinContrasts=xBetweenWithinContrasts ,
          saveName=fileNameRoot , saveType=graphFileType )

#------------------------------------------------------------------------------- 
# Make graphs in style of Cumming (2014) and do "simple" contrasts mentioned in
# Cumming (2014):
#.............................................................................
# Build a matrix of MCMC chains for cell means:
mcmcMat = as.matrix(mcmcCoda,chains=TRUE)
datFrm=myDataFrame
xBetweenLvl = levels(as.factor(datFrm[,xBetweenName]))
xWithinLvl = levels(as.factor(datFrm[,xWithinName]))
cellMean=NULL
for ( Bidx in 1:length(xBetweenLvl) ) {
  for ( Widx in 1:length(xWithinLvl) ) {
    m = ( mcmcMat[,"b0"]
          + mcmcMat[,paste0("bB[",Bidx,"]")]
          + mcmcMat[,paste0("bW[",Widx,"]")]
          + mcmcMat[,paste0("bBxW[",Bidx,",",Widx,"]")] )
    cellMean = cbind( cellMean , m )
    colnames(cellMean)[ncol(cellMean)] = paste0("m[",Bidx,",",Widx,"]") 
  }
}

#.............................................................................
# Make graph of cell means comparable to Cumming (2014)
hdiLow = hdiHigh = meanMean = rep(0,ncol(cellMean))
names(hdiLow) = names(hdiHigh) = names(meanMean) = colnames(cellMean)
for ( cIdx in 1:ncol(cellMean) ) {
  meanMean[cIdx] = mean( cellMean[,cIdx] )
  hdi = HDIofMCMC( cellMean[,cIdx] )
  hdiLow[cIdx] = hdi[1]
  hdiHigh[cIdx] = hdi[2]
}
openGraph(width=6,height=5)
par( mar=c(4,4,2,1) , mgp=c(2.0,0.7,0) )
plot(-10,-10, 
     xlim=c(0.7,length(xWithinLvl)+0.2) , 
     xlab=paste(xWithinName) , 
     xaxt="n" , ylab=yName ,
     ylim=c(min(hdiLow),max(hdiHigh)+10) ,
     main="Post. Mean with 95% HDI", cex.main=1.5 , cex.lab=1.5 )
axis( 1 , at=1:length(xWithinLvl) , tick=FALSE , 
      lab=paste( xWithinLvl  ) )
# Plot cell means connected by lines, and HDIs:
for ( bIdx in 1:length(xBetweenLvl) ) {
  xOffset = -0.1*(bIdx-2)
  cVec = grep( paste0("m\\[",bIdx) , colnames(cellMean) )
  lines( 1:length(xWithinLvl)+xOffset , meanMean[cVec] , type="o" , 
         pch=c(15:19)[bIdx] , cex=2.0 )
  segments( x0=1:length(xWithinLvl)+xOffset , y0=hdiLow[cVec]+xOffset ,
            x1=1:length(xWithinLvl)+xOffset , y1=hdiHigh[cVec]+xOffset ,
            lwd=3 , lend=1 )
}
legend( x="topright" , legend=xBetweenLvl , cex=1.5 ,
        pch=c(15:19)[1:length(xBetweenLvl)] , pt.cex=2.0 )
saveGraph( file=paste0(fileNameRoot,"CellMeans") , 
           type=graphFileType )

#.............................................................................
# Check that averages of means give same results as plotMCMC:

ControlIdx = which( levels( datFrm[,xBetweenName] ) == "Control" )
MedicationIdx = which( levels( datFrm[,xBetweenName] ) == "Medication" )
CounselingIdx = which( levels( datFrm[,xBetweenName] ) == "Counseling" )
PreIdx = which( levels( datFrm[,xWithinName] ) == "PreTreat" )
PostIdx = which( levels( datFrm[,xWithinName] ) == "PostTreat" )
Follow1Idx = which( levels( datFrm[,xWithinName] ) == "FollowUp1" )
Follow2Idx = which( levels( datFrm[,xWithinName] ) == "FollowUp2" )

openGraph(height=4*.8,width=4)
plotPost( 
  cellMean[,c( 
    paste0("m[",ControlIdx,",",Follow1Idx,"]") #,
  ) ] 
  - 
    cellMean[,c( 
      paste0("m[",MedicationIdx,",",Follow1Idx,"]") #,
    ) ] 
  ,
  main="Control - Medication @ FollowUp1" ,  xlab="Difference" ,
  ROPE=c(-2,2) , showCurve=TRUE )
saveGraph( file=paste0(fileNameRoot,"Control-Medication@FollowUp1") , 
           type=graphFileType )

openGraph(height=4*.8,width=4)
plotPost( 
  cellMean[,c( 
    paste0("m[",ControlIdx,",",PreIdx,"]") # ,
  ) ] 
  - 
    cellMean[,c( 
      paste0("m[",ControlIdx,",",PostIdx,"]") #,
    ) ] 
  ,
  main="PreTreat - PostTreat @ Control" ,  xlab="Difference" ,
  ROPE=c(-2,2) , showCurve=TRUE )
saveGraph( file=paste0(fileNameRoot,"PreTreat-PostTreat@Control") , 
           type=graphFileType )
