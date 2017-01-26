###############################################################################
###  For license and how to run this script, see 
###  RCT-README.txt
###############################################################################

# This program was supposed to generate data consistent with Figure 4 of Cumming
# (2014), which he refers to as a randomized control trial, and is specifically 
# a split-plot design. The means are from Cumming's ESCI-For-Excel-07-10 
# package, file: ESCI_chapters_14-15_Jul_4_2011.xlsm, sheet: Figure. But it
# turns out that generating such data is impossible if a reasonable constraint
# is placed on the pre-treatment scores, namely, that they are all above the
# clinical baseline for anxiety disorder. An email reply from Geoff Cumming, 25
# April 2015: "I'm embarrassed to say that I invented the figure without
# generating data then calculating means and CIs. Your requirement of all
# pretest values above 60 is perfectly reasonable, and I can see that that has
# large implications for variances, and thus all that follows.So you've blown my
# cover totally! And I hadn't realised the implications of the short cut I took.
# Oops. My apologies." So, this program will use different means and variances
# to make a similar point.

# This program is based on SplitPlotAgriData.R

graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

## Means used by Cumming 2014. Within-subject factor has four levels:
## PreTest, PostTest, FollowUp1, FollowUp2. Between subject factor has up to
## four levels: Treatment, Control, Treatment2, Treatment3. Here is the matrix
## of means:
# bBnames = c( "Treatment" , "Control" , "Treatment2" , "Treatment3" ) 
# bWnames = c( "PreTest" , "PostTest" , "FollowUp1" , "FollowUp2" )
# rctMeans = matrix( 
#   dimnames = list( Condition = bBnames , Time = bWnames ) ,
#   ncol=4 , nrow=4 , byrow=TRUE ,
#   data = c( 80 , 45 , 47 , 37 ,
#             87 , 80 , 72 , 84 ,
#             74 , 57 , 58 , 64 ,
#             92 , 93 , 70 , 67 )
# )
# # use only Treatment and Control:
# rctMeans = rctMeans[c("Treatment","Control"),]

bBnames = c( "Control" , "Medication" , "Counseling" ) 
bWnames = c( "PreTreat" , "PostTreat" , "FollowUp1" , "FollowUp2" )
rctMeans = matrix( 
  dimnames = list( Treatment = bBnames , Time = bWnames ) ,
  ncol=length(bWnames) , nrow=length(bBnames) , byrow=TRUE ,
  data = c( 80 , 78 , 76 , 74 ,
            80 , 58 , 72 , 72 ,
            80 , 55 , 58 , 58 ) 
)

b0 = mean(rctMeans)
bW = colMeans(rctMeans) - b0
bB = rowMeans(rctMeans) - b0
bBxW = rctMeans - ( b0 + outer(bB,bW,"+") )

set.seed(47405)
nS = c(15,12,13)*2 # number subjects in each level of bB
Svec = rep(1:length(bB),times=nS)

sigmaS = 9/2 # standard deviation of subject variance
bS = NULL
for ( Bidx in 1:length(bB) ) {
  bSb = rnorm( nS[Bidx] )
  bSb = (bSb-mean(bSb))/sd(bSb) * sigmaS
  bS = c( bS , bSb )
}

sigmaNoise = sigmaS*2

withinCellSigma = sqrt( sigmaS^2 + sigmaNoise^2 )
show( paste0("Within cell sigma = ",withinCellSigma) )

Subj = NULL
Treatment = NULL
Time = NULL
Anxiety = NULL
Sidx = 0
for ( Bidx in Svec ) {
  Sidx = Sidx+1
  for ( Widx in 1:length(bW) ) {
    Subj = c( Subj , Sidx )
    Treatment = c( Treatment , bBnames[Bidx] )
    Time = c( Time , bWnames[Widx] )
    Anxiety = c( Anxiety , 
               b0 + bB[Bidx] + bW[Widx] + bBxW[Bidx,Widx] + bS[Sidx] 
               + rnorm(1,0,sigmaNoise) )
  }
}

myDataFrame = data.frame( Subj=Subj , Treatment=Treatment , Time=Time , 
                          Anxiety=round(Anxiety) )
write.csv( myDataFrame , file="RCT-DataGen-Data.csv", 
           row.names=FALSE )

