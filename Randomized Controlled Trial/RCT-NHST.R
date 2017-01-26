###############################################################################
###  For license and how to run this script, see 
###  RCT-README.txt
###############################################################################

## Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!

library(nlme)

fileNameRoot="RCT-NHST"
graphFileType = "pdf" 
myDataFrame = read.csv("RCT-DataGen-Data.csv")
myDataFrame$Time = factor( 
  myDataFrame$Time , 
  levels=c("PreTreat","PostTreat","FollowUp1","FollowUp2") )
myDataFrame$Treatment = factor( 
  myDataFrame$Treatment , 
  levels=c("Control","Medication","Counseling") )
attach(myDataFrame)

## Set the contrasts
options(contrasts=c("contr.sum", "contr.poly"))

## General, omnibus ANOVA:
# lmeResults = lme( Anxiety ~ Treatment + Time + Treatment:Time , 
#                   random = ~ 1|Subj )
# anova(lmeResults)
# show(lmeResults)

###############
### Panel A ###
###############

## Between-subject factor marginal contrast comparing control to medication and
## counseling
##     Control  Medication  Counseling   
c1 = c(    1    ,    -.5   ,    -.5  )  # Cntl vs med and couseling
c2 = c(    0   ,      -1   ,      1  )  # med vs counsel (not needed here)
K = length(levels( Treatment )) 
cMat = cbind( rep(1/K,K) , c1 ,c2) 
contrasts( Treatment ) = solve(t(cMat))[,-1] # solve() is matrix inverse.
# contrasts( Treatment )
## See the results
panAs = summary( lme( Anxiety ~ Treatment + Time + Treatment:Time , 
                      random = ~ 1|Subj ) )
panelA = intervals( lme( Anxiety ~ Treatment + Time + Treatment:Time , 
                         random = ~ 1|Subj ) )
cat( "Panel A:\n")
show( round( c( 
  panelA$fixed["Treatmentc1",] , 
  p = panAs$tTable["Treatmentc1","p-value"] ) , 3 ) )
## Frequentist results are in Treatmentc1 in summary and interval above. The
## marginal contrasts on between subject factors have considerably higher
## uncertainty in frequentist version.


###############
### Panel B ###
###############

## Within-subject factor contrasts
##      PreTest PostTest FollowUp1 FollowUp2  
c1 = c(     0   ,    1   ,  -.5    ,  -.5  ) #post treatment versus follow ups 1 and 2
c2 = c(    -1   ,   1/3  ,  1/3    ,  1/3  ) #other comparisons we aren't concerned with
c3 = c(     0   ,   0    ,   -1    ,   1   )
K = length(levels( Time )) 
cMat = cbind( rep(1/K,K) , c1 , c2 , c3 ) 
contrasts( Time ) = solve(t(cMat))[,-1] # solve() is matrix inverse.
contrasts( Time )
## See the results:
panBs = summary( lme( Anxiety ~ Treatment + Time + Treatment:Time , 
                      random = ~ 1|Subj ) )
panelB = intervals( lme( Anxiety ~ Treatment + Time + Treatment:Time , 
                         random = ~ 1|Subj ) )
cat( "Panel B:\n")
show( round( c( 
  panelB$fixed["Timec1",] , 
  p = panBs$tTable["Timec1","p-value"] ) , 3 ) )
## Frequentist results are in Timec1 in summary and interval above. The marginal
## contrasts on within subject factors are similar across methods.

######################
### Panels C and D ###
######################

## Between-subject factor contrasts
##     Control  Medication  Counseling   
c1 = c(   -1    ,     1   ,      0  )  # Medication versus control
c2 = c(    0    ,    -1   ,      1  )  # Counseling versus medication
K = length(levels( Treatment )) 
cMat = cbind( rep(1/K,K) , c1 ,c2) 
contrasts( Treatment ) = solve(t(cMat))[,-1] # solve() is matrix inverse.
# contrasts( Treatment )

## Within-subject factor contrasts
##      PreTest PostTest FollowUp1 FollowUp2  
c1 = c(     1   ,    -1    ,    0    ,   0   ) #pre versus post treatment
c2 = c(     0   ,     1    ,  -.5    ,  -.5  ) #post treatment versus follow ups 1 and 2
c3 = c(     0   ,     0    ,   -1    ,   1   ) #other comparison we aren't concerned with
K = length(levels( Time )) 
cMat = cbind( rep(1/K,K) , c1 , c2 , c3 ) 
contrasts( Time ) = solve(t(cMat))[,-1] # solve() is matrix inverse.
# contrasts( Time )

panCDs = summary( lme( Anxiety ~ Treatment + Time + Treatment:Time , 
                       random = ~ 1|Subj ) )
panelCD = intervals( lme( Anxiety ~ Treatment + Time + Treatment:Time , 
                          random = ~ 1|Subj ) )
cat( "Panel C:\n")
show( round( c( 
  panelCD$fixed["Treatmentc1:Timec1",] , 
  p = panCDs$tTable["Treatmentc1:Timec1","p-value"] ) , 3 ) )

cat( "Panel D:\n")
show( round( c( 
  panelCD$fixed["Treatmentc2:Timec2",] , 
  p = panCDs$tTable["Treatmentc2:Timec2","p-value"] ) , 3 ) )
## Panel C, Medication versus Control by Pre versus Post treatment, this is
## Treatmentc1:Timec1.
## Panel D, Counseling versus Medication by Post treatment versus follow ups, 
## this is Treatmentc2:Timec2.

###############
### Panel E ###
###############

## The approach here is to just do an independent or repeated measures t-test on
## a subset of data (no pooling of variance from other cells).

## Extract the data for follow up 1
f1Data = myDataFrame[ myDataFrame[,'Time']=='FollowUp1' , ]
f1Cntrl   = f1Data[ f1Data[,'Treatment']=='Control'    , 'Anxiety' ]
f1Med     = f1Data[ f1Data[,'Treatment']=='Medication' , 'Anxiety' ]
f1Counsel = f1Data[ f1Data[,'Treatment']=='Counseling' , 'Anxiety' ]

## Compare Control versus Medication in first follow up:
panelE = t.test( f1Cntrl , f1Med )
cat( "Panel E:\n")
show( round( c( panelE$conf.int[1] , 
                unname(panelE$estimate[1] - panelE$estimate[2]) ,
                panelE$conf.int[2] , p=panelE$p.value ) , 3 ) )


###############
### Panel F ###
###############

## Extract the data for control
cntrlData = myDataFrame[ myDataFrame[,'Treatment']=='Control' , ]
preCntrl  = cntrlData[ cntrlData[,'Time']=='PreTreat'  , 'Anxiety' ]
postCntrl = cntrlData[ cntrlData[,'Time']=='PostTreat' , 'Anxiety' ]

## compare Pre versus Post in control condition:
panelF = t.test( preCntrl , postCntrl , paired=TRUE )
cat( "Panel F:\n")
show( round( c( panelF$conf.int[1] , unname(panelF$estimate) ,
                panelF$conf.int[2] , p=panelF$p.value ) , 3 ) )

## end .......................................................................