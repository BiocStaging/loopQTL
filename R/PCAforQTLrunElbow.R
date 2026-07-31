# Modified based on https://github.com/heatherjzhou/PCAForQTL/blob/b8f7d09e77cb5b50dd712ad99f3fbccab899b73d/R/22.01.04_main1.1_runElbow.R (retrieved 2026/07/31) - license: GPLv3

# Implements the elbow method for choosing K in PCA (based on distance to the diagonal line)

#Given X, the data matrix (must be observation by feature),
#and/or prcompResult, the output from running prcomp() on X,
#select the number of PCs by maximizing the distance to the diagonal line (see our 2022 paper or below for details).
#Either X or prcompResult must be provided. If both are provided, then only prcompResult will be used.
#We recommend only providing prcompResult because it's faster.
.PCAforQTLrunElbow<-function(X=NULL,
                   prcompResult=NULL){

  #Obtain prcompResult.
  if(is.null(prcompResult)){
    if(is.null(X)){
      stop("Please input X or prcompResult.")
    }else{
      cat("Running PCA...\n")
      prcompResult<-prcomp(X,center=TRUE,scale.=TRUE) #This may take a moment.
    }
  }else{
    if(class(prcompResult)!="prcomp"){
      stop("prcompResult must be a prcomp object returned by the function prcomp().")
    }
  }

  importanceTable<-summary(prcompResult)$importance
  x<-1:ncol(importanceTable) #PC indices.
  y<-importanceTable[2,] #PVEs.

  #Given x  and y, calculate the distance between each point and the diagonal line (the line connecting the first and last points).
  x1<-x[1] #First point.
  y1<-y[1] #First point.

  x2<-x[length(x)] #Last point.
  y2<-y[length(y)] #Last point.

  x0<-x
  y0<-y

  distancesDenominator<-sqrt((x2-x1)^2+(y2-y1)^2)
  distancesNumerator<-abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1))
  distances<-distancesNumerator/distancesDenominator

  numOfPCsChosen<-which.max(distances) #12.
  names(numOfPCsChosen)<-NULL
  return(numOfPCsChosen)
}
