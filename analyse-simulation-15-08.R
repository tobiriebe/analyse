
analyseSimulation <- function(dataFile) {
  
  library(pROC)
  library(Rmisc)
  
  data <- dataFile #save data file for later functions appearing
  
  
  Y <- dataFile$Y                          
  Ynoise <- dataFile$Ynoise
  yBin  <- dataFile$yBin
  yBinNoise <- dataFile$yBinNoise
  originalX <- dataFile$originalX
  originalXnoise <- dataFile$originalXnoise            
  coeffs <- dataFile$coeffs
  predictors <- dataFile$predictors
  kappa <- dataFile$kappa                    
  samples <- dataFile$samples
  simulations <- dataFile$simulations
  nOuterFolds <- dataFile$nOuterFolds               
  redSteps <- dataFile$redSteps
  sampleRatio <- dataFile$sampleRatio
  fixedMstop <- dataFile$fixedMstop                
  fixedNu <- dataFile$fixedNu
  offsetFinal <- dataFile$offsetFinal
  predModelList <- dataFile$predModelList           
  offsetFinalClass <- dataFile$offsetFinalClass
  predModelListClass <- dataFile$predModelListClass
  predictionVector <- dataFile$predictionVector          
  predictionVectorClass <- dataFile$predictionVectorClass
  offsetFinalNoise <- dataFile$offsetFinalNoise
  predModelListNoise <- dataFile$predModelListNoise        
  predictionVectorNoise <- dataFile$predictionVectorNoise
  offsetFinalClassNoise <- dataFile$offsetFinalClassNoise
  predModelListClassNoise <- dataFile$predModelListClassNoise   
  predictionVectorClassNoise <- dataFile$predictionVectorClassNoise
  
  
  
  
  
  #load(dataFile)
  #   if (grepl(pattern = "_Plain", x = dataFile)){
  #     fileType <- fileTypePlain
  #   } else {
  #     fileType <- fileTypeVols
  #   }
  
  if (!exists("redSteps")){
    redSteps <- dim(predictionVector)[2]
  }
  
  
  
  
  #initalize variables for absolute errors
  AbsErrors <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save absolute errors after each simulation
  sumAbsErrors <- rep(0, redSteps) #vector to calculate sum of absolut errors
  AbsErrorsCI97.5 <- rep(0, redSteps)#vectors to save confidence intervals
  AbsErrorsCI2.5 <- rep(0, redSteps)
  AbsErrorsNoiseX <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save absolute errors after each simulation
  sumAbsErrorsNoiseX <- rep(0, redSteps)#vector to calculate sum of absolut errors
  AbsErrorsNoiseXCI97.5 <- rep(0, redSteps)#vectors to save confidence intervals
  AbsErrorsNoiseXCI2.5 <- rep(0, redSteps)
  AbsErrorsNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save absolute errors after each simulation
  sumAbsErrorsNoiseY <- rep(0, redSteps)#vector to calculate sum of absolut errors
  AbsErrorsNoiseYCI97.5 <- rep(0, redSteps)#vectors to save confidence intervals
  AbsErrorsNoiseYCI2.5 <- rep(0, redSteps)
  AbsErrorsNoiseXNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save absolute errors after each simulation
  sumAbsErrorsNoiseXNoiseY <- rep(0, redSteps)#vector to calculate sum of absolut errors
  AbsErrorsNoiseXNoiseYCI97.5 <- rep(0, redSteps) #vectors to save confidence intervals
  AbsErrorsNoiseXNoiseYCI2.5 <- rep(0, redSteps)
  
  #initialize variables for AUC
  AUCsb <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save AUCafter each simulation
  sumAUC <- rep(0, redSteps) #vector for sum of AUC and 
  CIAUC2.5 <- rep(0, redSteps)#vectors to save confidence intervals
  CIAUC97.5 <- rep(0, redSteps)
  AUCNoiseX <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save AUCafter each simulation
  sumAUCNoiseX <- rep(0, redSteps) #vector for sum of AUC and 
  CIAUCNoiseX2.5 <- rep(0, redSteps)#vectors to save confidence intervals
  CIAUCNoiseX97.5 <- rep(0, redSteps)
  AUCNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save AUCafter each simulation
  sumAUCNoiseY <- rep(0, redSteps) #vector for sum of AUC and 
  CIAUCNoiseY2.5 <- rep(0, redSteps)#vectors to save confidence intervals
  CIAUCNoiseY97.5 <- rep(0, redSteps)
  AUCNoiseXNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix to save AUCafter each simulation
  sumAUCNoiseXNoiseY <- rep(0, redSteps) #vector for sum of AUC and 
  CIAUCNoiseXNoiseY2.5<- rep(0, redSteps)#vectors to save confidence intervals
  CIAUCNoiseXNoiseY97.5<- rep(0, redSteps)
  
  
  #initialize variables for variance
  variance <- matrix(0, nrow = redSteps, ncol = simulations)#matrix save variance after each reduction step
  sumvariance <- rep(0, redSteps) #vector for sum of variance
  CIvariance2.5 <- rep(0, redSteps) #vector to save confidence intervals of variance
  CIvariance97.5 <- rep(0, redSteps)
  varianceNoiseX <- matrix(0, nrow = redSteps, ncol = simulations)#matrix save variance after each reduction step
  sumvarianceNoiseX <- rep(0, redSteps) #vector for sum of variance
  CIvariance2.5NoiseX <- rep(0, redSteps) #vector to save confidence intervals of variance
  CIvariance97.5NoiseX <- rep(0, redSteps)
  varianceNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix save variance after each reduction step
  sumvarianceNoiseY <- rep(0, redSteps) #vector for sum of variance
  CIvariance2.5NoiseY <- rep(0, redSteps) #vector to save confidence intervals of variance
  CIvariance97.5NoiseY <- rep(0, redSteps)
  varianceNoiseXNoiseY <- matrix(0, nrow = redSteps, ncol = simulations)#matrix save variance after each reduction step
  sumvarianceNoiseXNoiseY <- rep(0, redSteps) #vector for sum of variance
  CIvariance2.5NoiseXNoiseY <- rep(0, redSteps) #vector to save confidence intervals of variance
  CIvariance97.5NoiseXNoiseY <- rep(0, redSteps)
  
  
  for (simulation in 1:simulations) {
    
    print(simulation)
    ########REGRESSION
    AbsErrors[,simulation] <- colSums(abs(t(predictionVector[simulation,,])-Y[,simulation]))/simulations
    AbsErrorsNoiseX[,simulation] <- colSums(abs(t(predictionVectorNoise[simulation,,])-Y[,simulation]))/simulations
    AbsErrorsNoiseY[,simulation] <- colSums(abs(t(predictionVector[simulation,,])-Ynoise[,simulation]))/simulations
    AbsErrorsNoiseXNoiseY[,simulation] <- colSums(abs(t(predictionVectorNoise[simulation,,])-Ynoise[,simulation]))/simulations
    
    #######CLASSIFICATION (VIA REGRESSION)
    for (reduction in 1:redSteps){
      #print(reduction)#DEBUG
      #sumAUC[reduction] <- sumAUC[reduction] + auc(yBin[,simulation], predictionVectorClass[simulation, reduction, ])
      #calculate AUC
      AUCsb[reduction, simulation] <- auc(yBin[,simulation], predictionVectorClass[simulation, reduction, ])
      AUCNoiseX[reduction, simulation] <- auc(yBin[,simulation], predictionVectorClassNoise[simulation, reduction, ])
      AUCNoiseY[reduction, simulation] <- auc(yBinNoise[,simulation], predictionVectorClass[simulation, reduction, ])
      AUCNoiseXNoiseY[reduction, simulation] <- auc(yBinNoise[,simulation], predictionVectorClassNoise[simulation, reduction, ])
      
      #calculate variance 
      variance[reduction, simulation] <- var(predictionVector[simulation , reduction ,])
      varianceNoiseX[reduction, simulation] <- var(predictionVectorNoise[simulation , reduction ,])
      varianceNoiseY[reduction, simulation] <- var(predictionVector[simulation , reduction ,] - Ynoise[,simulation])
      varianceNoiseXNoiseY[reduction, simulation] <- var(predictionVectorNoise[simulation , reduction ,] - Ynoise[,simulation] )
      
      
      #print(variance2.5) #DEBUG
    }
    
  } #end simulation loop
  
  
  
  for (reduction in 1:redSteps){
    #calculate confidence intervals for sum of absolut errors
    AbsErrorsCI <- CI(AbsErrors[reduction,]) #take all reductions steps over the simulations and calculate CI
    AbsErrorsCI97.5[reduction] <- AbsErrorsCI[1] #upper CI
    AbsErrorsCI2.5[reduction] <- AbsErrorsCI[3] #upper CI
    AbsErrorsCI <- CI(AbsErrorsNoiseX[reduction,]) #take all reductions steps over the simulations and calculate CI
    AbsErrorsNoiseXCI97.5[reduction] <- AbsErrorsCI[1] #upper CI
    AbsErrorsNoiseXCI2.5[reduction] <- AbsErrorsCI[3] #upper CI
    AbsErrorsCI <- CI(AbsErrorsNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    AbsErrorsNoiseYCI97.5[reduction] <- AbsErrorsCI[1] #upper CI
    AbsErrorsNoiseYCI2.5[reduction] <- AbsErrorsCI[3] #upper CI
    AbsErrorsCI <- CI(AbsErrorsNoiseXNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    AbsErrorsNoiseXNoiseYCI97.5[reduction] <- AbsErrorsCI[1] #upper CI
    AbsErrorsNoiseXNoiseYCI2.5[reduction] <- AbsErrorsCI[3] #upper CI
    
    #calculate confidence intervals for AUC
    CIAUC <- CI(AUCsb[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIAUC97.5[reduction] <- CIAUC[1] #upper CI
    CIAUC2.5[reduction] <- CIAUC[3] #lower CI
    CIAUC <- CI(AUCNoiseX[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIAUCNoiseX97.5[reduction] <- CIAUC[1] #upper CI
    CIAUCNoiseX2.5[reduction] <- CIAUC[3] #lower CI
    CIAUC <- CI(AUCNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIAUCNoiseY97.5[reduction] <- CIAUC[1] #upper CI
    CIAUCNoiseY2.5[reduction] <- CIAUC[3] #lower CI
    CIAUC <- CI(AUCNoiseXNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIAUCNoiseXNoiseY97.5[reduction] <- CIAUC[1] #upper CI
    CIAUCNoiseXNoiseY2.5[reduction] <- CIAUC[3] #lower CI
    
    #calculate confidence intervals for variance
    CIvariance <- CI(variance[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIvariance2.5[reduction] <-  CIvariance[3] #lower band
    CIvariance97.5[reduction] <- CIvariance[1] #upper band
    CIvariance <-  CI(varianceNoiseX[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIvariance2.5NoiseX[reduction] <-  CIvariance[3] #lower band
    CIvariance97.5NoiseX[reduction] <- CIvariance[1] #upper band
    CIvariance <-  CI(varianceNoiseY[reduction,]) #take all reductions steps over the simulations and calculate CI
    CIvariance2.5NoiseY[reduction] <-  CIvariance[3] #lower band
    CIvariance97.5NoiseY[reduction] <- CIvariance[1] #upper band
    CIvariance <-  CI(varianceNoiseXNoiseY[reduction,])  #take all reductions steps over the simulations and calculate CI
    CIvariance2.5NoiseXNoiseY[reduction] <-  CIvariance[3] #lower band
    CIvariance97.5NoiseXNoiseY[reduction] <- CIvariance[1] #upper band
  }
  
  sumAUC <- rowSums(AUCsb)/simulations
  sumAUCNoiseX <- rowSums(AUCNoiseX)/simulations
  sumAUCNoiseY <- rowSums(AUCNoiseY)/simulations
  sumAUCNoiseXNoiseY <- rowSums(AUCNoiseXNoiseY)/simulations
  
  sumAbsErrors <- rowSums(AbsErrors)/simulations
  sumAbsErrorsNoiseX <- rowSums(AbsErrorsNoiseX)/simulations
  sumAbsErrorsNoiseY <- rowSums(AbsErrorsNoiseY)/simulations
  sumAbsErrorsNoiseXNoiseY <- rowSums(AbsErrorsNoiseXNoiseY)/simulations
  
  sumvariance <- rowSums(variance)/simulations
  sumvarianceNoiseX <- rowSums(varianceNoiseX)/simulations
  sumvarianceNoiseY <- rowSums(varianceNoiseY)/simulations
  sumvarianceNoiseXNoiseY <- rowSums(varianceNoiseXNoiseY)/simulations
  #save these values
  #####convergenceIteration(data)#after how many redSteps do analyse plots converge, for plots 
  ####FOR LATER SAVE outcome of convergence####
  
  redStepsconv <- convergenceIteration(data) #after how many reduction steps does boosting iteration converge
  
  save(redStepsconv, AbsErrors, AbsErrorsNoiseX, AbsErrorsNoiseY, AbsErrorsNoiseXNoiseY,
       sumAbsErrors, AbsErrorsCI97.5, AbsErrorsCI2.5, 
       sumAbsErrorsNoiseX, AbsErrorsNoiseXCI97.5, AbsErrorsNoiseXCI2.5,
       sumAbsErrorsNoiseY, AbsErrorsNoiseYCI97.5, AbsErrorsNoiseYCI2.5, 
       sumAbsErrorsNoiseXNoiseY, AbsErrorsNoiseXNoiseYCI97.5, AbsErrorsNoiseXNoiseYCI2.5, 
       AUCsb, AUCNoiseX, AUCNoiseY, AUCNoiseXNoiseY,
       sumAUC, CIAUC97.5, CIAUC2.5,
       sumAUCNoiseX, CIAUCNoiseX97.5, CIAUCNoiseX2.5, 
       sumAUCNoiseY, CIAUCNoiseY97.5, CIAUCNoiseY2.5,
       sumAUCNoiseXNoiseY, CIAUCNoiseXNoiseY97.5, CIAUCNoiseXNoiseY2.5,
       variance, varianceNoiseX, varianceNoiseY, varianceNoiseXNoiseY,
       sumvariance, CIvariance2.5, CIvariance97.5, 
       sumvarianceNoiseX, CIvariance2.5NoiseX, CIvariance97.5NoiseX, 
       sumvarianceNoiseY, CIvariance2.5NoiseY, CIvariance97.5NoiseY,
       sumvarianceNoiseXNoiseY, CIvariance2.5NoiseXNoiseY, CIvariance97.5NoiseXNoiseY,
       file = "EVAL.rda")
  
  nred <- redStepsconv[1,] #convergence of Iteration for AbsError
  nredNoiseX <- redStepsconv[2,]        
  nredNoiseY  <- redStepsconv[3,]       
  nredNoiseXNoiseY  <- redStepsconv[4,] 
  meanred <- sum(nred)/simulations #mean of reduction Steps for Absolute Error over simulations
  meanredNoiseX <- sum(nredNoiseX)/simulations
  meanredNoiseY <- sum(nredNoiseY)/simulations
  meanredNoiseXNoiseY <- sum(nredNoiseXNoiseY)/simulations
  
  nredAUC <- redStepsconv[5,]#convergence of Iteration for AUC
  nredAUCNoiseX <- redStepsconv[6,]
  nredAUCNoiseY <- redStepsconv[7,]
  nredAUCNoiseXNoiseY <- redStepsconv[8,]
  meanredAUC <- sum(nredAUC)/simulations #mean of reduction Steps for AUC over simulations
  meanredAUCNoiseX <- sum(nredAUCNoiseX)/simulations
  meanredAUCNoiseY <- sum(nredAUCNoiseY)/simulations
  meanredAUCNoiseXNoiseY <- sum(nredAUCNoiseXNoiseY)/simulations
  
  #Plots
  pdf("plots.pdf") #save plot
  plot(sumAbsErrors, ylim=range(sumAbsErrors, AbsErrorsCI97.5, AbsErrorsCI2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of absolute Errors", main = "Sum of absolute Errors")
  lines(AbsErrorsCI97.5, col="red") #add CIs to plot
  lines(AbsErrorsCI2.5, col="red")
  #abline(h = sumAbsErrors[meanred], col = "green") #horizontal line that shows convergence
  plot(sumAbsErrorsNoiseX, ylim=range(sumAbsErrorsNoiseX, AbsErrorsNoiseXCI97.5, AbsErrorsNoiseXCI2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of absolute Errors", main = "Sum of absolute Errors for noisy X")
  lines(AbsErrorsNoiseXCI97.5, col="red") #add CIs to plot
  lines(AbsErrorsNoiseXCI2.5, col="red")
  #abline(h = sumAbsErrorsNoiseX[meanredNoiseX], col = "green") #horizontal line that shows convergence
  plot(sumAbsErrorsNoiseY, ylim=range(sumAbsErrorsNoiseY, AbsErrorsNoiseYCI97.5, AbsErrorsNoiseYCI2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of absolute Errors", main = "Sum of absolute Errors for noisy Y")
  lines(AbsErrorsNoiseYCI97.5, col="red") #add CIs to plot
  lines(AbsErrorsNoiseYCI2.5, col="red")
  #abline(h = sumAbsErrorsNoiseY[meanredNoiseY], col = "green") #horizontal line that shows convergence
  plot(sumAbsErrorsNoiseXNoiseY, ylim=range(sumAbsErrorsNoiseXNoiseY, AbsErrorsNoiseXNoiseYCI97.5, AbsErrorsNoiseXNoiseYCI2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of absolute Errors", main = "Sum of absolute Errors for noisy X and noisy Y")
  lines(AbsErrorsNoiseXNoiseYCI97.5, col="red") #add CIs to plot
  lines(AbsErrorsNoiseXNoiseYCI2.5, col="red")
  #abline(h = sumAbsErrorsNoiseXNoiseY[meanredNoiseXNoiseY], col = "green") #horizontal line that shows convergence
  
  # #relative values (to the basic boosting) for first and last simulation
  # plot(AbsErrors[,simulations]/AbsErrors[1, simulations], ylim=range(AbsErrors[,simulations]/AbsErrors[1, simulations], AbsErrors[,1]/AbsErrors[1, 1]), type="l", main = "Relative value of absolute error to basic boosting (first and last simulation)", xlab = "Iteration", ylab = "Relative value")
  # lines(AbsErrors[,1]/AbsErrors[1, 1], type="l", col = "blue")
  # abline(h = 1, col = "red") #if line lays above the line the Abs Error of basic boosting is higher
  # plot(AbsErrorsNoiseX[,simulations]/AbsErrorsNoiseX[1, simulations], ylim=range(AbsErrorsNoiseX[,simulations]/AbsErrorsNoiseX[1, simulations], AbsErrorsNoiseX[,1]/AbsErrorsNoiseX[1, 1]), type="l", main = "Relative value of absolute error to basic boosting for noisy X (first and last simulation)", xlab = "Iteration", ylab = "Relative value")
  # lines(AbsErrorsNoiseX[,1]/AbsErrorsNoiseX[1, 1], type="l", col = "blue")
  # abline(h = 1, col = "red") #if line lays above the line the Abs Error of basic boosting is higher
  # plot(AbsErrorsNoiseY[,simulations]/AbsErrorsNoiseY[1, simulations], ylim=range(AbsErrorsNoiseY[,simulations]/AbsErrorsNoiseY[1, simulations], AbsErrorsNoiseY[,1]/AbsErrorsNoiseY[1, 1]), type="l", main = "Relative value of absolute error to basic boosting for noisy Y (first and last simulation)", xlab = "Iteration", ylab = "Relative value")
  # lines(AbsErrorsNoiseY[,1]/AbsErrorsNoiseY[1, 1], type="l", col = "blue")
  # abline(h = 1, col = "red") #if line lays above the line the Abs Error of basic boosting is higher
  # plot(AbsErrorsNoiseXNoiseY[,simulations]/AbsErrorsNoiseXNoiseY[1, simulations], ylim=range(AbsErrorsNoiseXNoiseY[,simulations]/AbsErrorsNoiseXNoiseY[1, simulations], AbsErrorsNoiseXNoiseY[,1]/AbsErrorsNoiseXNoiseY[1, 1]), type="l", main = "Relative value of absolute error to basic boosting for noisy X and noisy Y (first and last simulation)", xlab = "Iteration", ylab = "Relative value")
  # lines(AbsErrorsNoiseXNoiseY[,1]/AbsErrorsNoiseXNoiseY[1, 1], type="l", col = "blue")
  # abline(h = 1, col = "red") #if line lays above the line the Abs Error of basic boosting is higher
  # #classification AUC
  
  
  #relative values (to the basic boosting)  
  plot(sumAbsErrors/sumAbsErrors[1], main = "Relative value of absolute error to basic boosting ")
  plot(sumAbsErrorsNoiseX/sumAbsErrorsNoiseX[1], main = "Relative value of absolute error to basic boosting ")
  plot(sumAbsErrorsNoiseY/sumAbsErrorsNoiseY[1], main = "Relative value of absolute error to basic boosting ")
  plot(sumAbsErrorsNoiseXNoiseY/sumAbsErrorsNoiseXNoiseY[1], main = "Relative value of absolute error to basic boosting ")
  
  
  plot(sumAUC, ylim=range(sumAUC, CIAUC97.5, CIAUC2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of AUC", main = "Sum aof AUC")
  lines(CIAUC2.5, col="red") #add CIs to plot
  lines(CIAUC97.5, col="red")
  #abline(h = sumAUC[meanredAUC], col = "green") #horizontal line that shows convergence
  plot(sumAUCNoiseX, ylim=range(sumAUCNoiseX, CIAUCNoiseX97.5, CIAUCNoiseX2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of AUC", main = "Sum aof AUC for noisy X")
  lines(CIAUCNoiseX2.5, col="red") #add CIs to plot
  lines(CIAUCNoiseX97.5, col="red")
  #abline(h = sumAUCNoiseX[meanredAUCNoiseX], col = "green") #horizontal line that shows convergence
  plot(sumAUCNoiseY, ylim=range(sumAUCNoiseY, CIAUCNoiseY97.5, CIAUCNoiseY2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of AUC", main = "Sum aof AUC for noisy Y")
  lines(CIAUCNoiseY2.5, col="red") #add CIs to plot
  lines(CIAUCNoiseY97.5, col="red")
  #abline(h = sumAUCNoiseY[meanredAUCNoiseY], col = "green") #horizontal line that shows convergence
  plot(sumAUCNoiseXNoiseY, ylim=range(sumAUCNoiseXNoiseY, CIAUCNoiseXNoiseY97.5, CIAUCNoiseXNoiseY2.5), col='black', type="l", xlab = "Iteration", ylab = "Sum of AUC", main = "Sum aof AUC for noisy X and noisy Y")
  lines(CIAUCNoiseXNoiseY2.5, col="red") #add CIs to plot
  lines(CIAUCNoiseXNoiseY97.5, col="red")
  #abline(h = sumAUCNoiseXNoiseY[meanredAUCNoiseXNoiseY], col = "green") #horizontal line that shows convergence
  #variance
  
  plot(sumvariance, ylim=range(sumvariance, CIvariance2.5, CIvariance97.5), col='black', type = "l", xlab = "Iteration", ylab = "Sum of Variance", main = "Sum of Variance")
  lines(CIvariance97.5, col="red") #add CIs to plot
  lines(CIvariance2.5, col="red")
  plot(sumvarianceNoiseX, ylim=range(sumvarianceNoiseX, CIvariance2.5NoiseX, CIvariance97.5NoiseX), col='black', type = "l", xlab = "Iteration", ylab = "Sum of Variance", main = "Sum of Variance for noisy X")
  lines(CIvariance2.5NoiseX, col="red") #add CIs to plot
  lines(CIvariance97.5NoiseX, col="red")
  plot(sumvarianceNoiseY, ylim=range(sumvarianceNoiseY, CIvariance2.5NoiseY, CIvariance97.5NoiseY), col='black', type = "l", xlab = "Iteration", ylab = "Sum of Variance", main = "Sum of Variance for noisy Y")
  lines(CIvariance2.5NoiseY, col="red") #add CIs to plot
  lines(CIvariance97.5NoiseY, col="red")
  plot(sumvarianceNoiseXNoiseY, ylim=range(sumvarianceNoiseXNoiseY, CIvariance2.5NoiseXNoiseY, CIvariance97.5NoiseXNoiseY), col='black', type = "l", xlab = "Iteration", ylab = "Sum of Variance" , main = "Sum of Variance for noisy Y and noisy X")
  lines(CIvariance2.5NoiseXNoiseY, col="red") #add CIs to plot
  lines(CIvariance97.5NoiseXNoiseY, col="red")
  dev.off() #end of saving plot
  #par(mfrow=c(1,1))
  
  
  
}

########################################################################


#######################################################################

convergenceIteration <- function(dataFile){
  ###FOR LATER -> give back maximum value for no convergence
  ###FOR LATER -> same for noisy X/Y
  library(pracma)
  library(pROC)
  #load(dataFile)
  
  Y <- dataFile$Y                          
  Ynoise <- dataFile$Ynoise
  yBin  <- dataFile$yBin
  yBinNoise <- dataFile$yBinNoise
  originalX <- dataFile$originalX
  originalXnoise <- dataFile$originalXnoise            
  coeffs <- dataFile$coeffs
  predictors <- dataFile$predictors
  kappa <- dataFile$kappa                    
  samples <- dataFile$samples
  simulations <- dataFile$simulations
  nOuterFolds <- dataFile$nOuterFolds               
  redSteps <- dataFile$redSteps
  sampleRatio <- dataFile$sampleRatio
  fixedMstop <- dataFile$fixedMstop                
  fixedNu <- dataFile$fixedNu
  offsetFinal <- dataFile$offsetFinal
  predModelList <- dataFile$predModelList           
  offsetFinalClass <- dataFile$offsetFinalClass
  predModelListClass <- dataFile$predModelListClass
  predictionVector <- dataFile$predictionVector          
  predictionVectorClass <- dataFile$predictionVectorClass
  offsetFinalNoise <- dataFile$offsetFinalNoise
  predModelListNoise <- dataFile$predModelListNoise        
  predictionVectorNoise <- dataFile$predictionVectorNoise
  offsetFinalClassNoise <- dataFile$offsetFinalClassNoise
  predModelListClassNoise <- dataFile$predModelListClassNoise   
  predictionVectorClassNoise <- dataFile$predictionVectorClassNoise
  
  grad <- 0 #variable for gradient of curve
  nred <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredNoiseX <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredNoiseY <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredNoiseXNoiseY <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredAUC <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredAUCNoiseX <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredAUCNoiseY <- rep(0, simulations) #variable to count how many reduction iterations were done
  nredAUCNoiseXNoiseY <- rep(0, simulations) #variable to count how many reduction iterations were done
  
  AbsError <- matrix(0, nrow = redSteps, ncol = simulations) #AbsError to stop reduction iteration
  AbsErrorNoiseX <- matrix(0, nrow = redSteps, ncol = simulations) #AbsError to stop reduction iteration
  AbsErrorNoiseY <- matrix(0, nrow = redSteps, ncol = simulations) #AbsError to stop reduction iteration
  AbsErrorNoiseXNoiseY <- matrix(0, nrow = redSteps, ncol = simulations) #AbsError to stop reduction iteration
  sumAUC <- matrix(0, nrow=redSteps, ncol = simulations) #sumAUC to stop reduction iteration
  sumAUCNoiseX <- matrix(0, nrow=redSteps, ncol = simulations) #sumAUC to stop reduction iteration
  sumAUCNoiseY <- matrix(0, nrow=redSteps, ncol = simulations) #sumAUC to stop reduction iteration
  sumAUCNoiseXNoiseY <- matrix(0, nrow=redSteps, ncol = simulations) #sumAUC to stop reduction iteration
  
  print(simulations)
  for (simulation in 1:simulations) {
    print(simulation)
    stopred <- 0
    #calculate number of iteration when Absolute Error converges
    for (reduction in 1:redSteps){ #for sumabsErrors with no noisy X or Y
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      AbsError[reduction, simulation] <- sum(abs(t(predictionVector[simulation,reduction,])-Y[,simulation]))
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(AbsError[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.07){ 
          nred[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }
    }#end iteration reduction
    
    stopred <- 0
    for (reduction in 1:redSteps){ #for sumabsErrors with no noisy X or Y
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      AbsErrorNoiseX[reduction, simulation] <- sum(abs(t(predictionVectorNoise[simulation,reduction,])-Y[,simulation]))
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(AbsErrorNoiseX[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.07){ 
          nredNoiseX[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }
    }#end iteration reduction
    
    stopred <- 0
    for (reduction in 1:redSteps){ #for sumabsErrors with no noisy X or Y
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      AbsErrorNoiseY[reduction ,simulation] <- sum(abs(t(predictionVector[simulation, reduction,])-Ynoise[,simulation]))
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(AbsErrorNoiseY[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.07){ 
          nredNoiseY[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      } 
    }#end iteration reduction
    
    stopred <- 0
    for (reduction in 1:redSteps){ #for sumabsErrors with no noisy X or Y
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      AbsErrorNoiseXNoiseY[reduction ,simulation] <- sum(abs(t(predictionVectorNoise[simulation,,])-Ynoise[,simulation]))
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(AbsErrorNoiseXNoiseY[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.01){ 
          nredNoiseXNoiseY[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }
    }#end iteration reduction
    
    
    #calculate number of iteration when AUC converges
    stopred <- 0 
    for (reduction in 1:redSteps){
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      sumAUC[reduction, simulation] <- auc(yBin[,simulation], predictionVectorClass[simulation, reduction, ])
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(sumAUC[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.0001){
          nredAUC[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }  
    }#end reduction
    
    stopred <- 0 
    for (reduction in 1:redSteps){
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      sumAUCNoiseX[reduction, simulation] <- auc(yBin[,simulation], predictionVectorClassNoise[simulation, reduction, ])
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(sumAUCNoiseX[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.0001){
          nredAUCNoiseX[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      } 
    }#end reduction
    
    stopred <- 0 
    for (reduction in 1:redSteps){
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      sumAUCNoiseY[reduction, simulation] <- auc(yBinNoise[,simulation], predictionVectorClass[simulation, reduction, ])
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(sumAUCNoiseY[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.0001){
          nredAUCNoiseY[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }
    }#end reduction
    
    stopred <- 0 
    for (reduction in 1:redSteps){
      if(stopred == 1 ){  #if grad was small enough leave reduction iteration
        break
      }
      sumAUCNoiseXNoiseY[reduction, simulation] <- auc(yBinNoise[,simulation], predictionVectorClassNoise[simulation, reduction, ])
      if(as.integer(reduction/10)==(reduction/10)){
        #print(reduction)
        grad <- mean(gradient(sumAUCNoiseXNoiseY[(reduction-9):reduction])) #mean of gradient of error curve
        #print(paste("grad:", grad)) #DEBUG
        if(abs(grad) < 0.0001){
          nredAUCNoiseXNoiseY[simulation] <- reduction #save for checking how many reduction iterations were made
          # print("stop") #DEBUG
          stopred <- 1 #if curve of errors is flat enough don't do the reduction iteration again
        }
      }
    }#end reduction
    
    
  }#end simulation
  
  rbind(nred, nredNoiseX, nredNoiseY, nredNoiseXNoiseY, nredAUC, nredAUCNoiseX, nredAUCNoiseY, nredAUCNoiseXNoiseY) #save values as matrix to use them analyse function
}#end function

setwd("/naslx/projects/ua341/di49suy/sampled-boosting-test-4/mytest3D-files/jobs/01")
load("result.RData")
analyseSimulation(result)

