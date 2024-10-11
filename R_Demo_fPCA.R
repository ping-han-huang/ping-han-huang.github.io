#R Demo for FDA fPCA
##This R demonstration is part of the presentation held by Ping-Han Huang.
##The code used is based on the following source:
##https://www4.stat.ncsu.edu/~staicu/FDAtutorial/index.html
#==============================================================================


#==============================Load Packages====================================
library(fda)
library(fdapace) 
library(refund)

#================================Load data======================================
data(DTI); attach(DTI)
names(DTI)

DTI.complete <- subset(DTI, complete.cases(DTI)) #removes rows containing NA
DTI.baseline <- subset(DTI.complete, visit == 1 & case == 1)[c("ID", "visit", "cca")] #MS patients at first visit
View(DTI.baseline)
m <- length(unique(DTI.baseline$ID))
dim(DTI.baseline$cca)
tract = 1:93

#==============================Spaghetti Plot===================================
set.seed(111)
matplot(tract, t(DTI.baseline$cca), 
        type='l', lty=1, col="light grey",
        main = "Diffusion Tensor Imaging : CCA",
        xlab="tract", ylab="Fractional anisotropy (FA)")
sel.crv <- sample(1:m, size = 3, replace = FALSE)
matlines(tract, t(DTI.baseline$cca[sel.crv,]), 
         type='l', lty=1, lwd=2, col = rainbow(3))

#=================================fPCA==========================================
DTI_fda <- MakeFPCAInputs(tVec = 1:93, yVec = DTI.baseline$cca)
DTI_fPCA <- FPCA(DTI_fda$Ly, DTI_fda$Lt, list(dataType = "Dense", 
                                                 FVEthreshold = 0.99,
                                                 #methodSelectK='AIC',
                                                 nRegGrid = 93))
round(DTI_fPCA$cumFVE,3)   #cumulative FVE
CreateScreePlot(DTI_fPCA)  #scree plot
CreatePathPlot(DTI_fPCA, subset=sel.crv, main = 'Estimated Paths',
               xlab = 'Tract', ylab = 'FA'); grid()


