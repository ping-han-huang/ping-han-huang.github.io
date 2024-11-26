#R Demo for FDA
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
DTI.baseline <- subset(DTI.complete, visit == 1 & case == 1)[c("ID", "visit", "cca")] #Patients at first visit
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

#=================================FPCA==========================================
DTI_fda <- MakeFPCAInputs(tVec = 1:93, yVec = DTI.baseline$cca)
DTI_fPCA <- FPCA(DTI_fda$Ly, DTI_fda$Lt, list(dataType = "Dense", 
                                                 FVEthreshold = 0.99,
                                                 #methodSelectK='AIC',
                                                 nRegGrid = 93))
round(DTI_fPCA$cumFVE,3)   #cumulative FVE
CreateScreePlot(DTI_fPCA)  #scree plot
CreatePathPlot(DTI_fPCA, subset=sel.crv, main = 'Estimated Paths',
               xlab = 'Tract', ylab = 'FA'); grid()


#=========================Curve Registration====================================
#generate curve registration plots used in prospectus
library(fda)
library(RColorBrewer) #for setting colors in plot

data(growth)
n <- 10

#setting colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


#before registration
age = growth$age
heightbasis = create.bspline.basis(c(1,18),35, 6, age)
heightPar = fdPar(heightbasis, 3, 10^(-0.5))
heightSmooth = smooth.basis(age,growth$hgtf[,1:10], heightPar)
accelUnreg = deriv.fd(heightSmooth$fd, 2)
plot(accelUnreg[,1], lwd = 2, xlab = "Age",
     ylab = "Acceleration", ylim = c(-6,2), 
     lty = rep(2,10),
     col = col_vector)
mean.accelUnreg = mean.fd(accelUnreg)
lines(mean.accelUnreg, lwd = 4, col = "black")


#after registration 
regList = register.fd(yfd = accelUnreg)
accelReg = regList$regfd
plot(accelReg, lwd = 2, xlab = "Age", ylab = "Acceleration",
     ylim = c(-6,2), lty = rep(2,10),
     col = col_vector)
mean.accelReg = mean.fd(accelReg)
lines(mean.accelReg, lwd = 4, col = "black")


#=========================FPCA for Sparse Data==================================
library(refund)

data(cd4)
n <- nrow(cd4)
month <- as.numeric(colnames(cd4)) # months -18 and 42 since seroconversion
m <- ncol(cd4)

#spaghetti plot
matplot(month, t(cd4), type='n', 
        main="CD4 cell counts", ylab="number of cd4 cells", 
        xlab="months since seroconversion" )
for(irow in 1:n){
  temp <- na.omit(data.frame(x = month, y = cd4[irow,]))
  #points(temp$x, temp$y, col="light grey", pch=16, cex=0.6)
  lines(temp$x, temp$y, col="light grey")
}
set.seed(9425)
n.crv <- 5
sel.crv <- sample(1:n, size = n.crv, replace = FALSE)
for(i in 1:n.crv){
  irow <- sel.crv[i]
  temp <- na.omit(data.frame(x = month, y = cd4[irow,]))
  points(temp$x, temp$y, col=rainbow(n)[sel.crv[i]], pch = 16, cex=1)
  lines(temp$x, temp$y, col=rainbow(n)[sel.crv[i]], lwd=2)
}

#sparse FPCA
fpca.res <- fpca.sc(cd4, argvals = month, pve = 0.95, var = TRUE)
m <- length(month)

efns <- fpca.res$efunctions*sqrt(m)
evals <- fpca.res$evalues/m

Yhat <- t(matrix(rep(fpca.res$mu, n), length(month))) + 
  fpca.res$scores %*% t(fpca.res$efunctions)
# Yhat <- fpca.res$Yhat

#plot after sFPCA
n.crv <- 2
sel.crv <- c(167, 170)
matplot(month, t(cd4), type='n', 
        main="CD4 cell counts", ylab="number of cd4 cells", 
        xlab="months since seroconversion" )
for(irow in 1:n){
  temp <- na.omit(data.frame(x = month, y = cd4[irow,]))
  lines(temp$x, temp$y, col="light grey")
}
color = matrix(c("blue", "cornflowerblue", "darkgoldenrod", "darkgoldenrod1"), ncol = 2)
for(i in 1:n.crv){
  irow <- sel.crv[i]
  temp <- na.omit(data.frame(x = month, y = cd4[irow,]))
  points(temp$x, temp$y, col=color[1,i], 
         pch = 16, cex=1)
  lines(temp$x, temp$y, col=color[1,i], lwd=2)
  lines(month, Yhat[sel.crv[i],], 
        col=color[2,i], lwd=2)
}





