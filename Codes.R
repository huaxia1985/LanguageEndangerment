#data = full data; Wphy, Wsp, Wnb = autocorrelative matrices for the full data
#variable = list of predictors; n.var = number of predictors
#var.group = groups of predictors + interaction terms with regions

library(ordinalNet)

###############
#divide the data into train data and test data
###############

N <- dim(data)[1]
idx <- sample.int(N,size=floor(2*N/3))
data1 <- data[idx,]
data2 <- data[-idx,]
Wphy <- phylomatrix
diag(Wphy) <- 0
tmp <- rowSums(Wphy)
Wphy <- t(Wphy)/rowSums(Wphy)
Wphy[tmp==0,] <- 1/(N-1)
diag(Wphy) <- 0
Wsp <- exp(-spmatrix)
diag(Wsp) <- 0
tmp <- rowSums(Wsp)
Wsp <- t(Wsp)/rowSums(Wsp)
Wsp[tmp==0,] <- 1/(N-1)
diag(Wsp) <- 0
Wnb <- nbmatrix
diag(Wnb) <- 0
tmp <- rowSums(Wnb)
Wnb <- t(Wnb)/rowSums(Wnb)
Wnb[tmp==0,] <- 1/(N-1)
diag(Wnb) <- 0
Wphy1 <- phylomatrix[idx,idx]
diag(Wphy1) <- 0
tmp <- rowSums(Wphy1)
Wphy1 <- t(Wphy1)/rowSums(Wphy1)
Wphy1[tmp==0,] <- 1/(length(idx)-1)
diag(Wphy1) <- 0
Wphy2 <- phylomatrix[-idx,-idx]
diag(Wphy2) <- 0
tmp <- rowSums(Wphy2)
Wphy2 <- t(Wphy2)/rowSums(Wphy2)
Wphy2[tmp==0,] <- 1/(N-length(idx)-1)
diag(Wphy2) <- 0
Wsp1 <- exp(-spmatrix[idx,idx])
diag(Wsp1) <- 0
tmp <- rowSums(Wsp1)
Wsp1 <- t(Wsp1)/rowSums(Wsp1)
Wsp1[tmp==0,] <- 1/(length(idx)-1)
diag(Wsp1) <- 0
Wsp2 <- exp(-spmatrix[-idx,-idx])
diag(Wsp2) <- 0
tmp <- rowSums(Wsp2)
Wsp2 <- t(Wsp2)/rowSums(Wsp2)
Wsp2[tmp==0,] <- 1/(N-length(idx)-1)
diag(Wsp2) <- 0
Wnb1 <- nbmatrix[idx,idx]
diag(Wnb1) <- 0
tmp <- rowSums(Wnb1)
Wnb1 <- t(Wnb1)/rowSums(Wnb1)
Wnb1[tmp==0,] <- 1/(length(idx)-1)
diag(Wnb1) <- 0
Wnb2 <- nbmatrix[-idx,-idx]
diag(Wnb2) <- 0
tmp <- rowSums(Wnb2)
Wnb2 <- t(Wnb2)/rowSums(Wnb2)
Wnb2[tmp==0,] <- 1/(N-length(idx)-1)
diag(Wnb2) <- 0

###############
#test parallel assumption for each independent variable: check p.value
###############

X <- data1[,variable[1:n.var,1]]
X <- as.matrix(X)
y <- data1$EGIDS_tr
P <- dim(X)[2]
N <- dim(X)[1]
p.value <- numeric(P)
names(p.value) <- colnames(X)
coefs <- matrix(NA,2*(length(levels(as.factor(y)))-1)+1,P)
colnames(coefs) <- colnames(X)
for (i in 1:P) {
	tmp <- which(!is.na(X[,i]))
	res1 <- ordinalNet(as.matrix(X[tmp,i],length(tmp),1),y[tmp],standardize=F,family="cumulative",link="probit",nonparallelTerms=F,lambdaVals=0)
	res2 <- ordinalNet(as.matrix(X[tmp,i],lentth(tmp),1),y[tmp],standardize=F,family="cumulative",link="probit",nonparallelTerms=T,lambdaVals=0)
	likratio <- res2$loglik-res1$loglik
	p.value[i] <- 1-pchisq(likratio,df=res2$nNonzero-res1$nNonzero)
	coefs[,i] <- res2$coefs
}

#group 1:6a to 6a
y1 <- y
levels(y1)[1:6] <- "6a"
p.value <- numeric(P)
names(p.value) <- colnames(X)
coefs <- matrix(NA,2*(length(levels(as.factor(y1)))-1)+1,P)
colnames(coefs) <- colnames(X)
for (i in 1:P) {
	tmp <- which(!is.na(X[,i]))
	res1 <- ordinalNet(as.matrix(X[tmp,i],length(tmp),1),y1[tmp],standardize=F,family="cumulative",link="probit",nonparallelTerms=F,lambdaVals=0)
	res2 <- ordinalNet(as.matrix(X[tmp,i],lentth(tmp),1),y1[tmp],standardize=F,family="cumulative",link="probit",nonparallelTerms=T,lambdaVals=0)
	likratio <- res2$loglik-res1$loglik
	p.value[i] <- 1-pchisq(likratio,df=res2$nNonzero-res1$nNonzero)
	coefs[,i] <- res2$coefs
}

###############
# get the ML estimates of the coefficients for the three autocorrelative matrices
###############

#custom function for optim to get the ML estimates of coefficients for autocorrelative matrices in ordinal probit regression. Data points with incomplete data are removed.
#a is the matrices coefficients
#y is the response vector
#X is the predictors' matrix
#Wsp, Wphy, Wnb are the autocorrelative matrices
autoord <- function (a,y,X,Wsp,Wphy,Wnb) {
		tmp <- which(!is.na(rowSums(X)))
        X.tmp <- X[tmp,]
        y.tmp <- y[tmp]
		W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
        W <- W[tmp,tmp]
        Wy <- W%*%as.numeric(y.tmp)
        X2 <- W%*%X.tmp
        X2 <- cbind(X.tmp,X2)
		res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X.tmp,res)
        out <- ordinalNet(X1,y.tmp,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,lambdaVals=0)
        -out$loglik
}

intercept <- variable[variable[,2]=="0",1]
X <- data1[,intercept]
X <- as.matrix(X)
res <- optim(c(0.5,0.5),autoord,method="L-BFGS-B",lower=c(0,0),upper=c(1,1),y=y1,X=X,Wsp=Wsp1,Wphy=Wphy1,Wnb=Wnb1)
a <- res$par

###############
#stepwise selection procedure for the best-fit models
###############

#custom function to fit an ordinal probit regression model given the autocorrelative structure that is inferred above.
#y is the response vector
#X is the predictors' matrix, which defines the regression model
#W is the autocorrelative structure
autoord <- function (y,X,W) {
        Wy <- W%*%as.numeric(y)
        X2 <- W%*%X
        X2 <- cbind(X,X2)
        res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X,res)
        ordinalNet(X1,y,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,lambdaVals=0)
}

X <- data1[,variable[1:n.var,1]]
X <- as.matrix(X)
loglik <- numeric(P)
names(loglik) <- colnames(X)
tmp <- which(!is.na(rowSums(X)))
X.tmp <- X[tmp,]
y.tmp <- y1[tmp]
W <- a[2]*(a[1]*Wsp1+(1-a[1])*Wnb1)+(1-a[2])*Wphy1
W <- W[tmp,tmp]
for (i in 1:P) {
    loglik[i] <- autoord(y=y.tmp,X=as.matrix(cbind(X.tmp[,i],data1[tmp,intercept])),W=W)$loglik
}
#the best fit single predictor is lang_L1.POP_lang_tr

library(foreach)
library(doParallel)
registerDoParallel(8)
X <- data1[,variable[,1]]
X <- as.matrix(X)
tmp <- which(!is.na(rowSums(X)))
X.tmp <- X[tmp,]
y.tmp <- y1[tmp]
W <- a[2]*(a[1]*Wsp1+(1-a[1])*Wnb1)+(1-a[2])*Wphy1
W <- W[tmp,tmp]
aa <- 1  #repeat aa=1:50
set.seed(aa)
basemodel <- c("lang_L1.POP_lang_tr",intercept)
baseloglik <- min(loglik)
var.group1 <- var.group
var.group1[[2]] <- var.group1[[2]][-1]
while (!setequal(basemodel,basemodel1)) {
    basemodel1 <- basemodel
    for (j in sample(1:length(var.group),length(var.group),replace=F)) {
    	mark <- 1
        vars=var.group1[[j]]
        if (length(vars)>0) {
    		res <- vector("list",length(vars))
    		names(res) <- vars
    		logLik <- numeric(length(vars))
    		names(logLik) <- vars
    		res <- foreach (i = 1:length(vars), .packages='ordinalNet') %dopar% {
       			 altermodel <- c(basemodel,vars[i])
       			 autoord(y=y.tmp,X=X.tmp[,altermodel],W=W)
    		}
    		logLik <- sapply(res,function (z) z$loglik)
    		names(logLik) <- vars
    		names(res) <- vars
    		likratio.test <- 1-pchisq(2*(logLik-baseloglik),df=1)
    		tmp1 <- vars[likratio.test==min(likratio.test)]
    		if (length(tmp1)>1) {
        		tmp1 <- tmp1[which(logLik[tmp1]==max(logLik[tmp1]))[1]]
    		}
    		if (min(likratio.test)<0.05) {
        		basemodel <- c(basemodel,tmp1)
        		baseloglik <- logLik[tmp1]
        		res <- list(res[[tmp1]])
        		var <- tmp1
        		var.group1[[j]] <- vars[-which(vars==var)]
        	}
        }
        mark <- 2
    	vars <- basemodel[!is.element(basemodel,intercept)]
    	if (length(vars)>0) {
    		res <- vector("list",length(vars))
    		names(res) <- vars
    		logLik <- numeric(length(vars))
    		names(logLik) <- vars
    		names(res) <- vars
    		res <- foreach (i = 1:length(vars), .packages='ordinalNet') %dopar% {
        		altermodel <- basemodel[-which(basemodel==vars[i])]
        		autoord(y=y.tmp,X=X.tmp[,altermodel],W=W)
    		}
    		logLik <- sapply(res,function (z) z$loglik)
    		names(logLik) <- vars
    		names(res) <- vars
    		likratio.test <- 1-pchisq(2*(-logLik+baseloglik),df=1)
    		tmp1 <- vars[likratio.test==max(likratio.test)]
    		if (length(tmp1)>1) {
       			 tmp1 <- tmp1[which(logLik[tmp1]==max(logLik[tmp1]))[1]]
    		}
    		if (max(likratio.test)>0.05) {
        		basemodel <- basemodel[-which(basemodel==tmp1)]
        		baseloglik <- logLik[tmp1]
        		res <- list(res[[tmp1]])
        		var <- tmp1
        		var.group1[[j]] <- c(var.group1[[j]],var)
    		}
    	}
    }
}
filename <- paste0("basemodel",aa,".csv")
write.csv(basemodel,filename)

###############
#apply cross-validation and lasso to each selected model
###############

#custom function to fit an ordinal probit regression model with lasso given the autocorrelative structure.
#y is the response vector
#X is the predictors' matrix, which defines the regression model
#W is the autocorrelative structure
autoord <- function (y,X,W) {
        Wy <- W%*%as.numeric(y)
        X2 <- W%*%X
        X2 <- cbind(X,X2)
        res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X,res)
        ordinalNetCV(X1,y,nFolds=10,alpha=0.5,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,lambdaVals=0,maxiterOut=1000)
}
#custom function to fit an ordinal probit regression model with cross-validation given the autocorrelative structure.
#y is the response vector
#X is the predictors' matrix, which defines the regression model
#W is the autocorrelative structure
autoord.lasso <- function (y,X,W) {
        Wy <- W%*%as.numeric(y)
        X2 <- W%*%X
        X2 <- cbind(X,X2)
        res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X,res)
        ordinalNetCV(X1,y,nFolds=10,alpha=0.5,standardize=F,family="cumulative",link="probit",nonparallelTerms=F,tuneMethod="cvLoglik",maxiterOut=1000)
}

res.ind1 <- vector("list",50)
res.ind2 <- vector("list",50)
for (a in 1:50) {
    model <- read.csv(paste0("basemodel",i,".csv"),header=T,row.names=1)
    res.ind1[[i]] <- autoord.lasso(y=y.tmp,X=X.tmp[,model[[1]]],W=W)
    res.ind2[[i]] <- autoord(y=y.tmp,X=X.tmp[,model[[1]]],W=W)
}

i <- 1
bes <- 20
tmp <- -res.ind1[[i]]$fit$coefs[bes,-c(1:7,length(res.ind1[[i]]$fit$coefs[bes,]))]
coe <- numeric(length(variable[,1]))
names(coe) <- variable[,1]
coe[names(tmp)] <- tmp
coefs.ind <- c(coe,res.ind1[[i]]$fit$loglik[bes],res.ind1[[i]]$fit$aic[bes],res.ind1[[i]]$fit$bic[bes],res.ind1[[i]]$fit$devPct[bes],mean(res.ind1[[i]]$brier),sd(res.ind1[[i]]$brier),mean(res.ind1[[i]]$devPct),sd(res.ind1[[i]]$devPct),mean(res.ind1[[i]]$misclass),sd(res.ind1[[i]]$misclass))
for (i in 2:50) {
	bes <- 20
	tmp <- -res.ind1[[i]]$fit$coefs[bes,-c(1:7,length(res.ind1[[i]]$fit$coefs[bes,]))]
	coe <- numeric(length(variable[,1]))
	names(coe) <- variable[,1]
	coe[names(tmp)] <- tmp
	coefs.ind <- rbind(coefs.ind,c(coe,res.ind1[[i]]$fit$loglik[bes],res.ind1[[i]]$fit$aic[bes],res.ind1[[i]]$fit$bic[bes],res.ind1[[i]]$fit$devPct[bes],mean(res.ind1[[i]]$brier),sd(res.ind1[[i]]$brier),mean(res.ind1[[i]]$devPct),sd(res.ind1[[i]]$devPct),mean(res.ind1[[i]]$misclass),sd(res.ind1[[i]]$misclass)))
}
write.csv(coefs.ind,file="coefs.ind.lasso.csv")

i <- 1
bes <- 1
tmp <- -res.ind2[[i]]$fit$coefs[bes,-c(1:7,length(res.ind2[[i]]$fit$coefs[bes,]))]
coe <- numeric(length(variable[,1]))
names(coe) <- variable[,1]
coe[names(tmp)] <- tmp
coefs.ind <- c(coe,res.ind2[[i]]$fit$loglik[bes],res.ind2[[i]]$fit$aic[bes],res.ind2[[i]]$fit$bic[bes],res.ind2[[i]]$fit$devPct[bes],mean(res.ind2[[i]]$brier),sd(res.ind2[[i]]$brier),mean(res.ind2[[i]]$devPct[is.finite(res.ind2[[i]]$devPct)]),sd(res.ind2[[i]]$devPct[is.finite(res.ind2[[i]]$devPct)]),mean(res.ind2[[i]]$misclass),sd(res.ind2[[i]]$misclass))
for (i in 2:50) {
	bes <- which(res.ind2[[i]]$fit$bic==min(res.ind2[[i]]$fit$bic))
	tmp <- -res.ind2[[i]]$fit$coefs[bes,-c(1:7,length(res.ind2[[i]]$fit$coefs[bes,]))]
	coe <- numeric(length(variable[,1]))
	names(coe) <- variable[,1]
	coe[names(tmp)] <- tmp
	coefs.ind <- rbind(coefs.ind,c(coe,res.ind2[[i]]$fit$loglik[bes],res.ind2[[i]]$fit$aic[bes],res.ind2[[i]]$fit$bic[bes],res.ind2[[i]]$fit$devPct[bes],mean(res.ind2[[i]]$brier),sd(res.ind2[[i]]$brier),mean(res.ind2[[i]]$devPct[is.finite(res.ind2[[i]]$devPct)]),sd(res.ind2[[i]]$devPct[is.finite(res.ind2[[i]]$devPct)]),mean(res.ind2[[i]]$misclass),sd(res.ind2[[i]]$misclass)))
}
write.csv(coefs.ind,file="coefs.ind.csv")

###############
#calculate predictive power of each selected model on data2
###############

yFactorToMatrix <- function(y)
{
    nObs <- length(y)
    nLev <- length(levels(y))
    yMat <- matrix(0, nrow=nObs, ncol=nLev, dimnames=list(NULL, levels(y)))
    yInt <- as.integer(y)
    yMat[cbind(1:nObs, yInt)] <- 1
    yMat
}
getLoglik <- function(pMat, yMat)
{
    pkplusone <- 1 - rowSums(pMat)
    pkplusone[pkplusone<0] <- 0
    pMatFull <- cbind(pMat, pkplusone)
    if (any(pMatFull < 0)) return(-Inf)
    llMat <- yMat * log(pMatFull)
    llMat[yMat==0] <- 0 
    llik <- sum(llMat)
    llik
}

getMisclass <- function(pMat, yMat)
{
    pkplusone <- 1 - rowSums(pMat)
    pMatFull <- cbind(pMat, pkplusone)
    predClass <- apply(pMatFull, 1, which.max)
    nMisclass <- sapply(1:nrow(yMat), function(i) sum(yMat[i, -predClass[i]]))
    misclass <- sum(nMisclass) / sum(yMat)
    misclass
}

getBrier <- function(pMat, yMat)
{
    pkplusone <- 1 - rowSums(pMat)
    pMatFull <- cbind(pMat, pkplusone)
    n <- rowSums(yMat)
    brier <- sum(yMat * (1 - pMatFull)^2 + (n - yMat) * pMatFull^2) / sum(n)
    brier
}
getLoglikNull <- function(yMat)
{
    pHatNull <- colSums(yMat) / sum(yMat)
    llmat0 <- yMat * rep(log(pHatNull), each=nrow(yMat))
    llmat0[yMat == 0] <- 0 
    loglik0 <- sum(llmat0)
    loglik0
}

#custom function to use a fitted ordinal probit regression model given the autocorrelative structure to new data.
#y is the response vector in the new data
#X is the predictors' matrix in the new data
#W is the autocorrelative structure
#object is a fitted ordinal probit regression model
#whichLambda indicates which coefficient for lasso to use
autoord.pred <- function (y,X,W,object,whichLambda) {
        Wy <- W%*%as.numeric(y)
        X2 <- W%*%X
        X2 <- cbind(X,X2)
        res <- lm(Wy~X2)$residuals
        X1 <- cbind(Wy,X,res)
        yMat <- yFactorToMatrix(y)
        y_mat_full <- predict(object=object,newx=X1,whichLambda=whichLambda,type="response")
        y_mat <- y_mat_full[,-ncol(y_mat_full),drop=FALSE]
        y_mat[y_mat==0] <- min(y_mat[y_mat!=0])/10
        y_hat <- predict(object=object,newx=X1,whichLambda=whichLambda,type="class")
        loglik <- getLoglik(y_mat, yMat)
        misclass <- getMisclass(y_mat, yMat)
        brier <- getBrier(y_mat, yMat)
        loglikNull <- getLoglikNull(yMat)
        devPct <- 1 - loglik / loglikNull
        list(loglik=loglik,misclass=misclass,brier=brier,devPct=devPct,y_mat=y_mat_full,y_hat=y_hat)
}

X <- data2[,variable[,1]]
X <- as.matrix(X)
y <- data2$EGIDS_tr
y2 <- y
levels(y2)[1:6] <- "6a"
tmp <- which(!is.na(rowSums(X)))
y.tmp <- y2[tmp]
X.tmp <- X[tmp,]
W <- a[2]*(a[1]*Wsp2+(1-a[1])*Wnb2)+(1-a[2])*Wphy2
W <- W[tmp,tmp]

predict.ind <- vector("list",50)
predict.ind.lasso <- vector("list",50)
for (i in 1:50) {
	predict.ind[[i]] <- autoord.pred(y=y.tmp,X=X.tmp[,colnames(res.ind2[[i]]$fit$coefs)[-c(1:7,length(res.ind2[[i]]$fit$coefs))]],W=W,object=res.ind2[[i]]$fit,whichLambda=1)
	predict.ind.lasso[[i]] <- autoord.pred(y=y.tmp,X=X.tmp[,colnames(res.ind1[[i]]$fit$coefs)[-c(1:7,dim(res.ind1[[i]]$fit$coefs)[2])]],W=W,object=res.ind1[[i]]$fit,whichLambda=20)
}
write.csv(cbind(sapply(1:50,function (i) predict.ind[[i]]$loglik[is.finite(predict.ind[[i]]$loglik)]),sapply(1:50,function (i) predict.ind[[i]]$brier),sapply(1:50,function (i) predict.ind[[i]]$devPct[is.finite(predict.ind[[i]]$devPct)]),sapply(1:50,function (i) predict.ind[[i]]$misclass)),file="predict.ind.cv.csv")
write.csv(cbind(sapply(1:50,function (i) predict.ind.lasso[[i]]$loglik),sapply(1:50,function (i) predict.ind.lasso[[i]]$brier),sapply(1:50,function (i) predict.ind.lasso[[i]]$devPct),sapply(1:50,function (i) predict.ind.lasso[[i]]$misclass)),file="predict.ind.lasso.csv")

###############
#fit the best model to full data
###############

ind <- read.csv("coefs.ind.csv")
lasso <- read.csv("coefs.ind.lasso.csv")
model <- rbind(lasso,ind[-rm,])
model <- colnames(model[,c(2:615)])[which(colSums(model[,c(2:615)]!=0)>28)]
model <- gsub("South.Eastern","South-Eastern",model)
X <- data[,model]
X <- as.matrix(X)
y <- data$EGIDS_tr
y2 <- y
levels(y2)[1:6] <- "6a"
W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
tmp <- which(!is.na(rowSums(X)))
y.tmp <- y2[tmp]
X.tmp <- X[tmp,]
W <- W[tmp,tmp]
best.model <- autoord.lasso(y=y.tmp,X=X.tmp,W=W)
predict.all <- autoord.pred(y=y.tmp,X=X.tmp,W=W,object=best.model$fit,whichLambda=20)

###############
#use the best model to predict 40 years
###############

model_40 <- model
model_40[grep("lang_L1.POP_lang",model)] <- gsub("_tr","_40_tr",model_40[grep("lang_L1.POP_lang",model)])
model_40[grep("pop.density_polygon",model)] <- gsub("_tr","_40_tr",model_40[grep("pop.density_polygon",model)])
model_40[grep("landuse_hfp_polygon",model)] <- gsub("_tr","_40_tr",model_40[grep("landuse_hfp_polygon",model)])
model_40[grep("landuse_croplands_polygon",model)] <- gsub("_tr","_40_tr",model_40[grep("landuse_croplands_polygon",model)])
model_40[grep("landuse_built_polygon",model)] <- gsub("_tr","_40_tr",model_40[grep("landuse_built_polygon",model)])
model_40[grep("landuse_pasture_polygon",model)] <- gsub("_tr","_40_tr",model_40[grep("landuse_pasture_polygon",model)])
model_40[grep("enviro_temperature.seasonality_polygon",model)] <- gsub("_tr","_40_tr",model_40[grep("enviro_temperature.seasonality_polygon",model)])
X <- data[,model_40]
X <- as.matrix(X)
y <- data$EGIDS_tr
y2 <- y
levels(y2)[1:6] <- "6a"
W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
tmp <- which(!is.na(rowSums(X)))
y.tmp <- y2[tmp]
X.tmp <- X[tmp,]
W <- W[tmp,tmp]
predict.40 <- autoord.pred(y=y.tmp,X=X.tmp,W=W,object=best.model$fit,whichLambda=20)
dat <- cbind(predict.40$y_mat,y_hat=predict.40$y_hat,y_exp=predict.40$y_mat%*%c(1:7),y_obs=y.tmp)
colnames(dat)[9] <- "y_avg"

data$AES_tr <- factor(data$response_AES_lang)
data$AES_tr <- factor(data$AES_tr,levels(data$AES_tr)[c(4,6,5,2,3,1)])
data$EGIDS_tr_40 <- NA
data$EGIDS_tr_40[tmp] <- dat[,8]
data$EGIDS_tr_40[which(data$EGIDS_tr=="9")] <- 7
data$EGIDS_tr_40[which(data$EGIDS_tr=="10")] <- 7
data$EGIDS_tr_40[which(data$EGIDS_tr=="8a")] <- 7
data$EGIDS_tr_40[which(data$EGIDS_tr=="8b")] <- 7
data$EGIDS_tr_40[which((data$EGIDS_tr=="7")*(as.numeric(data$EGIDS_tr_40<4))==1)] <- 4
data$EGIDS_tr_40[which(((data$EGIDS_tr=="6b")*(data$lang_L1.POP_lang<1000)*(as.numeric(data$AES_tr)>2)*(data$EGIDS_tr_40<3))==1)] <- 3
dat[,8] <- data$EGIDS_tr_40[tmp]
dat <- as.data.frame(dat)
rownames(dat) <- data[tmp,1]
write.csv(dat,file="project40.csv")

###############
#use the best model to predict 80 years
###############

model_80 <- model
model_80[grep("lang_L1.POP_lang",model)] <- gsub("_tr","_80_tr",model_80[grep("lang_L1.POP_lang",model)])
model_80[grep("pop.density_polygon",model)] <- gsub("_tr","_80_tr",model_80[grep("pop.density_polygon",model)])
model_80[grep("landuse_hfp_polygon",model)] <- gsub("_tr","_80_tr",model_80[grep("landuse_hfp_polygon",model)])
model_80[grep("landuse_croplands_polygon",model)] <- gsub("_tr","_80_tr",model_80[grep("landuse_croplands_polygon",model)])
model_80[grep("landuse_built_polygon",model)] <- gsub("_tr","_80_tr",model_80[grep("landuse_built_polygon",model)])
model_80[grep("landuse_pasture_polygon",model)] <- gsub("_tr","_80_tr",model_80[grep("landuse_pasture_polygon",model)])
model_40[grep("enviro_temperature.seasonality_polygon",model)] <- gsub("_tr","_80_tr",model_40[grep("enviro_temperature.seasonality_polygon",model)])
X <- data[,model_80]
X <- as.matrix(X)
y <- data$EGIDS_tr
y2 <- y
levels(y2)[1:6] <- "6a"
W <- a[2]*(a[1]*Wsp+(1-a[1])*Wnb)+(1-a[2])*Wphy
tmp <- which(!is.na(rowSums(X)))
y.tmp <- y2[tmp]
X.tmp <- X[tmp,]
W <- W[tmp,tmp]
predict.80 <- autoord.pred(y=y.tmp,X=X.tmp,W=W,object=best.model$fit,whichLambda=20)
dat <- cbind(predict.80$y_mat,y_hat=predict.80$y_hat,y_exp=predict.80$y_mat%*%c(1:7),y_obs=y.tmp)
colnames(dat)[9] <- "y_avg"

data$EGIDS_tr_80 <- NA
data$EGIDS_tr_80[tmp] <- dat[,8]
data$EGIDS_tr_80[which(data$EGIDS_tr_40>=4)] <- 7
data$EGIDS_tr_80[which((data$EGIDS_tr_40==3)*(data$EGIDS_tr_80<4)==1)] <- 4
data$EGIDS_tr_80[which(((data$EGIDS_tr_40==2)*(data$lang_L1.POP_lang_40<1000)*(as.numeric(data$AES_tr)>2)*(data$EGIDS_tr_80<3))==1)] <- 3
dat[,8] <- data$EGIDS_tr_80[tmp]
dat <- as.data.frame(dat)
rownames(dat) <- data[tmp,1]
write.csv(dat,file="project80.csv")
