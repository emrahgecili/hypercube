rm(list=ls(all=TRUE))
#install libraries
#to install package lmenssp, install from source( .tar.gz): 
#install.packages(path_to_file, repos = NULL, type="source")
library(pROC)
library(MASS)
library(lmenssp)
require(dplyr)
require(tidyr)
library(future.apply)
library(tableHTML)
library(Hmisc)
library(xtable)

options("scipen"=100, "digits"=4)

setwd('~/hypercube/ModelSelection.Phase')

#load pre-processed data
load("input_data.RData")

load('Reduction.Phase_out.RData')
var_list_full <- Reduction.Phase_out$List.Selection$`Hypercube with dim 2`$numSelected2

load('Exploratory.Phase_out.RData')
sq_var_full <- Exploratory.Phase_out[["mat.select.SQ"]]
inter_var_full <- Exploratory.Phase_out[["mat.select.INTER"]]

###########################################################################################################################
#define ModelSelection.Phase function
ModelSelection.Phase = function(data, var_list, sq.terms=NULL, in.terms=NULL, signif=0.01, modelSize=NULL){
  
  lmenssp2 <-
    function(formula, data = NULL, id, process = "bm", timeVar, init = NULL, tol = 1e-5, maxiter = 100, silent = TRUE){ 
      
      #library(nlme)
      #library(MASS)
      #library(geoR)
      
      mf <- model.frame(formula = formula, data = data)
      y  <- as.matrix(model.extract(mf, "response"))
      x  <- as.matrix(model.matrix(attr(mf, "terms"), data = mf))
      colnames(x)[1] <- gsub("[[:punct:]]", "", colnames(x)[1])
      
      nsubj  <- length(unique(id)) # number of subjects
      ntotal <- nrow(y)            # total number of observations
      #idlist <- unique(id)
      
      Time  <- tapply(timeVar, id, function(x) x)
      Intercept <- rep(1, nrow(data))
      data2 <- data.frame(cbind(id, x))
      DM    <- split(data2[, -1], data2$id)
      YM    <- tapply(y, id, function(x) x)
      nobs  <- tapply(id, id, function(x) length(x)) # number of observations per patient
      
      
      ##########################################################
      ################### INTEGRATED BROWNIAN MOTION ###########
      ##########################################################
      
      if(process == "ibm"){
        
        if (length(init) == 0){
          data.init         <- data
          data.init$timeVar <- timeVar
          data.init$id      <- id
          init <- as.numeric(VarCorr(lme(formula, random = ~ timeVar|id, method = "ML", data = data.init))[,1])
        }
        
        theta.new <- init
        theta.old <- init * 20
        tol       <- tol
        iter      <- 1
        Niter     <- maxiter
        
        
        while (sqrt((theta.old - theta.new) %*% (theta.old - theta.new)) > tol & iter <= Niter){
          
          # print(iter)
          theta.old <- theta.new
          
          a <- theta.old[1]  ## omegasq
          b <- theta.old[2]  ## sigmasq
          c <- theta.old[3]  ## nu
          
          
          ### betahat
          
          sum.left.beta <- sum.right.beta <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv  <- new_solve(Vi)
            Xi.transp <- t(Xi)
            
            sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
            sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi
            
          }#for (i in 1 : nsubj)
          
          b.hat <- new_solve(sum.left.beta) %*% sum.right.beta
          
          ####### score for theta
          
          ### a = omegasq
          
          sum.a <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            ri       <- Yi - Xi %*% b.hat 
            
            derVa <- Ji
            
            sum.a <- sum.a + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVa))
            
          }#for (i in 1 : nsubj)
          
          
          ### b = sigmasq
          
          sum.b <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            ri       <- Yi - Xi %*% b.hat 
            
            derVb <- Ri
            
            sum.b <- sum.b + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVb))
            
          }#for (i in 1 : nsubj)
          
          
          # c = tausq
          
          sum.c <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            ri       <- Yi - Xi %*% b.hat 
            
            derVc <- Ii
            
            sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVc))
            
          }#for (i in 1 : nsubj)
          
          score.theta <- 0.5 * matrix(c(sum.a, sum.b, sum.c))
          
          
          ######### INFORMATION MATRIX
          
          #### 1) a = omegasq, 2) a = omegasq
          
          sum.a.a <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            
            derVa <- Ji
            
            sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))
            
          }#for (i in 1 : nsubj)
          
          
          #### 1) a = omegasq, 2) b = sigmasq
          
          sum.a.b <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            
            derVa <- Ji
            derVb <- Ri
            
            sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))
            
          }#for (i in 1 : nsubj)
          
          
          #### 1) a = omegasq, 2) c = tausq
          
          sum.a.c <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            
            derVa <- Ji
            #derVc <- Ii
            #sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))
            
            sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv))
            
          }#for (i in 1 : nsubj)
          
          
          #### 1) b = sigmasq, 2) b = sigmasq
          
          sum.b.b <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            
            derVb <- Ri
            
            sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))
            
          }#for (i in 1 : nsubj)
          
          
          #### 1) b = sigmasq, 2) c = tausq
          
          sum.b.c <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            
            derVb <- Ri
            #derVc <- Ii
            #sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))
            
            sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv))
            
          }#for (i in 1 : nsubj)
          
          
          #### 1) c = tausq, 2) c = tausq
          
          sum.c.c <- 0
          
          for (i in 1 : nsubj){
            
            Timei <- Time[[i]]
            ni    <- nobs[i]
            Xi    <- as.matrix(DM[[i]], nrow = ni)
            Yi    <- as.matrix(YM[[i]], nrow = ni)
            
            Ji <- matrix(1, ni, ni)
            Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
            Ii <- diag(ni)
            
            Vi <- a * Ji + b * Ri + c * Ii
            
            Vi.inv <- new_solve(Vi)
            
            #derVc <- Ii
            #sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))
            sum.c.c  <- sum.c.c + sum(Vi.inv * t(Vi.inv))
            
          }#for (i in 1 : nsubj)
          
          expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c,
                                              sum.a.b, sum.b.b, sum.b.c,
                                              sum.a.c, sum.b.c, sum.c.c), 
                                            ncol = 3, byrow = T)
          
          theta.new <- abs(as.numeric(theta.old - ginv(expected.hessian) %*% score.theta))
          
          if(silent == FALSE){
            
            cat("iteration = ", iter, "\n")
            cat("theta.old = ", theta.old, "\n")
            cat("beta.hat  = ", b.hat, "\n")
            cat("theta.new = ", theta.new, "\n")
            cat("sqrt.diff=", sqrt((theta.old - theta.new) %*% (theta.old - theta.new)), "\n")
            cat("score=", score.theta, "\n")      
            print("-----------------------------")
            
          }
          
          iter <- iter + 1
          
        }#while
        
        
        theta.old <- theta.new
        
        a <- theta.old[1]  ## omegasq
        b <- theta.old[2]  ## sigmasq
        c <- theta.old[3]  ## tausq
        
        ### betahat
        
        sum.left.beta <- sum.right.beta <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          Xi.transp <- t(Xi)
          
          sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
          sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi
          
        }#for (i in 1 : nsubj)
        
        b.hat    <- new_solve(sum.left.beta) %*% sum.right.beta
        b.varcov <- new_solve(sum.left.beta)
        
        ####### score for theta
        
        ### a = omegasq
        
        sum.a <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          ri       <- Yi - Xi %*% b.hat 
          
          derVa <- Ji
          
          sum.a <- sum.a + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVa))
          
        }#for (i in 1 : nsubj)
        
        
        ### b = sigmasq
        
        sum.b <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          ri       <- Yi - Xi %*% b.hat 
          
          derVb <- Ri
          
          sum.b <- sum.b + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVb))
          
        }#for (i in 1 : nsubj)
        
        
        # c = tausq
        
        sum.c <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          ri       <- Yi - Xi %*% b.hat 
          
          derVc <- Ii
          
          sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVc))
          
        }#for (i in 1 : nsubj)
        
        score.theta <- 0.5 * matrix(c(sum.a, sum.b, sum.c))
        
        ######### INFORMATION MATRIX
        
        #### 1) a = omegasq, 2) a = omegasq
        
        sum.a.a <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVa <- Ji
          
          sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) a = omegasq, 2) b = sigmasq
        
        sum.a.b <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVa <- Ji
          derVb <- Ri
          
          sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) a = omegasq, 2) c = tausq
        
        sum.a.c <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVa <- Ji
          #derVc <- Ii
          #sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))
          
          sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) b = sigmasq, 2) b = sigmasq
        
        sum.b.b <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVb <- Ri
          
          sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) b = sigmasq, 2) c = tausq
        
        sum.b.c <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          derVb <- Ri
          #derVc <- Ii
          #sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))
          
          sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv))
          
        }#for (i in 1 : nsubj)
        
        
        #### 1) c = tausq, 2) c = tausq
        
        sum.c.c <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          #derVc <- Ii
          #sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))
          
          sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% Vi.inv))
          
        }#for (i in 1 : nsubj)
        
        expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c,
                                            sum.a.b, sum.b.b, sum.b.c,
                                            sum.a.c, sum.b.c, sum.c.c), 
                                          ncol = 3, byrow = T)
        
        sd.theta <- sqrt(diag(ginv(-expected.hessian)))
        
        ## loglik
        
        sum.loglik <- 0
        
        for (i in 1 : nsubj){
          
          Timei <- Time[[i]]
          ni    <- nobs[i]
          Xi    <- as.matrix(DM[[i]], nrow = ni)
          Yi    <- as.matrix(YM[[i]], nrow = ni)
          
          Ji <- matrix(1, ni, ni)
          Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
          Ii <- diag(ni)
          
          Vi <- a * Ji + b * Ri + c * Ii
          
          Vi.inv <- new_solve(Vi)
          
          ri     <- Yi - Xi %*% b.hat
          
          sum.loglik  <- sum.loglik + as.numeric(determinant(Vi, logarithm=TRUE)$modulus) +  t(ri) %*% Vi.inv %*% ri
          
        }#for (i in 1 : nsubj)
        
        max.loglik <- as.numeric(- 0.5 * ntotal * log(2 * pi) - 0.5 * sum.loglik)
        
        
        result <- cbind(c(c(b.hat), theta.new) , c(sqrt(diag(b.varcov)), sd.theta),
                        c(c(b.hat)/sqrt(diag(b.varcov)), NA, NA, NA))
        result <- cbind(result, 2 * (1 - pnorm(abs(result[, 3]))))
        colnames(result) <- c("Estimate", "Standard error", "Z-estimate", "p-value")
        rownames(result)[(nrow(result) - 2) : nrow(result)] <- c("omegasq", "sigmasq", "tausq")
        
        score.theta <- matrix(score.theta, nrow = 1)
        colnames(score.theta) <- c("omegasq", "sigmasq", "tausq")
        
        rand.varcov <- ginv(-expected.hessian) 
        colnames(rand.varcov) <- rownames(rand.varcov) <- c("omegasq", "sigmasq", "tausq")
        
        output           <- list()
        output$title     <- "Mixed effects model with random intercept and integrated Brownian motion"
        output$date      <- date()
        output$estimates <- result
        output$maxloglik <- max.loglik 
        output$score     <- score.theta
        output$fix.varcov  <- b.varcov 
        output$rand.varcov <- rand.varcov
        output
        
      }#ibm
      
      output
      
    }
  
  new_solve <- function(A){
    return(chol2inv(chol(A)))}
  
  f1 <- function(x) if (is.factor(x)) {as.numeric(levels(x))[x]} else {x}
  
  modelselectionRoutine = function(var_name){
    
    print(var_name)
    
    f_name <- function(input_char, id_char) {
      if (!is.null(id_char)){
        a = grep(id_char, input_char, value = TRUE)
      } else{
        a = input_char
      }
      if (identical(a, character(0))) {a = NULL} else {a}
      return(a)
    }
    
    sq_var_name = f_name(var_name, '_square')
    inter_var_name = f_name(var_name, '_inter_')
    main_var_name = f_name(var_name[!var_name%in%c(sq_var_name,inter_var_name)], NULL)
    age_inter_name = if(!is.null(main_var_name)) {paste("age_",main_var_name, sep = "")} else {NULL}
    
    fm <- c('age', main_var_name, age_inter_name, sq_var_name, inter_var_name)
    
    #formula
    temp <- as.formula(paste("FEV1 ~ ", paste(fm, collapse= "+")))
    
    initial.var<- c(8.273,4.798, 78.837)
    
    model <- lmenssp2(formula = temp, data = test,
                      id = test$ID, process = "ibm",
                      timeVar = test$age, init=initial.var, silent = TRUE)
    
    return(list('logLik' = model$maxloglik, 'df' = dim(model$estimates)[1]))
  }
  
  d1 = data.frame(data$id,
                  data$fev,
                  data$age,
                  data[,var_list+3],
                  data$age * data.frame(sapply(data[,var_list+3], f1)))
  
  colnames(d1) <- c('ID', 'FEV1', 'age', colnames(data[var_list+3]), 
                    paste("age_",colnames(data[var_list+3]), sep = ""))
  
  for (i in 1:(2*length(var_list))){
    d1 <- d1[!is.na(d1[,i+3]),]
    if (is.numeric(d1[,i+3])){
      d1[,i+3] <- (d1[,i+3] - mean(d1[,i+3]))/sd(d1[,i+3])
    }
  }
  
  SQ.names = NULL
  if (!is.null(sq.terms)){
    for (i in sq.terms){
      nm  = colnames(data[i+3])
      d1[,paste(nm, '_square', sep = "")] <- d1[,nm]^2
      
      SQ.names = c(SQ.names,paste(nm, '_square', sep = ""))
    }
  }
  
  inter.names = NULL
  if (!is.null(in.terms)){
    for (i in 1: nrow(in.terms)){
      nm_1  = colnames(data[in.terms[i,1]+3])
      nm_2  = colnames(data[in.terms[i,2]+3])
      d1[,paste(nm_1, '_inter_', nm_2, sep = "")] <- d1[,nm_1] * d1[,nm_2]
      
      inter.names = c(inter.names, paste(nm_1, '_inter_', nm_2, sep = ""))
    }
  }
  
  idname<-unique(d1$ID)
  subject.choose2<-idname[1:round(length(idname)*0.8)]
  # train=d1[d1$ID%in%subject.choose2,]
  test=d1[!d1$ID%in%subject.choose2,]
  
  ###########################Select subset models, perform likelihood ratio test######################
  setSelected = c(colnames(data)[var_list+3], SQ.names, inter.names)
                  
  if(is.null(modelSize)){ modelSize=min(5,length(setSelected)) }
  
  ### Error message
  if(modelSize>7){
    stop('Sorry, this version only support model sizes<8')
  }
  
  goodModels = list()
  
  modelBig = modelselectionRoutine(setSelected)
  LBig = modelBig$logLik
  dfBig = modelBig$df
  
  for (j in 1:modelSize){

    if (!is.null(setSelected)){
      if (length(setSelected) == 1){
        combinationMatrixNames=as.matrix(setSelected)
      } else{
        combinationMatrixNames=cbind(t(combn(setSelected,j)))
      }
      
      logicFitVectorF=matrix(0,nrow(combinationMatrixNames))
      
      for (l in 1:nrow(combinationMatrixNames)){
        XSelect = combinationMatrixNames[l,]
        modelSmall = modelselectionRoutine(XSelect)
        LSmall = modelSmall$logLik
        dfSmall = modelSmall$df
        
        chiCrit = qchisq(1-signif,df = dfBig - dfSmall)
        
        if(is.nan(chiCrit)){
          logicFitVectorF[l]=1
        } else if (-2*(LSmall-LBig)<=chiCrit){
          logicFitVectorF[l]=1
        }
      }
      
      if (j==1) {
        goodModels[[paste('Model Size', j)]] = 
          as.matrix(combinationMatrixNames[which(logicFitVectorF>0),])
      } else{
        goodModels[[paste('Model Size', j)]] = 
          combinationMatrixNames[which(logicFitVectorF>0),]
      }
      
      # help reduce the running time: if a model includes a set of variable is 
      # significant, then a model includes the super set of it should always be significant
      exclude_var = unique(c(combinationMatrixNames[which(logicFitVectorF>0),]))
      setSelected <- setSelected[!setSelected%in%exclude_var]
      setSelected <- if (identical(setSelected, character(0))) {NULL} else {setSelected}
      
    } else{
      goodModels[[paste('Model Size', j)]] = NULL
    }

  }
  
  return(goodModels)
  
}





###########################################################################################################################

ModelSelection.Phase_out <- ModelSelection.Phase(data = input_data, 
                                                 var_list = var_list_full,
                                                 sq.terms = sq_var_full, 
                                                 in.terms = inter_var_full,
                                                 signif = 0.01
                                                 )

save(ModelSelection.Phase_out, file = 'ModelSelection.Phase_out.RData')










