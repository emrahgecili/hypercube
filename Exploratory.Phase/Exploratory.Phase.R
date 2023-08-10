rm(list=ls(all=TRUE))

#install libraries
library(pROC)
# library(lmenssp)
library(future.apply)
library(MASS)

#set working directory
setwd('~/hypercube/Exploratory.Phase')

#load pre-processed data

load('input_data.RData')

#localgeo_data: 3099 obs. of 68 variables

################################################################################################
#Exploratory.Phase function, signif controls significant level
Exploratory.Phase = function(list.reduction, data, signif=0.01){
  
  library(ggplot2) 
  library(dplyr) 
  
  idx.combinations = t(combn(list.reduction,2))
  
  mat.select.SQ = mat.select.INTER = NULL
  
  #Define function for lmenssp model fitting of checking significant square term
  lmensspRoutine_sq=function(a){
    
    new_solve <- function(A){
      return(chol2inv(chol(A)))}
    
    #lmenssp2
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
    
    #f1: transfer factor variable to numeric
    f1 <- function(x) if (is.factor(x)) {as.numeric(levels(x))[x]} else {x}
    
    d1 = data.frame(data$id,
                    data$fev,
                    data$age,
                    data[,a+3],
                    data$age * data.frame(sapply(data[,a+3], f1)))
    
    colnames(d1) <- c('ID', 'FEV1', 'age', colnames(data[a+3]), 
                      paste("age_",colnames(data[a+3]), sep = ""))
    
    for (i in 1:(2*length(a))){
      d1 <- d1[!is.na(d1[,(i+3)]),]
      if (is.numeric(d1[,(i+3)])){
        d1[,(i+3)] <- (d1[,(i+3)] - mean(d1[,(i+3)]))/sd(d1[,(i+3)])
      }
    }
    
    d1[,paste(colnames(data[a+3]), '_square', sep = "")] <- d1[,4]^2
    
    # for (i in 1:length(a)){
    #     d1[,(i+3)] <- f1(d1[,(i+3)])
    # }
    
    #choose for training data: 80% for training
    idname<-unique(d1$ID)
    subject.choose2<-idname[1:round(length(idname)*0.8)]
    train=d1[d1$ID%in%subject.choose2,]
    #test=d1[!d1$ID%in%subject.choose2,]
    
    ####   second part of the code #####################################
    ###### Model specification ####
    
    #create model formula "temp1"
    fm1<-c('age')
    fm2<-colnames(train[4:ncol(train)])
    fm<-c(fm1,fm2)
    temp1 <- as.formula(paste("FEV1 ~ ", paste(fm, collapse= "+")))
    initial.var<- c(8.273,4.798, 78.837)
    
    #fit the model "lmemssp2", get p-values for all variables
    model1 <- lmenssp2(formula = temp1, data = train,
                       id = train$ID, process = "ibm",init=initial.var,
                       timeVar = train$age, silent = TRUE)
    pVals=model1$estimates[5,4]
    
    return(pVals)
  }
  
  #checking significant square term
  for(ii in 1:length(list.reduction)){
    ### square term, only apply on continuous variable
    if (!is.factor(data[,list.reduction[ii]+3])){
      pVals = lmensspRoutine_sq(list.reduction[ii])
    } else {
      pVals = 1000
    }
    
    if(pVals<signif){
      mat.select.SQ = c(mat.select.SQ,list.reduction[ii])
    }
  }
  
  #Define function for lmenssp model fitting of checking significant interaction term
  lmensspRoutine_inter=function(a){
    
    new_solve <- function(A){
      return(chol2inv(chol(A)))}
    
    #lmenssp2
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
    
    # # f1: transfer factor variable to numeric
    # f1 <- function(x) if (is.factor(x)) {as.numeric(levels(x))[x]} else {x}
    
    d1 = data.frame(data$id,
                    data$fev,
                    data$age,
                    data[,a[1]+3],
                    data[,a[2]+3])
    
    colnames(d1) <- c('ID', 'FEV1', 'age', colnames(data[a[1]+3]), 
                      colnames(data[a[2]+3]))
    
    for (i in 1:length(a)){
      d1 <- d1[!is.na(d1[,(i+3)]),]
      if (is.numeric(d1[,(i+3)])){
        d1[,(i+3)] <- (d1[,(i+3)] - mean(d1[,(i+3)]))/sd(d1[,(i+3)])
      }
    }
    
    d1[,paste(colnames(data[a[1]+3]), '_inter_', colnames(data[a[2]+3]), sep = "")] <- d1[,4] * d1[,5]
    
    # for (i in 1:length(a)){
    #     d1[,(i+3)] <- f1(d1[,(i+3)])
    # }
    
    #choose for training data: 80% for training
    idname<-unique(d1$ID)
    subject.choose2<-idname[1:round(length(idname)*0.8)]
    train=d1[d1$ID%in%subject.choose2,]
    #test=d1[!d1$ID%in%subject.choose2,]
    
    ####   second part of the code #####################################
    ###### Model specification ####
    
    #create model formula "temp1"
    fm1<-c('age')
    fm2<-colnames(train[4:ncol(train)])
    fm<-c(fm1,fm2)
    temp1 <- as.formula(paste("FEV1 ~ ", paste(fm, collapse= "+")))
    initial.var<- c(8.273,4.798, 78.837)
    
    #fit the model "lmemssp2", get p-values for all variables
    model1 <- lmenssp2(formula = temp1, data = train,
                       id = train$ID, process = "ibm",init=initial.var,
                       timeVar = train$age, silent = TRUE)
    pVals=model1$estimates[5,4]
    
    return(pVals)
  }
  
  #checking significant interaction term
  for(ii in 1:nrow(idx.combinations)){
    ### interaction term, only apply on continuous variable
    if ((!is.factor(data[,idx.combinations[ii,1]+3])&(!is.factor(data[,idx.combinations[ii,2]+3])))){
      print(idx.combinations[ii,])
      tryCatch(                       
        expr = {                     
          pVals = lmensspRoutine_inter(idx.combinations[ii,])
        },
        
        error = function(e){
          message(paste("There was an error for: (variable id) "), paste(idx.combinations[ii,1], idx.combinations[ii,2]))
          return(pVals = 1000)
        }
      )
    } else {
      pVals = 1000
    }
    if(pVals<signif){
      mat.select.INTER = rbind(mat.select.INTER,c(idx.combinations[ii,1],idx.combinations[ii,2]))
    }
  }
  
  return(list("mat.select.SQ"=mat.select.SQ,"mat.select.INTER"=mat.select.INTER))  
  
}

################################################################################################

#the function input list.reduction should correspond to the column number of input_data - 3
#the selected variables list (get from Reduction.Phase)

load('Reduction.Phase_out.RData')

ls = Reduction.Phase_out$List.Selection$`Hypercube with dim 2`$numSelected2

Exploratory.Phase_out <- Exploratory.Phase(ls, input_data)

save(Exploratory.Phase_out, file = 'Exploratory.Phase_out.RData')
