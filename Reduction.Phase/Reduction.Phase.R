rm(list=ls(all=TRUE))

#install libraries
library(pROC)
library(lmenssp)
library(future.apply)
library(MASS)

setwd('~/hypercube/Reduction.Phase')

load("input_data.RData")

n_pred = 50

###################################################################################################
#vector.signif:(2,0.01)
Reduction.Phase = function(data, dmHC=3, vector.signif=NULL, seed.HC=NULL, n=152, d=n_pred){
  
  new_solve <- function(A){
    return(chol2inv(chol(A)))}
  
  lmensspRoutine=function(a,significance=NULL){
    
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
    pVals=model1$estimates[3:(3+length(a)-1),4]
    
    #decide which procedure to use for variable selection
    if(significance==2){
      idxSelected=which(pVals %in% sort(pVals)[1:2])
    } else if(significance==1){
      idxSelected=which(pVals %in% sort(pVals)[1])
    }  else{
      idxSelected=which(pVals<significance)
    }
    
    return(list("idxSelected"=idxSelected))
  }
  
  #vector.signif=NULL: (2,0.01,0.01)
  if(is.null(vector.signif)){
    signif.Default =TRUE
  } else{
    signif.Default =FALSE
    #rev:reverse the order of the elements
    vector.signif = rev(vector.signif)
    vector.signif = c(NA,vector.signif)
  }
  
  highest.dmHC = dmHC
  
  ### Outputs
  Matrix.Selection = list()
  List.Selection = list()
  
  aux.dmHC4 = aux.dmHC3 = aux.dmHC2 = 'N'
  
  ########## case in which dmHC=4 ##########
  if(dmHC==4 ){
    
    if(signif.Default==TRUE & highest.dmHC==4){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>4){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[4]
    }
    
    if(dmHC==4){
      dimHC = ceiling(d^(1/4))
      nearestHypercube=dimHC^4
      remainderHC=nearestHypercube-d
      if(!is.null(seed.HC)){
        set.seed(seed.HC)
        hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC,dimHC))
      } else{
        hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC,dimHC))
      }}
    
    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC,dimHC,dimHC))
    
    if(all(dim(hypercube)==0)){
      stop(paste('No variables selected at stage 4, increase p-value.'))
    }
    
    for(ind4 in 1:dim(hypercube)[4]){
      for(ind3 in 1:dim(hypercube)[3]){
        for(indR in 1:dim(hypercube)[1]){
          if(length(which(hypercube[indR,,ind3,ind4]!=0))>0){
            index = which(hypercube[indR,,ind3,ind4]!=0)
            idx.aux = hypercube[indR,index,ind3,ind4]
            idxSelected <- tryCatch(                       
              expr = {                     
                idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
              },
              
              error = function(e){
                message(paste("There was an error for: (variable id) "),paste(hypercube[indR,,ind3,ind4],collapse = ","))
                return(idxSelected = list('idxSelected' = NULL))
              }
            )
            
            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,index[idxSelected$idxSelected],ind3,ind4]=hypercubeSelect[indR,index[idxSelected$idxSelected],ind3,ind4]+1
            }
          }
        } # indR
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[,indC,ind3,ind4]!=0))>0){
            index = which(hypercube[,indC,ind3,ind4]!=0)
            idx.aux = hypercube[index,indC,ind3,ind4]
            idxSelected <- tryCatch(                       
              expr = {                     
                idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
              },
              
              error = function(e){
                message(paste("There was an error for: (variable id) "),paste(hypercube[,indC,ind3,ind4],collapse = ","))
                return(idxSelected = list('idxSelected' = NULL))
              }
            )
            
            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[index[idxSelected$idxSelected],indC,ind3,ind4]=hypercubeSelect[index[idxSelected$idxSelected],indC,ind3,ind4]+1
            }
          }
        } # indC
      } # ind3
      for(indR in 1:dim(hypercube)[1]){
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[indR,indC,,ind4]!=0))>0){
            index = which(hypercube[indR,indC,,ind4]!=0)
            idx.aux = hypercube[indR,indC,index,ind4]
            idxSelected <- tryCatch(                       
              expr = {                     
                idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
              },
              
              error = function(e){
                message(paste("There was an error for: (variable id) "),paste(hypercube[indR,indC,,ind4],collapse = ","))
                return(idxSelected = list('idxSelected' = NULL))
              }
            )
            
            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,indC,index[idxSelected$idxSelected],ind4]=hypercubeSelect[indR,indC,index[idxSelected$idxSelected],ind4]+1
            }
          }
        } # indC
      } # indR
    } # ind4
    # now traverse in the 4th dimension
    for(ind3 in 1:dim(hypercube)[3]){
      for(indR in 1:dim(hypercube)[1]){
        for(indC in 1:dim(hypercube)[2]){
          if(length(which(hypercube[indR,indC,ind3,]!=0))>0){
            index = which(hypercube[indR,indC,ind3,]!=0)
            idx.aux = hypercube[indR,indC,ind3,index]
            idxSelected <- tryCatch(                       
              expr = {                     
                idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
              },
              
              error = function(e){
                message(paste("There was an error for: (variable id) "),paste(hypercube[indR,indC,ind3,],collapse = ","))
                return(idxSelected = list('idxSelected' = NULL))
              }
            )
            
            if(length(idxSelected$idxSelected)>0){
              hypercubeSelect[indR,indC,ind3,index[idxSelected$idxSelected]]=hypercubeSelect[indR,indC,ind3,index[idxSelected$idxSelected]]+1
            }
          }
        } # indC
      } # indR
    } # ind3
    
    setSelected4Times = hypercube[which(hypercubeSelect>3, arr.ind = TRUE)]
    setSelected3Times = hypercube[which(hypercubeSelect>2, arr.ind = TRUE)]
    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]
    
    numSelected4times5=length(setSelected4Times)
    numSelected3times5=length(setSelected3Times)
    numSelected2times5=length(setSelected2Times)
    numSelected1times5=length(setSelected1Times)
    
    Matrix.Selection[[paste('Hypercube with dim',4)]] = c(numSelected1times5,numSelected2times5,numSelected3times5,numSelected4times5)
    names(Matrix.Selection[[paste('Hypercube with dim',4)]]) = c('numSelected1','numSelected2','numSelected3','numSelected4')
    
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected1')]] = setSelected1Times
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected2')]] = setSelected2Times
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected3')]] = setSelected3Times
    List.Selection[[paste('Hypercube with dim',4)]][[paste('numSelected4')]] = setSelected4Times
    
    aux.dmHC4 <- 'Y' 
    #readline(cat("Reduction of dimension 4 done!", "\n", length(setSelected4Times5),
    #        "Variables selected at least 4 times","\n",length(setSelected3Times5),
    #        "Variables selected at least 3 times","\n",length(setSelected2Times5),
    #       "Variables selected at least 2 times","\n", "Wanna proceed with reduction?[Y/N]"))
  }
  
  ########## case in which dmHC=3 ##########
  if(dmHC==3 | aux.dmHC4=='Y'){
    
    if(signif.Default==TRUE & highest.dmHC==3){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>3){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[3]
    }
    
    if(dmHC==3){
      dimHC = ceiling(d^(1/3))
      nearestHypercube=dimHC^3
      remainderHC=nearestHypercube-d
      if(!is.null(seed.HC)){
        set.seed(seed.HC)
        #arrange variable id in hypercube in a random arrange manner
        hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC))
      } else{
        hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC))
      }} else if(dmHC>3){
        dimHC = ceiling(length(setSelected2Times)^(1/3))
        nearestHypercube=dimHC^3
        remainderHC=nearestHypercube-length(setSelected2Times)
        if(!is.null(seed.HC)){
          set.seed(seed.HC)
          hypercube=array(sample(c(setSelected2Times,rep(0,remainderHC))),dim=c(dimHC,dimHC,dimHC))
        } else{
          hypercube=array(c(setSelected2Times,rep(0,remainderHC)),dim=c(dimHC,dimHC,dimHC))
        }
      }
    
    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC,dimHC))
    
    if(all(dim(hypercube)==0)){
      stop(paste('No variables selected at stage 3, increase p-value.'))
    }
    
    for(indL in 1:dim(hypercube)[3]){
      for(indR in 1:dim(hypercube)[1]){
        if(length(which(hypercube[indR,,indL]!=0))>0){
          index = which(hypercube[indR,,indL]!=0)
          idx.aux = hypercube[indR,index,indL]
          if (length(idx.aux[idx.aux == 0])>0) {print(idx.aux)}
          idxSelected <- tryCatch(                       
            expr = {                     
              idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
            },
            
            error = function(e){
              message(paste("There was an error for: (variable id) "),paste(hypercube[indR,,indL],collapse = ","))
              return(idxSelected = list('idxSelected' = NULL))
            }
          )
          
          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,index[idxSelected$idxSelected],indL]=hypercubeSelect[indR,index[idxSelected$idxSelected],indL]+1
          }
        }
      } # indR
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[,indC,indL]!=0))>0){
          index = which(hypercube[,indC,indL]!=0)
          idx.aux = hypercube[index,indC,indL]
          idxSelected <- tryCatch(                       
            expr = {                     
              idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
            },
            
            error = function(e){
              message(paste("There was an error for: (variable id) "),paste(hypercube[,indC,indL],collapse = ","))
              return(idxSelected = list('idxSelected' = NULL))
            }
          )
          
          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[index[idxSelected$idxSelected],indC,indL]=hypercubeSelect[index[idxSelected$idxSelected],indC,indL]+1
          }
        }
      } # indC
    } # indL
    for(indR in 1:dim(hypercube)[1]){
      for(indC in 1:dim(hypercube)[2]){
        if(length(which(hypercube[indR,indC,]!=0))>0){
          index = which(hypercube[indR,indC,]!=0)
          idx.aux = hypercube[indR,indC,index]
          idxSelected <- tryCatch(                       
            expr = {                     
              idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
            },
            
            error = function(e){
              message(paste("There was an error for: (variable id) "),paste(hypercube[indR,indC,],collapse = ","))
              return(idxSelected = list('idxSelected' = NULL))
            }
          )
          
          if(length(idxSelected$idxSelected)>0){
            hypercubeSelect[indR,indC,index[idxSelected$idxSelected]]=hypercubeSelect[indR,indC,index[idxSelected$idxSelected]]+1
          }
        }
      } # indC
    } # indR
    
    setSelected3Times = hypercube[which(hypercubeSelect>2, arr.ind = TRUE)]
    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]
    
    numSelected3times=length(setSelected3Times)
    numSelected2times=length(setSelected2Times)
    numSelected1times=length(setSelected1Times)
    
    Matrix.Selection[[paste('Hypercube with dim',3)]] = c(numSelected1times,numSelected2times,numSelected3times)
    names(Matrix.Selection[[paste('Hypercube with dim',3)]]) = c('numSelected1','numSelected2','numSelected3')
    
    if(numSelected1times>1){
      List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected1')]] = setSelected1Times
    } else{
      List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected1')]] = c(0,ifelse(length(setSelected1Times)==0,0,c(setSelected1Times)))
      List.Selection$`Hypercube with dim 3`$numSelected1 = List.Selection$`Hypercube with dim 3`$numSelected1[-1]
    }
    
    if(numSelected2times>1){
      List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected2')]] = setSelected2Times
    } else{
      List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected2')]] = c(0,ifelse(length(setSelected2Times)==0,0,c(setSelected2Times)))
      List.Selection$`Hypercube with dim 3`$numSelected2 = List.Selection$`Hypercube with dim 3`$numSelected2[-1]
    }
    
    if(numSelected3times>1){
      List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected3')]] = setSelected3Times
    } else{
      List.Selection[[paste('Hypercube with dim',3)]][[paste('numSelected3')]] = c(0,ifelse(length(setSelected3Times)==0,0,c(setSelected3Times)))
      List.Selection$`Hypercube with dim 3`$numSelected3 = List.Selection$`Hypercube with dim 3`$numSelected3[-1]
    }
    
    aux.dmHC3 <- 'Y'
    
  }
  
  ########## case in which dmHC=2 ##########
  if(dmHC==2 | aux.dmHC3=='Y'){
    
    if(signif.Default==TRUE & highest.dmHC==2){
      aux.signif = 2
    } else if(signif.Default==TRUE & highest.dmHC>2){
      aux.signif = 0.01
    } else{
      aux.signif = vector.signif[2]
    }
    
    if(dmHC==2){
      dimHC = ceiling(d^(1/2))
      nearestHypercube=dimHC^2
      remainderHC=nearestHypercube-d
      if(!is.null(seed.HC)){
        set.seed(seed.HC)
        hypercube=array(sample(c((1:d),rep(0,remainderHC))),dim=c(dimHC,dimHC))
      } else{
        hypercube=array(c((1:d),rep(0,remainderHC)),dim=c(dimHC,dimHC))
      }} else if(dmHC>2){
        dimHC = ceiling(numSelected1times^(1/2))
        nearestHypercube=dimHC^2
        remainderHC=nearestHypercube-length(setSelected2Times)
        if(!is.null(seed.HC)){
          set.seed(seed.HC)
          hypercube=array(sample(c(setSelected2Times,rep(0,remainderHC))),dim=c(dimHC,dimHC))
        } else{
          hypercube=array(c(setSelected2Times,rep(0,remainderHC)),dim=c(dimHC,dimHC))
        }
      }
    
    ### intermediate step to handle 0`s here
    hypercubeSelect=array(0,dim=c(dimHC,dimHC))
    
    if(all(dim(hypercube)==0)){
      stop(paste('No variables selected at stage 2, increase p-value.'))
    }
    
    for(indR in 1:dim(hypercube)[1]){
      if(length(which(hypercube[indR,]!=0))>0){
        index = which(hypercube[indR,]!=0)
        idx.aux = hypercube[indR,index]
        idxSelected <- tryCatch(                       
          expr = {                     
            idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
          },
          
          error = function(e){
            message(paste("There was an error for: (variable id) "),paste(hypercube[indR,],collapse = ","))
            return(idxSelected = list('idxSelected' = NULL))
          }
        )
        
        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[indR,index[idxSelected$idxSelected]]=hypercubeSelect[indR,index[idxSelected$idxSelected]]+1
        }
      }
    } # indR
    for(indC in 1:dim(hypercube)[2]){
      if(length(which(hypercube[,indC]!=0))>0){
        index = which(hypercube[,indC]!=0)
        idx.aux = hypercube[index,indC]
        idxSelected <- tryCatch(                       
          expr = {                     
            idxSelected = lmensspRoutine(idx.aux,significance=aux.signif)
          },
          
          error = function(e){
            message(paste("There was an error for: (variable id) "),paste(hypercube[,indC],collapse = ","))
            return(idxSelected = list('idxSelected' = NULL))
          }
        )
        
        if(length(idxSelected$idxSelected)>0){
          hypercubeSelect[index[idxSelected$idxSelected],indC]=hypercubeSelect[index[idxSelected$idxSelected],indC]+1
        }
      }
    } # indC
    
    setSelected2Times = hypercube[which(hypercubeSelect>1, arr.ind = TRUE)]
    setSelected1Times = hypercube[which(hypercubeSelect>0, arr.ind = TRUE)]
    
    numSelected2times=length(setSelected2Times)
    numSelected1times=length(setSelected1Times)
    
    Matrix.Selection[[paste('Hypercube with dim',2)]] = c(numSelected1times,numSelected2times)
    names(Matrix.Selection[[paste('Hypercube with dim',2)]]) = c('numSelected1','numSelected2')
    
    if(numSelected1times>1){
      List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected1')]] = setSelected1Times
    } else{
      List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected1')]] = c(0,ifelse(length(setSelected1Times)==0,0,c(setSelected1Times)))
      List.Selection$`Hypercube with dim 2`$numSelected1 = List.Selection$`Hypercube with dim 2`$numSelected1[-1]
    }
    
    if(numSelected2times>1){
      List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected2')]] = setSelected2Times
    } else{
      List.Selection[[paste('Hypercube with dim',2)]][[paste('numSelected2')]] = c(0,ifelse(length(setSelected2Times)==0,0,c(setSelected2Times)))
      List.Selection$`Hypercube with dim 2`$numSelected2 = List.Selection$`Hypercube with dim 2`$numSelected2[-1]
    }
    
    aux.dmHC2 <- 'Y' #readline(cat("Reduction of dimension 2 done!", "\n", length(setSelected2Times),
    
  }
  
  return(list("seed.HC" = seed.HC, "Matrix.Selection" = Matrix.Selection, "List.Selection" = List.Selection))

}

###################################################################################################


Reduction.Phase_out <- Reduction.Phase(input_data)

save(Reduction.Phase_out, file = 'Reduction.Phase_out.RData')
   