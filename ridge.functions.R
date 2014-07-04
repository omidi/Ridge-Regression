# -----------------------------------------------------
# Collection of ridge regression functions
# -----------------------------------------------------

# fast.svd: Taken from corpcor package
# source('fast.svd.R')

ridge.qr = function(X,Y,lambda,center.cols=FALSE,center.rows=FALSE,return.Yhat=FALSE,return.Yvar=FALSE,return.Winv=FALSE) {
  # ridge regression by QR decomposition of X^tX (faster than svd)
  
  # control centering for genes/samples independently
  if (center.cols)
    {
      Y = sweep(Y,2,colMeans(Y))
      X = sweep(X,2,colMeans(X))
    }
  if (center.rows & ncol(Y)>1)
    {
      Y = sweep(Y,1,rowMeans(Y))
    }

  ngene = nrow(X)
  nmotif = ncol(X)  
  nsample = ncol(Y)
  
  XXl = crossprod(X) + diag(ngene*lambda,nmotif,nmotif)
  QR = qr(XXl)
  XY = crossprod(X,Y)
  Ahat = qr.coef(QR,XY) # predicted activities
  Yhat = X %*% Ahat # predicted expression values
  resid = Y-Yhat
  Chi2 = colSums(resid^2)   
  Yvar = Chi2/ngene # sigma_sample; un-explained variance
  
  # log marginal likelihood P(E|N)
  f = nsample*(ngene-nmotif)/2
  log.det = sum(log(abs(diag(QR$qr)))) # determinant of NtN+I*lambda based on QR 
  log.marg.lik = nsample*0.5*(nmotif*log(2*pi) - log.det) + (f-1)*log(2) + lgamma(f) - f*log(sum(Chi2)) 
      
  # compute Z-score
  R = qr.R(QR)
  Rinv = backsolve(R,diag(1,nrow(R),nrow(R))) # R inverse
  Winv = Rinv %*% t(qr.Q(QR)) # (NtN + I*lambda)^-1
  AhatSE = sqrt(diag(Winv) %x% t(Yvar))
  Zscore = Ahat/AhatSE
  combined.Zscore = sqrt(rowMeans(Zscore^2))
  
  fit = list(Ahat=Ahat,Zscore=Zscore,combined.Zscore=combined.Zscore,log.marg.lik=log.marg.lik,stats=ridge.stats(Y,resid),np=ngene,ns=ncol(Y))
  if (return.Yhat)
    fit$Yhat = Yhat
  if (return.Yvar)
    fit$Yvar = Yvar
  if (return.Winv)
    {
      rownames(Winv) = colnames(X)
      colnames(Winv) = colnames(X)
      fit$Winv = Winv
    }
  
  return(fit)
}

# -----------------------------------------------------

ridge.qr.coef = function(X,Y,lambda) {
  # ridge regression by qr decomposition of X^tX (faster than svd)
  # coefficients only
  ngene = nrow(X)
  nmotif = ncol(X)
  XXl = crossprod(X) + diag(ngene*lambda,nmotif,nmotif)
  QR = qr(XXl)
  Ahat = qr.coef(QR,crossprod(X,Y)) # predicted activities
  return(Ahat)
}

# -----------------------------------------------------

ridge.qr.coef2 = function(XtX,XtY,ngeneLambda) {
  # fastest ridge regression with pre-computed cross products
  # NOTE: lambda has to be already scaled by the number of genes!
  QR = qr(XtX + diag(ngeneLambda,ncol(XtX),ncol(XtX)))
  Ahat = qr.coef(QR,XtY) # predicted activities
  return(Ahat)
}

# -----------------------------------------------------
  
ridge.svd = function(X,Y,lambda,center.cols=FALSE,center.rows=FALSE,return.Yhat=FALSE,return.Yvar=FALSE,return.Winv=FALSE) {
  # Ridge regression by SVD
  
  # control centering for genes/samples independently
  if (center.cols) {
    Y = sweep(Y,2,colMeans(Y))
    X = sweep(X,2,colMeans(X))
  }
  if (center.rows & ncol(Y)>1) {
    Y = sweep(Y,1,rowMeans(Y))
  }

  nmotif = ncol(X)
  ngene = nrow(Y)
  nsample = ncol(Y)
  
  Xs = fast.svd(X)
  rhs = crossprod(Xs$u,Y)
  dia = Xs$d/(Xs$d^2 + ngene*lambda)
  Ahat = Xs$v %*% drop(dia * rhs)  
  Yhat = X %*% Ahat
  resid = Y-Yhat
  Chi2 = colSums(resid^2)   
  Yvar = Chi2/ngene # sigma_sample; un-explained variance  

  # log marginal likelihood P(E|N)
  f = nsample*(ngene-nmotif)/2
  log.det = sum(log(Xs$d^2+ngene*lambda)) # determinant of NtN+I*lambda based on SVD 
  log.marg.lik = nsample*0.5*(nmotif*log(2*pi) - log.det) + (f-1)*log(2) + lgamma(f) - f*log(sum(Chi2))
  
  # compute Z-score
  Winv = tcrossprod(sweep(Xs$v,2,1/(Xs$d^2 + ngene*lambda),FUN='*'),Xs$v)
  AhatSE = sqrt(diag(Winv) %x% t(Yvar))
  Zscore = Ahat/AhatSE
  rownames(Ahat) = colnames(X)
  rownames(Zscore) = colnames(X)
  combined.Zscore = sqrt(rowMeans(Zscore^2))

  fit = list(Ahat=Ahat,Zscore=Zscore,combined.Zscore=combined.Zscore,log.marg.lik=log.marg.lik,stats=ridge.stats(Y,resid),np=ngene,ns=ncol(Y))
  
  if (return.Yhat)
    fit$Yhat = Yhat
  if (return.Yvar)
    fit$Yvar = Yvar
  if (return.Winv) {
    rownames(Winv) = colnames(X)
    colnames(Winv) = colnames(X)
    fit$Winv = Winv
  }
  
  return(fit)
}

# -----------------------------------------------------

ridge.svd.coef2 = function(X,y,lambda) {
  # Ridge regression by SVD: coefficients only
  # Allows for single y and multiple lambda
  Xs <- fast.svd(X) 
  rhs <- t(Xs$u) %*% y 
  d <- Xs$d
  div <- d^2 + rep(length(y)*lambda, rep(length(d),length(lambda))) # lambda*dx x 1
  a <- drop(d * rhs)/div
  dim(a) <- c(length(d), length(lambda))
  Ahat <- Xs$v %*% a
  ## Ehat <- X %*% Ahat
  ## resid <- sweep(Ehat,1,y)
  return(Ahat)
}

# -----------------------------------------------------

ridge.cv.bisect = function(X,Y,center.cols=FALSE,center.rows=FALSE,lambda.bnd=NULL,k.groups=NULL,k=10,compute.optimum=TRUE,print.it=TRUE) {
  # Estimate coefficients by ridge regression;
  # choose optimal lambda by k-fold cross-validation
  # bisection method: chooses an optimal lambda for MULTIPLE conditions!

  # this pre-processing can be done outside of the cross-validation loop
  if (center.rows & ncol(Y)>1) {
    Y = sweep(Y,1,rowMeans(Y))
  }
    
  # guess sensible upper and lower boundaries for lambda: taken from parcor package 
  if (is.null(lambda.bnd)) {
    nr = nrow(X) - floor(nrow(X)/k) # nrows in Xtrain
    lambda.bnd = 10^c(-10,-1) * nr * ncol(X)
  }

  if (is.null(k.groups)) {
    # randomly choose k groups
    k.groups = split(sample(1:nrow(X)), rep(1:k, length = nrow(X)))
  }
  
  # pre-compute cross products to speed up optimization
  if (print.it)
    cat(sprintf('Pre-computing cross-products of training sets...\n'))
  XtX = list(); XtY = list(); cmx = list(); cmy = list()
  for (i in 1:length(k.groups)) {
    omit = k.groups[[i]]
    Xtrain = X[-omit, , drop = FALSE]
    Ytrain = Y[-omit, , drop = FALSE]
    # column centering: this pre-processing has to be done for every
    # training set independently (i.e. "within" the crossvalidation
    # loop)
    if (center.cols) {
      cmy[[i]] = colMeans(Ytrain)
      cmx[[i]] = colMeans(Xtrain)
      Ytrain = sweep(Ytrain,2,cmy[[i]])
      Xtrain = sweep(Xtrain,2,cmx[[i]])
    } 
    XtX[[i]] = crossprod(Xtrain)
    XtY[[i]] = crossprod(Xtrain,Ytrain)
  }

  # now optimize over lambda
  opt = optimize(
    function(lambda,X,Y,XtX,XtY,cmx,cmy,k.groups) { 
      # perform k-fold cross validation
      cv = 0
      for (i in 1:length(k.groups)) {
      	# ridge based on QR
        QR = qr(XtX[[i]] + diag(nrow(X)*lambda,ncol(X),ncol(X)))
        Ahat = qr.coef(QR,XtY[[i]]) # predicted activities        
        # apply to test data
        omit = k.groups[[i]]
        Xtest = X[omit, , drop = FALSE]
        Ytest = Y[omit, , drop = FALSE]
        if (center.cols) {
          Xtest = sweep(Xtest,2,cmx[[i]])
          Ytest = sweep(Ytest,2,cmy[[i]])
        } 
        resid = (Xtest %*% Ahat) - Ytest
        cv = cv + mean(resid^2)
      }
      if (print.it)
        cat(sprintf('Lambda: %.4g, CV: %.6g\n',lambda,cv))
      return(cv/length(k.groups))
    }, interval=lambda.bnd, X, Y, XtX, XtY, cmx, cmy, k.groups)
  lambda.opt = opt$minimum
  cv.opt = opt$objective
 
  # compute activities for optimal lambda over all data
  if (compute.optimum) {
    r = ridge.qr(X,Y,lambda.opt,center.cols=center.cols,center.rows=center.rows)
    return(list(Ahat = r$Ahat,
                Zscore = r$Zscore,
                combined.Zscore = r$combined.Zscore,
                Ehat = r$Yhat,
                Evar = r$Yvar,
                stats = r$stats,
                lambda.opt = lambda.opt,
                cv.opt = cv.opt))
  } else {
    return(list(lambda.opt = lambda.opt,
                cv.opt = cv.opt))
  }
}

# -----------------------------------------------------

ridge.gcv = function(X,Y,XtX=NULL,Xs=NULL,center.cols=FALSE,center.rows=FALSE,lambda.bnd=NULL,compute.optimum=TRUE,return.resid=FALSE,print.it=TRUE) {
  # Choose optimal lambda by GENERALIZED cross-validation
  if (center.cols) {
    Y = sweep(Y,2,colMeans(Y))
    X = sweep(X,2,colMeans(X))
  }
  if (center.rows & ncol(Y)>1) {
    Y = sweep(Y,1,rowMeans(Y))
  }

  # guess sensible upper and lower boundaries for lambda: taken from parcor package 
  if (is.null(lambda.bnd)) {
    lambda.bnd = 10^c(-12,-6) * nrow(X) * ncol(X)
  }

  # pre-compute: fast SVD + RHS
  if (is.null(Xs)) {
    if (print.it)
      cat(sprintf('SVD of XtX...\n'))      
    if (is.null(XtX)) { 
      Xs = svd(crossprod(X),nu=0)
      Xs$d = sqrt(Xs$d)
      Xs$u = X %*% sweep(Xs$v,2,1/Xs$d,FUN='*') #Xs$v %*% diag(1/Xs$d, nrow=ncol(XtX))
    } else {
      Xs = svd(XtX,nu=0)
      Xs$d = sqrt(Xs$d)
      Xs$u = X %*% sweep(Xs$v,2,1/Xs$d,FUN='*')
    }
  } 
  rhs = crossprod(Xs$u,Y)

  # now optimize over lambda
  opt = optimize(
    function(lambda,Y,Xs,rhs,print.it) {
      dia = Xs$d^2/(Xs$d^2 + nrow(Y)*lambda) # of Hat matrix
      DF = sum(dia) # degrees of freedom
      resid = Y - sweep(Xs$u,2,dia,FUN='*') %*% rhs
      GCV = sum((resid/(nrow(Y)-DF))^2) # generalized cross validation error
      if (print.it)
        cat(sprintf('Lambda=%g\tDoF=%g\tGCV=%g\n',lambda,DF,GCV))      
      return(GCV)
    },
    lambda.bnd,Y,Xs,rhs,print.it)
  lambda.opt = opt$minimum
  gcv.opt = opt$objective
  
  # compute activities for optimal lambda over all data
  if (compute.optimum) {
     dia = Xs$d/(Xs$d^2 + nrow(Y)*lambda.opt)
     Ahat = Xs$v %*% drop(dia * rhs)  
     resid = Y - X %*% Ahat
     Chi2 = colSums(resid^2)   
     Yvar = Chi2/nrow(Y) # sigma_sample; un-explained variance  
     # compute Z-score  
     Winv = tcrossprod(sweep(Xs$v,2,1/(Xs$d^2 + nrow(Y)*lambda.opt),FUN='*'),Xs$v)
     AhatSE = sqrt(diag(Winv) %x% t(Yvar))
     Zscore = Ahat/AhatSE
     rownames(Ahat) = colnames(X)
     rownames(Zscore) = colnames(X)
     combined.Zscore = sqrt(rowMeans(Zscore^2))
     res = list(Ahat = Ahat,
                AhatSE = AhatSE,
                Zscore = Zscore,
                combined.Zscore = combined.Zscore,
                Evar = Yvar,
                stats = ridge.stats(Y,resid),
                lambda.opt = lambda.opt,
                gcv.opt = gcv.opt)
     if (return.resid)
       res$resid = resid
     return(res)
  } else {
    return(list(lambda.opt = lambda.opt,
                gcv.opt = gcv.opt))
  }
}

# -----------------------------------------------------

ridge.gcv.2lambda = function(X,Y,lambda.idx,XtX=NULL,Xs=NULL,center.cols=FALSE,center.rows=FALSE,lambda.bnd=NULL,compute.optimum=TRUE,print.it=TRUE) {
  # Choose two optimal lambdas for miRNA and TFs seperately by GENERALIZED cross-validation
  if (center.cols) {
    Y = sweep(Y,2,colMeans(Y))
    X = sweep(X,2,colMeans(X))
  }
  if (center.rows & ncol(Y)>1) {
    Y = sweep(Y,1,rowMeans(Y))
  }

  # pre-compute: fast SVD + RHS
  if (is.null(Xs)) {
    if (print.it)
      cat(sprintf('SVD of XtX...\n'))      
    if (is.null(XtX)) { 
      Xs = svd(crossprod(X),nu=0)
      Xs$d = sqrt(Xs$d)
      Xs$u = X %*% sweep(Xs$v,2,1/Xs$d,FUN='*') #Xs$v %*% diag(1/Xs$d, nrow=ncol(XtX))
    } else {
      Xs = svd(XtX,nu=0)
      Xs$d = sqrt(Xs$d)
      Xs$u = X %*% sweep(Xs$v,2,1/Xs$d,FUN='*')
    }
  } 
  rhs = crossprod(Xs$u,Y)

  # now optimize over lambda
  lambda.bnd = 10^c(-10,-6) * nrow(X) * ncol(X)
  lambda.init = rep(lambda.bnd[1] + (lambda.bnd[2]-lambda.bnd[1])/2,2)
  opt = optim(lambda.init,
    function(lambda,X,Y,Xs,rhs,lambda.idx,print.it) {
      lambda.vec = rep(0,ncol(X))
      lambda.vec[lambda.idx==1] = nrow(Y)*lambda[1]
      lambda.vec[lambda.idx==2] = nrow(Y)*lambda[2]
      dia = Xs$d^2/(Xs$d^2 + lambda.vec) # of Hat matrix
      DF = sum(dia) # degrees of freedom
      resid = Y - sweep(Xs$u,2,dia,FUN='*') %*% rhs
      GCV = sum((resid/(nrow(X)-DF))^2) # generalized cross validation error
      if (print.it)
        cat(sprintf('Lambda1=%g\tLambda2=%g\tDoF=%g\tGCV=%g\n',lambda[1],lambda[2],DF,GCV))      
      return(GCV)
    },NULL,X,Y,Xs,rhs,lambda.idx,print.it)
  lambda.opt = opt$par
  gcv.opt = opt$value
  
  # compute activities for optimal lambda over all data
  if (compute.optimum) {
    lambda.vec = rep(0,ncol(X))
    lambda.vec[lambda.idx==1] = nrow(Y)*lambda.opt[1]
    lambda.vec[lambda.idx==2] = nrow(Y)*lambda.opt[2]
    dia = Xs$d/(Xs$d^2 + lambda.vec)
    Ahat = Xs$v %*% drop(dia * rhs)  
    resid = Y - X %*% Ahat
    Chi2 = colSums(resid^2)   
    Yvar = Chi2/nrow(Y) # sigma_sample; un-explained variance  
                                        # compute Z-score  
    Winv = tcrossprod(sweep(Xs$v,2,1/(Xs$d^2 + lambda.vec),FUN='*'),Xs$v)
    AhatSE = sqrt(diag(Winv) %x% t(Yvar))
    Zscore = Ahat/AhatSE
    rownames(Ahat) = colnames(X)
    rownames(Zscore) = colnames(X)
    combined.Zscore = sqrt(rowMeans(Zscore^2))
    return(list(Ahat = Ahat,
                Zscore = Zscore,
                combined.Zscore = combined.Zscore,
                Evar = Yvar,
                stats = ridge.stats(Y,resid),
                lambda.opt = lambda.opt,
                gcv.opt = gcv.opt))
  } else {
    return(list(lambda.opt = lambda.opt,
                gcv.opt = gcv.opt))
  }
}

# -----------------------------------------------------

target.prediction = function(Ahat,centered.N,centered.E,N.colMeans) {
    # For every promotor p, motif m and N_pm > 0, calculate the
    # likelihood ratio of E_p with and without site m 
    # assuming Ahat stays un-changed by droping site m (checked, OK)   

    Ehat = centered.N %*% Ahat # predicted expression using all motif activities
    chi2 = rowSums((centered.E-Ehat)^2)
    chi2.mean = sum(chi2)/prod(dim(centered.E))
    # likelihood ratio: LR < 0 drop site, LR > 0 keep site
    LLR = matrix(NA,nrow(centered.E),ncol(centered.N),dimnames=list(rownames(centered.E),colnames(centered.N)))
    for (m in 1:ncol(centered.N)) {
      Nm = centered.N[,m]+N.colMeans[m]
      tidx = which(Nm>0) # target genes of motif m
      chi2.drop = rowSums((centered.E[tidx,] - (Ehat[tidx,] - Nm[tidx] %o% Ahat[m,]))^2) # the same using outer product
      LLR[tidx,m] = (chi2.drop-chi2[tidx])/chi2.mean
    }
    return(LLR) 
}

# -----------------------------------------------------
# helper functions
# -----------------------------------------------------

ridge.stats = function(Y,resid) {
  # assumes Y is centered across columns
  n = nrow(Y)
  SSE = colSums(resid^2) # sum of squared error around estimate
  if (is.null(dim(Y)))
    SEE = sum(Y^2) 
  else
    SEE = colSums(Y^2) # sum of squared error around mean (=0)
  R2 = (SEE-SSE)/SEE # % explained variance
  return(list(MSE=SSE/n,R2=R2))
}

# -----------------------------------------------------

group.zscore = function(Ahat,M,Yvar,idx)
{
  # compute marginal Z-score of a group of TFs defined by a logical index
  
  # Application of Jaynes formula:
  # Given the posterior P ~ exp(a M a); with M = (NtN+lambda)
  # Split activity vector a up into:
  # 1. u the interesting marginals (a[idx])
  # 2. w the unintersting variables to integrate over (a[!idx])
  # Then write M = (U0 V
  #                 Vt W0)
  # Now the maribal posterior is P' ~ exp(u U u) with U = U0 - V * W0^-1 * V^t
  
  # compute the inverse of W0
  QR = qr(M[!idx,!idx,drop=FALSE]) #
  R = qr.R(QR)
  Rinv = backsolve(R,diag(1,nrow(R),nrow(R)))
  W0inv = Rinv %*% t(qr.Q(QR))
  
  V = M[idx,!idx]
  U = M[idx,idx] - V %*% W0inv %*% t(V)
  
  # compute sample Z-scores of the group
  Zscore = rep(NA,ncol(Ahat))
  for (s in 1:ncol(Ahat)) {
    Zscore[s] = sum(sqrt(diag(U/Yvar[s])) * Ahat[idx,s,drop=FALSE])
  }
  return(Zscore)
}

# -----------------------------------------------------

## group.zscore = function(Ahat,Winv,Yvar,idx)
## {
##   # compute marginal Z-score of a group of TFs defined by idx
##   # for each sample s
    
##   # compute the inverse marginal covariance matrix
##   QR = qr(Winv[idx,idx,drop=FALSE])
##   R = qr.R(QR)
##   Rinv = backsolve(R,diag(1,nrow(R),nrow(R)))
##   Cinv = Rinv %*% t(qr.Q(QR))
  
##   # compute sample Z-scores of the group
##   Zscore = rep(NA,length(Yvar))
##   for (s in 1:ncol(Ahat)) {
## #    Zscore[s] = sqrt(t(Ahat[idx,s,drop=FALSE]) %*% (Cinv/Yvar[s]) %*% Ahat[idx,s,drop=FALSE])
##       Zscore[s] = sum(sqrt(Cinv/Yvar[s]) %*% Ahat[idx,s,drop=FALSE])
##   }
##   return(Zscore)
## }

# -----------------------------------------------------

all.groups.zscore = function(Ahat,M,Yvar,groups) {
  # compute marginal Z-score of all groups
  # groups is a list

  Zscore = matrix(0,length(groups),ncol(Ahat))
  rownames(Zscore) = names(groups)
  colnames(Zscore) = colnames(Ahat)
  idx = rep(FALSE,nrow(Ahat))
  names(idx) = rownames(Ahat)
  for (g in 1:length(groups)) {
    print(g)
    idx = FALSE # reset
    idx[groups[[g]]] = TRUE
    Zscore[g,] = group.zscore(Ahat,M,Yvar,idx)
  }

  return(Zscore)
}

# -----------------------------------------------------
# experimental ridge versions, generalized SVD, feature selection based on log marginal likelihood
# -----------------------------------------------------

ridge.edf = function(X,Y,EDF) {
  # perfom a ridge regression for a given effective degree of freedom (EDF)
  # i.e. optimize lambda such that rank(H) machtes the EDF, where Yhat = H * Y
    
  nm = ncol(X)
  nr = nrow(Y)
  ns = ncol(Y)
    
  Xs = fast.svd(X)
  lambda.bnd = c(0,10^2)
  opt = optimize(
    function(lambda,EDF,d,nr,nm) { 
      DF = sum(d^2/(d^2 + nr*lambda))
      return((DF-EDF)^2)
    },
    lambda.bnd,EDF,Xs$d,nr,nm)
  lambda = opt$minimum
  # solve with optimal lambda
  dia = Xs$d/(Xs$d^2 + nr*lambda)
  rhs = crossprod(Xs$u,Y)
  Ahat = Xs$v %*% drop(dia * rhs)
  DF = sum(Xs$d*dia)

  # compute residuals
  Yhat = X %*% Ahat
  resid = Y-Yhat
  Yvar = colSums(resid^2)/nr # sigma_sample
  
  # compute Z-score  
  W = Xs$v %*% diag(1/(Xs$d^2 + nr*lambda),nm,nm) %*% t(Xs$v)
  AhatSE = sqrt(diag(W) %x% t(Yvar))
  Zscore = Ahat/AhatSE
  rownames(Ahat) = colnames(X)
  rownames(Zscore) = colnames(X)
  combined.Zscore = sqrt(rowMeans(Zscore^2))

  fit = list(Ahat=Ahat,Zscore=Zscore,combined.Zscore=combined.Zscore,
    stats=ridge.stats(Y,resid),lambda=lambda,DF=DF)
  
  return(fit)
}

# -----------------------------------------------------

pca.ridge = function(X,Y,lambda,sY=NULL) {
  # do ridge regression on the principal components of Y
  X = sweep(X,2,colMeans(X)) 
  Y = sweep(Y,2,colMeans(Y))
  Y = sweep(Y,1,rowMeans(Y)) 

  if (is.null(sY))
    sY = svd(Y)
  
  r = ridge.qr(X,sY$u,lambda,center.cols=FALSE,center.rows=FALSE)

  return(list(s=sY,r=r))
}


# -----------------------------------------------------

decorrelate = function(X)
  {
    # decorrelate columns of X using the inverse correlation matrix
    X = scale(X)
    s = svd(crossprod(X))
    Pi12 = s$v %*% diag(sqrt(1/s$d),ncol(s$v),ncol(s$v)) # P^-1/2
    Xd = X %*% Pi12 
    return(Xd)
  }

# -----------------------------------------------------

compute.lml = function(X,Y,XtX,XtY,lambda,active)
{
  m = sum(active)
  ## print(active)
  ## print(dim(XtX[active,active,drop=FALSE]))
  ## print(dim(diag(rep(nrow(X)*lambda,m),m,m)))
  # compute the log marginal likelihood P(E|N)
  QR = qr(XtX[active,active,drop=FALSE] + diag(rep(nrow(X)*lambda,m),m,m))
  Ahat = qr.coef(QR,XtY[active,,drop=FALSE])
  Chi2 = sum((Y - tcrossprod(X[,active,drop=FALSE],t(Ahat)))^2) 
  f = ncol(Y)*(nrow(Y)-nrow(Ahat))/2
  log.det = sum(log(abs(diag(QR$qr)))) # determinant based on QR 
  log.marg.lik = ncol(Y)*0.5*(nrow(Ahat)*log(2*pi) - log.det) + (f-1)*log(2) + lgamma(f) - f*log(Chi2) 
  return(log.marg.lik)
}

# -----------------------------------------------------

lml.select = function(X,Y,center.cols=FALSE,center.rows=FALSE,start.full=FALSE) {
  # do a Bayesian stepwise selection (add/drop features) untill no features can be added (log[BF]<0)
  if (center.cols) {
    Y = scale(Y,scale=FALSE)
    X = scale(X,scale=FALSE)
  }
  if (center.rows & ncol(Y)>1) {
    Y = t(scale(t(Y),scale=FALSE))
  }
  XtX = crossprod(X)
  XtY = crossprod(X,Y)
  # initialize
  if (start.full) {
    active = rep(TRUE,ncol(X)) # inital model is empty
    r = ridge.gcv(X,Y,XtX,print.it=FALSE)
    lambda = r$lambda.opt
  } else {
    active = rep(FALSE,ncol(X)) # inital model is empty
    lambda = 0
  }
  lml = compute.lml(X,Y,XtX,XtY,lambda,active) 
  ok = TRUE
  Zscore = c()
  log.BF = c()
  lambdas = c(lambda)
  gcv.error = c()
  while (ok) {
    # determine lml of all neighboring models
    alml = rep(-Inf,ncol(X))
    for (n in 1:ncol(X)) {
      active[n] = !active[n] # add/drop
      alml[n] = compute.lml(X,Y,XtX,XtY,lambda,active)
      active[n] = !active[n] # flip back
    }
    # which is the best neighboring model?
    lbf = alml-lml
    midx = which.max(lbf)
    ok = lbf[midx]>0
    if (ok) {
      # update model
      active[midx] = !active[midx]
      lml = alml[midx]
      cat(sprintf("%s %s,\tlog[BF] = %g\n",colnames(X)[midx],ifelse(active[midx],"added","dropped"),lbf[midx]))
      # update lambda
      r = ridge.gcv(X[,active,drop=FALSE],Y,XtX[active,active,drop=FALSE],print.it=FALSE)
      lambda = r$lambda.opt
      gcv.error = c(gcv.error,r$gcv.opt)
      # remember lambda, z-score, log.BF
      lambdas = c(lambdas,lambda)
      z = rep(0,ncol(X))
      z[active] = r$combined.Zscore
      Zscore = cbind(Zscore,z)
      log.BF = cbind(log.BF,lbf)
    } else {
      cat(sprintf("\nIteration stoped at max(log[BF]) = %g\nwith %d features included\n",lbf[midx],sum(active)))
    }
  }
  rownames(Zscore) = colnames(X)
  rownames(log.BF) = colnames(X)
  active = colnames(X)[active]
  return(list(active=active,lambda=lambdas,gcv.error=gcv.error,Zscore=Zscore,log.BF=log.BF))
}

plot.selection.trace = function(Zscore,log.BF) {
  idx = which(rowSums(Zscore)>0)
  par(mfrow=c(3,1))
  plot(0,0,type='n',xlim=c(0,ncol(Zscore)),ylim=c(0,max(Zscore[idx,])),xlab='Steps',ylab='Combined Z-Score')
  for (n in 1:length(idx)) {
    lines(1:ncol(Zscore),Zscore[idx[n],],lty=ifelse(Zscore[idx[n],ncol(Zscore)]>0,1,2))
  }
  barplot(log.BF[Zscore[,ncol(log.BF)]>0,ncol(log.BF)],las=2)
  plot(0,0,type="n",axes=FALSE,xlab='',ylab='')
}
