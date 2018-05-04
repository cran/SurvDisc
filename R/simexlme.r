

simexlme= function (model, model.model, SIMEXvariable, respvar,grpvar,corform,measurement.error, measurement.error.resp,
                     lambda = c(0.5, 1, 1.5, 2),B = 100, fitting.method = "quadratic", jackknife.estimation = "quadratic")
{
#SIMEX algorithm for lme models where the response is the change from baseline and the baseline value is the SIMEXvariable

#these first 3 functions construct.s, extract.covmat, and fitnls (there, named fit.nls) are copied directly from the SIMEX package

  construct.s = function (ncoef, lambda, fitting.method, extrapolation = NULL)
  {
    nl <- length(lambda)
    switch(fitting.method, quad = ngamma <- 3, line = ngamma <- 2,
           nonl = ngamma <- 3, logl = ngamma <- 2, log2 = ngamma <- 2)
    null.mat <- matrix(0, nrow = ngamma, ncol = ncoef)
    s <- list()
    for (j in 1:ncoef) {
      switch(fitting.method, quad = gamma.vec <-
               c(-1, -lambda[1],-lambda[1]^2), line = gamma.vec <- c(-1, -lambda[1]),
             nonl = gamma.vec <- c(1, 1/(coef(extrapolation[[j]])[3] +lambda[1]),
                                   -coef(extrapolation[[j]])[2]/(coef(extrapolation[[j]])[3] +
                                                                   lambda[1])^2), logl = gamma.vec <-
               c(exp(coef(extrapolation)[1,j] + coef(extrapolation)[2, j] * lambda[1]),exp(coef(extrapolation)[1, j] +
                                                                                             coef(extrapolation)[2,j] * lambda[1]) * lambda[1]),
             log2 = gamma.vec <- c(exp(coef(extrapolation[[j]])[1] +
                                         coef(extrapolation[[j]])[2] * lambda[1]), exp(coef(extrapolation[[j]])[1] +
                                                                                         coef(extrapolation[[j]])[2] * lambda[1]) * lambda[1]))
      a <- null.mat
      a[, j] <- gamma.vec
      for (i in 2:nl) {
        switch(fitting.method, quad = gamma.vec <- c(-1,-lambda[i], -lambda[i]^2),
               line = gamma.vec <- c(-1,-lambda[i]), nonl = gamma.vec <-
                 c(1, 1/(coef(extrapolation[[j]])[3] +lambda[i]), -coef(extrapolation[[j]])[2]/(coef(extrapolation[[j]])[3] +lambda[i])^2),
               logl = gamma.vec <- c(exp(coef(extrapolation)[1,j] + coef(extrapolation)[2, j] * lambda[i]),
                                     exp(coef(extrapolation)[1, j] +
                                           coef(extrapolation)[2,j] * lambda[i]) * lambda[i]),
               log2 = gamma.vec <- c(exp(coef(extrapolation[[j]])[1] +
                                           coef(extrapolation[[j]])[2] * lambda[i]), exp(coef(extrapolation[[j]])[1] +
                                                                                           coef(extrapolation[[j]])[2] * lambda[i]) * lambda[i]))
        b <- null.mat
        b[, j] <- gamma.vec
        a <- cbind(a, b)
      }
      s[[j]] <- a
    }
    s <- t(matrix(unlist(lapply(s, t), recursive = FALSE), nrow = nl *
                    ncoef, ncol = ngamma * ncoef))
    return(s)
  }


  extract.covmat = function (model)
  {
    type <- class(model)[1]
    sum.model <- summary(model)
    switch(type, glm = covmat <- sum.model$cov.scaled, lm = covmat <- sum.model$cov.unscaled *
             sum.model$sigma^2, gam = covmat <- model$Vp, nls = covmat <- sum.model$cov.unscaled *
             sum.model$sigma^2, lme = covmat <- model$apVar, nlme = covmat <- model$apVar)
    return(covmat)
  }

  fitnls = function (lambda, p.names, estimates)
  {
    extrapolation <- list()
    lambdastar <- c(0, max(lambda)/2, max(lambda))
    for (d in p.names) {
      extrapolation.quad <- lm(estimates[, d] ~ lambda + I(lambda^2))
      a.nls <- predict(extrapolation.quad, newdata = data.frame(lambda = lambdastar))
      gamma.est.3 <- ((a.nls[2] - a.nls[3]) * lambdastar[3] *
                        (lambdastar[2] - lambdastar[1]) - lambdastar[1] *
                        (a.nls[1] - a.nls[2]) * (lambdastar[3] - lambdastar[2]))/
        ((a.nls[1] -a.nls[2]) * (lambdastar[3] - lambdastar[2]) - (a.nls[2] -a.nls[3]) * (lambdastar[2] - lambdastar[1]))
      gamma.est.2 <- ((a.nls[2] - a.nls[3]) * (gamma.est.3 +
                                                 lambdastar[2]) * (gamma.est.3 + lambdastar[3]))/(lambdastar[3] -lambdastar[2])
      gamma.est.1 <- a.nls[1] - (gamma.est.2/(gamma.est.3 +
                                                lambdastar[1]))
      extrapolation[[d]] <- nls(estimates[, d] ~ gamma.1 +
                                  gamma.2/(gamma.3 + lambda), start = list(gamma.1 = gamma.est.1,
                                                                           gamma.2 = gamma.est.2, gamma.3 = gamma.est.3))
    }
    return(extrapolation)
  }





  model.model=model.model[order(model.model[,grpvar]),]
  unqgrp=unique(model.model[,grpvar])
  N.group =length(unqgrp)
  ni.group=rep(0,N.group)
  for (i in 1:N.group) ni.group[i]=sum(model.model[,grpvar]==unqgrp[i])
  asymptotic = F #implemented only for this case with lme
  fitting.method <- substr(fitting.method, 1, 4)
  if (!any(fitting.method == c("quad", "line", "nonl"))) {
    warning("Fitting Method not implemented. Using: quadratic",
            call. = FALSE)
    fitting.method <- "quad"
  }
  if (jackknife.estimation[1] != FALSE)
    jackknife.estimation <- substr(jackknife.estimation,
                                   1, 4)
  if (!any(jackknife.estimation == c("quad", "line", "nonl",
                                     FALSE))) {
    warning("Fitting Method (jackknife) not implemented. Using: quadratic",
            call. = FALSE)
    jackknife.estimation <- "quad"
  }
  if (!is.character(SIMEXvariable[1]))
    stop("SIMEXvariable must be character", call. = FALSE)
  if (any(lambda <= 0)) {
    warning("lambda should not contain 0 or negative values. 0 or negative values will be ignored",
            call. = FALSE)
    lambda <- lambda[lambda >= 0]
  }
  if (!any(names(model) == "x") && asymptotic)
    stop("The option x must be enabled in the naive model for asymptotic variance estimation",
         call. = FALSE)
#  measurement.error <- as.matrix(measurement.error)
  SIMEXdata = model.model
#  if (NROW(measurement.error) != NROW(SIMEXdata) && NROW(measurement.error) ==
#      1) {
#    measurement.error <- matrix(measurement.error, nrow = NROW(SIMEXdata),
#                                ncol = NCOL(measurement.error), byrow = TRUE)
#  }
#  if (NROW(measurement.error) != NROW(SIMEXdata) && NROW(measurement.error) !=
#      1) {
#    stop("NROW(measurement.error) must be either 1 or must take the number of rows of the data used.",
#         call. = FALSE)
#  }
#  if (any(measurement.error < 0))
#    stop("measurement.error is negative", call. = FALSE)
#  any0 <- (apply(measurement.error, 2, all.equal, current = rep(0,
#                                                                times = NROW(measurement.error))) == TRUE)
#  if (sum(any0) > 0)
#    stop("measurement.error is constant 0 in column(s) ",
#         which(any0), call. = FALSE)
#  if (asymptotic == TRUE & ((class(model)[1] != "glm" & class(model)[1] !=
#                             "lm") | dim(unique(measurement.error))[1] != 1))
#    stop("Asymptotic is only implemented for naive models of class lm or glm with homoscedastic measurement error.")
  cl <- match.call()
#  ncoef <- length(model$coefficients)
#  ndes <- dim(model$model)[1]
#  p.names <- names(coef(model))
  ncoef <- length(model$coefficients$fixed)
  ndes <- dim(model.model)[1]
  p.names <- names(model$coefficients$fixed)
  nlambda <- length(lambda)
  estimates <- matrix(data = NA, nlambda + 1, ncoef)
  theta <- matrix(data = NA, B, ncoef)
  colnames(theta) <- p.names
  theta.all <- vector(mode = "list", nlambda)
  if (jackknife.estimation[1] != FALSE) {
    var.exp <- list()
    var.exp[[1]] <- extract.covmat(model)
  }
  if (asymptotic[1]) {
    psi <- matrix(rep(0, ndes * ncoef), ncol = ncoef, nrow = ndes)
    psi <- resid(model, type = "response") * model$x
    PSI <- psi
    am <- list()
    a <- list()
    xi <- model$x
    dh <- rep(1, ndes)
    if (class(model)[1] == "glm")
      dh <- model$family$mu.eta(model$linear.predictors)
    for (k in 1:ndes) a[[k]] <- dh[k] * xi[k, ] %*% t(xi[k,
                                                         ])
    a.mat <- matrix(unlist(a), nrow = length(a), byrow = TRUE)
    ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]), byrow = FALSE)
    am[[1]] <- -ab/ndes
    a <- list()
  }
  estimates[1, ] <- model$coefficients$fixed
  for (i in 1:length(lambda)) {
    if (jackknife.estimation[1] != FALSE)
      variance.est <- matrix(0, ncol = ncoef, nrow = ncoef)
    if (asymptotic[1]) {
      psi <- matrix(0, ncol = ncoef, nrow = ndes)
      a <- list()
      for (k in 1:ndes) a[[k]] <- matrix(0, nrow = ncoef,
                                         ncol = ncoef)
    }
    for (j in 1:B) {
      SIMEXdata <- model.model
      epsilon =rnorm(N.group)
      epsilon = rep(epsilon* measurement.error,ni.group)
      epsilon = matrix((sqrt(lambda[i]) * epsilon ) ,ncol=1)
      SIMEXdata[, SIMEXvariable] <- SIMEXdata[, SIMEXvariable] + epsilon
      SIMEXdata[, respvar] <- SIMEXdata[, respvar] - epsilon #response is change from baseline incl. measurement error
      varb=(1+lambda[i])*measurement.error^2
      varpb=measurement.error.resp+lambda[i]*measurement.error^2
      model.SIMEX <- update(model, correlation = corCompSymm(varb/(varb+varpb),
                                                             form = as.formula(corform), fixed = T), data = data.frame(SIMEXdata))
      theta[j, ] <- model.SIMEX$coefficients$fixed
      if (jackknife.estimation[1] != FALSE) {
        variance.est <- variance.est + extract.covmat(model.SIMEX)
      }
      if (asymptotic[1]) {
        xi <- model.SIMEX$x
        psi <- psi + (resid(model.SIMEX, type = "response") *
                        xi)
        dh <- rep(1, ndes)
        if (class(model)[1] == "glm")
          dh <- model$family$mu.eta(model.SIMEX$linear.predictors)
        for (k in 1:ndes) a[[k]] <- a[[k]] - dh[k] *
          xi[k, ] %*% t(xi[k, ])
      }
    }
    estimates[i + 1, ] <- colMeans(theta)
    theta.all[[i]] <- theta
    if (jackknife.estimation[1] != FALSE) {
      variance.est <- variance.est/B
      s2 <- cov(theta)
      var.exp[[i + 1]] <- variance.est - s2
    }
    if (asymptotic[1]) {
      xiB <- psi/B
      PSI <- cbind(PSI, xiB)
      a.mat <- matrix(unlist(a), nrow = length(a), byrow = TRUE)
      ab <- matrix(colSums(a.mat), nrow = NROW(a[[1]]),
                   byrow = FALSE)
      am[[i + 1]] <- ab/(B * ndes)
    }
  }
  SIMEX.estimate <- vector(mode = "numeric", length = ncoef)
  colnames(estimates) <- p.names
  lambda <- c(0, lambda)
  switch(fitting.method, quad = extrapolation <- lm(estimates ~
                                                      lambda + I(lambda^2)), `line` = extrapolation <- lm(estimates ~
                                                                                                             lambda), nonl = extrapolation <- fitnls(lambda, p.names,
                                                                                                                                                      estimates))
  if (fitting.method[1] == "nonl") {
    for (i in 1:length(p.names)) SIMEX.estimate[i] <- predict(extrapolation[[p.names[i]]],
                                                              newdata = data.frame(lambda = -1))
  }  else {
    SIMEX.estimate <- predict(extrapolation, newdata = data.frame(lambda = -1))
  }
  if (jackknife.estimation[1] != FALSE) {
    variance.jackknife <- matrix(unlist(var.exp), ncol = ncoef^2,
                                 byrow = TRUE)
    switch(jackknife.estimation, quad = extrapolation.variance <- lm(variance.jackknife ~
                                                                       lambda + I(lambda^2)), line = extrapolation.variance <- lm(variance.jackknife ~
                                                                                                                                    lambda), nonl = extrapolation.variance <- fitnls(lambda,
                                                                                                                                                                                      1:NCOL(variance.jackknife), variance.jackknife))
    variance.jackknife2 <- vector("numeric", ncoef^2)
    switch(jackknife.estimation, nonl = for (i in 1:NCOL(variance.jackknife)) variance.jackknife2[i] <- predict(extrapolation.variance[[i]],
                                                                                                                newdata = data.frame(lambda = -1)), quad = variance.jackknife2 <- predict(extrapolation.variance,
                                                                                                                                                                                          newdata = data.frame(lambda = -1)), line = variance.jackknife2 <- predict(extrapolation.variance,
                                                                                                                                                                                                                                                                    newdata = data.frame(lambda = -1)))
    variance.jackknife <- rbind(variance.jackknife2, variance.jackknife)
    variance.jackknife.lambda <- cbind(c(-1, lambda), variance.jackknife)
    variance.jackknife <- matrix(variance.jackknife[1, ],
                                 nrow = ncoef, ncol = ncoef, byrow = TRUE)
    dimnames(variance.jackknife) <- list(p.names, p.names)
  }
  if (asymptotic[1]) {
    c11 <- cov(PSI)
    a11 <- diag.block(am)
    a11.inv <- solve(a11)
    sigma <- a11.inv %*% c11 %*% t(a11.inv)
    s <- construct.s(ncoef, lambda, fitting.method,
                             extrapolation)
    d.inv <- solve(s %*% t(s))
    sigma.gamma <- d.inv %*% s %*% sigma %*% t(s) %*% d.inv
    g <- list()
    switch(fitting.method, quad = g <- c(1, -1, 1), line = g <- c(1,
                                                                  -1), nonl = for (i in 1:ncoef) g[[i]] <- c(-1, -(coef(extrapolation[[i]])[3] -
                                                                                                                     1)^-1, coef(extrapolation[[i]])[2]/(coef(extrapolation[[i]])[3] -
                                                                                                                                                           1)^2))
    g <- diag.block(g, ncoef)
    variance.asymptotic <- (t(g) %*% sigma.gamma %*% g)/ndes
    dimnames(variance.asymptotic) <- list(p.names, p.names)
  }
  theta <- matrix(unlist(theta.all), nrow = B)
  theta.all <- list()
  for (i in 1:ncoef) theta.all[[p.names[i]]] <- data.frame(theta[,
                                                                 seq(i, ncoef * nlambda, by = ncoef)])
  z <- cbind(lambda, estimates)
  z <- rbind(c(-1, SIMEX.estimate), z)
#  colnames(z) <- c("lambda", names(coef(model))
  colnames(z) <- c("lambda", names(model$coefficients$fixed))
  erg <- list(coefficients = z[1, -1], SIMEX.estimates = z,
              lambda = lambda, model = model, measurement.error = measurement.error,
              B = B, extrapolation = extrapolation, fitting.method = fitting.method,
              SIMEXvariable = SIMEXvariable, theta = theta.all, call=cl)
  class(erg) <- ("simex")
#  fitted.values1 <- predict(erg, newdata = model.model[, -1,
#                                                      drop = FALSE], type = "response")
#  erg$fitted.values <- fitted.values
#  if (is.factor(model.model[, 1])) {
#    erg$residuals <- as.numeric(levels(model.model[, 1]))[model.model[,
#                                                                      1]] - fitted.values
#  }  else {
#    erg$residuals <- model.model[, 1] - fitted.values
#  }
  if (jackknife.estimation[1] != FALSE) {
    erg$extrapolation.variance <- extrapolation.variance
    erg$variance.jackknife <- variance.jackknife
    erg$variance.jackknife.lambda <- variance.jackknife.lambda
  }
  if (asymptotic[1]) {
    erg$PSI <- PSI
    erg$c11 <- c11
    erg$a11 <- a11
    erg$sigma <- sigma
    erg$sigma.gamma <- sigma.gamma
    erg$g <- g
    erg$s <- s
    erg$variance.asymptotic <- variance.asymptotic
  }
  return(erg)
}




