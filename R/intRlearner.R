# library(glmnet)
# library(mvtnorm)
#
# library(stringr)
# library(caret)
# library(ncvreg)
# source("rlasso.R")
# source("utils.R")
#' @include utils.R
#'
#' @title R-learner, implemented via glmnet (lasso)
#'
#' @description  R-learner, as proposed by Nie and Wager (2017), implemented via glmnet (lasso)
#'
#' @param x the input features
#' @param w the treatment variable (0 or 1)
#' @param y the observed response (real valued)
#' @param alpha tuning parameter for the elastic net
#' @param k_folds number of folds for cross-fitting
#' @param foldid user-supplied foldid. Must have length equal to length(w). If provided, it overrides the k_folds option.
#' @param lambda_y user-supplied lambda sequence for cross validation in learning E[y|x]
#' @param lambda_w user-supplied lambda sequence for cross validation in learning E[w|x]
#' @param lambda_tau user-supplied lambda sequence for cross validation in learning the treatment effect E[y(1) - y(0) | x]
#' @param lambda_choice how to cross-validate for learning the treatment effect tau; choose from "lambda.min" or "lambda.1se"
#' @param rs whether to use the RS-learner (logical).
#' @param p_hat user-supplied estimate for E[W|X]
#' @param m_hat user-supplied estimte for E[Y|X]
#' @param penalty_factor user-supplied penalty factor, a vector of length the same as the number of covariates in x.
#' @return an rlasso object
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rlasso_fit = rlasso(x, w, y)
#' rlasso_est = predict(rlasso_fit, x)
#' }
#' @export
rlasso = function(x, w, y, p,
                  alpha = 1,
                  k_folds = NULL,
                  foldid = NULL,
                  lambda_y = NULL,
                  lambda_w = NULL,
                  lambda_tau = NULL,
                  lambda_choice = c("lambda.min"),
                  rs = FALSE,
                  p_hat = NULL,
                  m_hat = NULL,
                  penalty_factor = NULL,
                  ps_const = FALSE){

  input = sanitize_input(x,w,y)
  x = input$x
  w = input$w
  y = input$y

  standardization = caret::preProcess(x, method=c("center", "scale")) # get the standardization params
  x_scl = predict(standardization, x)							 # standardize the input
  x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]

  lambda_choice = match.arg(lambda_choice)
  # print(lambda_choice)
  nobs = nrow(x_scl)
  pobs = ncol(x_scl)

  if (is.null(foldid) || length(foldid) != length(w)) {

    if (!is.null(foldid) && length(foldid) != length(w)) {
      warning("supplied foldid does not have the same length ")
    }

    if (is.null(k_folds)) {
      k_folds = floor(max(3, min(10,length(w)/4)))
    }

    # fold ID for cross-validation; balance treatment assignments
    foldid = sample(rep(seq(k_folds), length = length(w)))

  }

  # penalty factor for nuisance and tau estimators
  if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
    if (!is.null(penalty_factor) && length(penalty_factor) != pobs) {
      warning("penalty_factor supplied is not of the same length as the number of columns in x after removing NA columns. Using all ones instead.")
    }
    penalty_factor_nuisance = rep(1, pobs)
    if (rs) {
      penalty_factor_tau = c(0, rep(1, 2 * pobs))
    }
    else {
      penalty_factor_tau = c(0, rep(1, pobs))
    }
  } else {
    penalty_factor_nuisance = penalty_factor
    if (rs) {
      penalty_factor_tau = c(0, penalty_factor, penalty_factor)
    }
    else {
      penalty_factor_tau = c(0, penalty_factor)
    }
  }

  if (is.null(m_hat)){
    if (p > 2) {
      y_fit = glmnet::cv.glmnet(x, y,
                                foldid = foldid,
                                keep = TRUE,
                                lambda = lambda_y,
                                alpha = alpha,
                                penalty.factor = penalty_factor_nuisance)

      y_lambda_min = y_fit$lambda[which.min(y_fit$cvm[!is.na(colSums(y_fit$fit.preval))])]
      m_hat = y_fit$fit.preval[,!is.na(colSums(y_fit$fit.preval))][, y_fit$lambda[!is.na(colSums(y_fit$fit.preval))] == y_lambda_min]
    } else {
      y_fit <- lm(y ~ x)
      m_hat <- y_fit$fitted.values
    }
  }
  else {
    y_fit = NULL
  }

  if (is.null(p_hat)){

    if (is.logical(w)) {
      if (p > 2) {
        if (ps_const) {
          ps <- mean(w)
          p_hat <- rep(ps, length(w))
          p_hat[w==0] <- 1-ps
          w_fit <- NA
        } else {
          w_fit = glmnet::cv.glmnet(x, w,
                                    foldid = foldid,
                                    family="binomial",
                                    type.measure="deviance",
                                    keep = TRUE,
                                    lambda = lambda_w,
                                    alpha = alpha,
                                    penalty.factor = penalty_factor_nuisance)

          w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
          theta_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
          p_hat = 1/(1 + exp(-theta_hat))
        }
      } else {
        w_fit <- glm(w ~ x, family = "binomial")
        p_hat <- w_fit$fitted.value
      }
    } else {
      w_fit = glmnet::cv.glmnet(x, w,
                                foldid = foldid,
                                lambda = lambda_w,
                                keep = TRUE,
                                alpha = alpha,
                                penalty.factor = penalty_factor_nuisance)

      w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
      p_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
    }
  }
  else{
    w_fit = NULL
  }

  y_tilde = y - m_hat

  if (rs){
    x_scl_tilde = cbind(as.numeric(w - p_hat) * cbind(1, x_scl), x_scl)
    x_scl_pred = cbind(1, x_scl, x_scl * 0)
  }
  else{
    x_scl_tilde = cbind(as.numeric(w - p_hat) * cbind(1, x_scl))
    x_scl_pred = cbind(1, x_scl)
  }

  if (p > 2) {
    tau_fit = glmnet::cv.glmnet(x_scl_tilde,
                                y_tilde,
                                foldid = foldid,
                                alpha = alpha,
                                lambda = lambda_tau,
                                penalty.factor = penalty_factor_tau,
                                standardize = FALSE)
  } else {
    tau_fit <- lm(y_tilde ~ x_scl_tilde)
  }

  tau_beta = as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))

  tau_hat = x_scl_pred %*% tau_beta

  ret = list(tau_fit = tau_fit,
             tau_beta = tau_beta,
             w_fit = w_fit,
             y_fit = y_fit,
             p_hat = p_hat,
             m_hat = m_hat,
             tau_hat = tau_hat,
             rs = rs,
             standardization = standardization,
             x_scl = x_scl)
  class(ret) <- "rlasso"
  ret
}


#' predict for rlasso
#'
#' get estimated tau(x) using the trained rlasso model
#'
#' @param object an rlasso object
#' @param newx covariate matrix to make predictions on. If null, return the tau(x) predictions on the training data
#' @param ... additional arguments (currently not used)
#'
#' @examples
#' \dontrun{
#' n = 100; p = 10
#'
#' x = matrix(rnorm(n*p), n, p)
#' w = rbinom(n, 1, 0.5)
#' y = pmax(x[,1], 0) * w + x[,2] + pmin(x[,3], 0) + rnorm(n)
#'
#' rlasso_fit = rlasso(x, w, y)
#' rlasso_est = predict(rlasso_fit, x)
#' }
#'
#'
#' @return vector of predictions
#' @export
predict.rlasso <- function(object,
                           newx = NULL,
                           ...) {
  if (!is.null(newx)) {

    newx = sanitize_x(newx)
    # newx_scl <- newx
    newx_scl = predict(object$standardization, newx) # standardize the new data using the same standardization as with the training data
    newx_scl = newx_scl[,!is.na(colSums(newx_scl)), drop = FALSE]

    if (object$rs){
      newx_scl_pred = cbind(1, newx_scl, newx_scl * 0)
    }
    else{
      newx_scl_pred = cbind(1, newx_scl)
    }
    tau_hat = newx_scl_pred %*% object$tau_beta
  }
  else {
    tau_hat = object$tau_hat
  }
  return(tau_hat)
}


intrlearner <- function (x, w, y, s, p, alpha = 1, k_folds = NULL, foldid = NULL,
                         lambda_y = NULL, lambda_w = NULL, lambda_tau = NULL,
                         lambda_choice = c("lambda.min"),
                         p_hat = NULL, m_hat = NULL,
                         penalty_factor = NULL) {
  input = sanitize_input(x, w, y)
  x = input$x
  w = input$w
  y = input$y
  standardization = caret::preProcess(x, method = c("center",  "scale"))
  x_scl = predict(standardization, x)
  x_scl = x_scl[, !is.na(colSums(x_scl)), drop = FALSE]
  lambda_choice = match.arg(lambda_choice)
  nobs = nrow(x_scl)
  pobs = ncol(x_scl)
  if (is.null(foldid) || length(foldid) != length(w)) {
    if (!is.null(foldid) && length(foldid) != length(w)) {
      alpha = 1
      warning("supplied foldid does not have the same length ")
    }
    if (is.null(k_folds)) {
      k_folds = floor(max(3, min(5, length(w)/4)))  # change 10 to 5
    }
    foldid = sample(rep(seq(k_folds), length = length(w)))
  }
  if (is.null(penalty_factor) || (length(penalty_factor) != pobs)) {
    if (!is.null(penalty_factor) && length(penalty_factor) != pobs) {
      warning("penalty_factor supplied is not of the same length as the number of columns in x after removing NA columns. Using all ones instead.")
    }
    penalty_factor_nuisance = rep(1, pobs)
    penalty_factor_tau = c(0, rep(1, 2 * pobs + 1))  # intercept of tau(X) is not penalized
  } else {
    penalty_factor_nuisance = penalty_factor
    penalty_factor_tau = c(0, penalty_factor, 1, penalty_factor)
  }

  if (is.null(m_hat)) {
    m_hat <- rep(0, length(w))
    for (si in c(0, 1)) {
      if (ncol(x) > 1) {
        nonzero_col_idx <- which(apply(x[s==si, ], 2, sd) != 0)
        xs <- x[s==si, nonzero_col_idx]
        y_fit <- glmnet::cv.glmnet(xs, y[s==si],
                                   foldid = foldid[s==si],
                                   keep = TRUE,
                                   lambda = lambda_y,
                                   alpha = alpha,
                                   penalty.factor = penalty_factor_nuisance[nonzero_col_idx])
        y_lambda_min <- y_fit$lambda[which.min(y_fit$cvm[!is.na(colSums(y_fit$fit.preval))])]
        m_hat[s == si] <- y_fit$fit.preval[,!is.na(colSums(y_fit$fit.preval))][, y_fit$lambda[!is.na(colSums(y_fit$fit.preval))] == y_lambda_min]
      } else {
        df_si <- data.frame(x = x[s == si], y = y[s == si])
        foldid_si <- foldid[s == si]
        for (i in unique(foldid_si)) {
          y_fit <- lm(y ~ x, data = df_si[foldid_si != i, ])
          m_hat[s == si][foldid_si == i] <- predict(y_fit, newdata = df_si[foldid_si == i, ])
        }
        nonzero_col_idx <- NULL
      }
    }
  } else {
    y_fit = NULL
  }
  if (is.null(p_hat)) {
    p_hat <- rep(0, length(w))
    for (si in c(0, 1)) {
      if (ncol(x) > 1) {
        nonzero_col_idx <- which(apply(x[s==si, ], 2, sd) != 0)
        xs <- x[s==si, nonzero_col_idx]
        w_fit <- glmnet::cv.glmnet(xs, w[s==si],
                                   foldid = foldid[s==si],
                                   family="binomial",
                                   type.measure="deviance",
                                   keep = TRUE,
                                   lambda = lambda_w,
                                   alpha = alpha,
                                   penalty.factor = penalty_factor_nuisance[nonzero_col_idx])
        w_lambda_min = w_fit$lambda[which.min(w_fit$cvm[!is.na(colSums(w_fit$fit.preval))])]
        theta_hat = w_fit$fit.preval[,!is.na(colSums(w_fit$fit.preval))][, w_fit$lambda[!is.na(colSums(w_fit$fit.preval))] == w_lambda_min]
        p_hat[s == si] = 1/(1 + exp(-theta_hat))
      } else {
        # One-dimensional covariates, use linear regression rather than lasso.
        df_si <- data.frame(x = x[s == si], y = w[s == si])
        foldid_si <- foldid[s == si]
        for (i in unique(foldid_si)) {
          w_fit <- glm(y ~ x, data = df_si[foldid_si != i, ], family = "binomial")
          p_hat[s == si][foldid_si == i] <- predict(w_fit, newdata = df_si[foldid_si == i, ], type = "response")
        }
        nonzero_col_idx <- NULL
      }
    }
  } else {
    w_fit <- NULL
  }

  y_tilde <- y - m_hat
  eps_w <- as.numeric(w - p_hat)
  x_scl_tilde <- cbind(eps_w * cbind(1, x_scl), (1-s) * eps_w * cbind(1, x_scl))
  x_scl_pred <- cbind(1, x_scl, x_scl * 0)

  temp_penalty_factor_tau <- penalty_factor_tau[1:(1+pobs)]
  # No standardization
  x_tilde <-  cbind(eps_w * cbind(1, x), (1-s) * eps_w * cbind(1, x))

  if (is.null(penalty_factor)) {
    print("Parameter searching ...")
    temp_cvm <- Inf
    for (sc in seq(0, 1.5, 0.5)) {
      penalty_factor_tau[1:(1+pobs)] <- sc * temp_penalty_factor_tau
      for (ga in seq(3, 8, 0.5)) {
        temp_tau_fit_scad = ncvreg::cv.ncvreg(x_scl_tilde, y_tilde, fold = foldid,
                                              penalty = "SCAD",
                                              gamma = ga,
                                              penalty.factor = penalty_factor_tau)#,

        cve <- min(temp_tau_fit_scad$cve)
        if (cve < temp_cvm) {
          temp_cvm <- cve
          tau_fit <- temp_tau_fit_scad
          sc_chosen <- sc
          ga_chosen <- ga
        }
      }
    }
    cat("penalty scale parameter is chosen as", sc_chosen, "\n")
  } else {
    tau_fit = glmnet::cv.glmnet(x_scl_tilde, y_tilde, foldid = foldid,
                                alpha = alpha, lambda = lambda_tau,
                                penalty.factor = penalty_factor_tau,
                                standardize = FALSE)
  }
  coef_tau_fit <- coef(tau_fit)
  # Delete the intercept of the regression (index = 1) and
  #   the intercept in the confounding function (index = pobs + 3)
  tau_beta = as.vector(t(coef_tau_fit[-c(1, (pobs+3))]))
  # Index of the selected variables in the tau(x)
  idx_nonzero_tau <- which(tau_beta[1:(pobs + 1)] != 0)

  tau_hat = x_scl_pred %*% tau_beta
  ret = list(tau_fit = tau_fit, tau_beta = tau_beta,
             w_fit = w_fit, y_fit = y_fit, rs = TRUE, # rs is set to be TRUE in order to be used in predict.rlasso in the origianl Nie's rlasso function
             p_hat = p_hat, m_hat = m_hat, tau_hat = tau_hat,
             standardization = standardization, x_scl = x_scl,
             sc_chosen = sc_chosen,  ga_chosen = ga_chosen,
             nonzero_col_idx = nonzero_col_idx)
  class(ret) <- "rlasso"
  ret
}


poly_interaction <- function(x, d=2, interact=TRUE) {
  # To expand the feature matrix with polynomial terms of each column of x and interactions between each other.
  # x: n by p matrix.
  # d: the polynomial degree >= 2.
  # interact: a boolean type variable to decide whether to add interaction terms.
  x <- x_poly <- as.matrix(x)
  p <- ncol(x)
  if (p == 1) interact = FALSE  # No interactions naturally when p = 1.
  if (d > 1) {
    for (dd in 2:d) {
      poly_d <- x^dd
      x_poly <- cbind(x_poly, poly_d)
      if (interact) {
        x_interact <- x_poly[, 1:(ncol(x_poly)-p)]
        combn2 <- combn(1:ncol(x_interact), 2)
        interactions <- x_interact[, combn2[1, ]] * x_interact[, combn2[2, ]]
        x_poly  <- cbind(x_poly, interactions)
      }
    }
  }
  x_poly
}

# Same data generating setting as KPS's paper
gen_y <- function(x, u, a, sample_size, cf) {
  # When coef = 0: tau(X) = 1 + 2X; o.w., the confounding function c(X) is more restrictive than tau(X).
  eps <- rnorm(sample_size)
  1 + a + x + 2*a*x + 0.5*x^2 + cf*a*x^2 + u + 0.5*eps
  # 1 + x + 0.5*x^2 + u + 0.5*eps  # tau(x) = 0
}

gen_data <- function(m, N, b=NULL, p=10, lambda="lambda.min", n0=NULL, n=NULL,
                     KPS_data=FALSE, multidim=NA, cf=0, p0=1) {
  # m, N: sample size of RWE, and population data, respectively
  # b: indicates the different strengths of unmeasured confounding
  # p: the dimension of covariates x
  # d: the order of the polynomial basis regression
  # lambda: a string variable used in cv.glmnet, {"lambda.min", "lambda1se"}
  # n0: the intercept in the logit{P(A=1|RCT)} when KPS_data=FALSE
  # n: sample size of RCT when KPS_data=TRUE
  # cf \in {0, 0.75}, no use when KPS_data=FALSE
  # p0: the number of non-zero coefs in tau(X), an integer \in {1, p-1}
  if (KPS_data) {
    # Data generated similar as the experiments in NK's paper
    x_rct <- matrix(runif(n, min = -1, max = 1), ncol = 1)
    u_rct <- matrix(rnorm(n), ncol = 1)
    a_rct <- matrix(rbinom(n, 1, 0.5), ncol = 1)
    y_rct <- gen_y(x_rct, u_rct, a_rct, n, cf)
    a_rwe <- matrix(rbinom(m, 1, 0.5), ncol = 1)
    xu_rwe <- matrix(NA, m, 2)
    for (i in 1:m) {
      xu_rwe[i, ] <- rmvnorm(n=1, mean=c(0, 0), sigma=matrix(c(1, a_rwe[i]-0.5, a_rwe[i]-0.5, 1), 2, 2))
    }
    x_rwe <- matrix(xu_rwe[, 1], ncol=1)
    u_rwe <- matrix(xu_rwe[, 2], ncol=1)
    y_rwe <- gen_y(x_rwe, u_rwe, a_rwe, m, cf)
    x <- matrix(runif(N, min = -1, max = 1), ncol=1)  # test data should have the same support as rct
    tau <- 0
    true_tau_coef <- c(1, 2, cf)
    if (multidim) {
      x_rct <- cbind(x_rct, matrix(runif(n*(p-2), min=-1, max=1), ncol=p-2))
      x_rwe <- cbind(x_rwe, matrix(rnorm(m*(p-2)), ncol=p-2))
      x <- cbind(x, matrix(rnorm(N*(p-2)), ncol=p-2))
      true_tau_coef <- c(rep(1, p-2-p0), rep(0, p0))
      tau <- c(x[, 2:(p-1)] %*% true_tau_coef) + tau
    }
    x_test <- x
  } else {
    x <-  rmvnorm(n=N, mean=rep(0, p))
    # Xp is the unmeasured confounder
    mu_x <- x[, 1] + x[, p]
    true_tau_coef <- c(rep(1, p-p0), rep(0, p0))
    tau <- c(cbind(1, x[, 1:(p-1)]) %*% true_tau_coef)
    true_nonzero_ind <- which(true_tau_coef[-1] != 0)
    y1_star <-  mu_x + tau + rnorm(N)
    y0_star <-  mu_x + rnorm(N)
    # Generate RCT
    exp_term <- exp(-n0 - x[, 1] + x[, 2])
    delta_prob <- exp_term / (1 + exp_term)
    delta1 <- rbinom(N, 1, delta_prob)
    (n <- sum(delta1))
    delta1_ind <- which(delta1 == 1)
    x_rct <- x[delta1_ind, 1:(p-1)]
    exp_term <- exp(-1 - x_rct[, 1] + x_rct[, 2])
    a_rct_prob <-  exp_term / (1 + exp_term)
    a_rct <- rbinom(n, 1, a_rct_prob)
    y_rct <- y1_star[delta1_ind] * a_rct + y0_star[delta1_ind] * (1 - a_rct)
    # Generate X in RWE
    delta0_ind <- sample(N, size=m)
    x_rwe <- x[delta0_ind, 1:(p-1)]
    x_test <- x[, 1:(p-1)]
  }
  # Polynomial basis
  x_rct <- poly_interaction(x_rct, interact = T)
  x_rwe <- poly_interaction(x_rwe, interact = T)
  x_test <- poly_interaction(x_test, interact = T)
  newx_rlearner <- matrix(x_test, ncol=ncol(x_test))
  if (!KPS_data) {
    # Generate A and Y in RWE data
    exp_term <- exp(-1 - x[, 1] + x[, 2] + b * x[, p])
    a_rwe_prob <-  exp_term / (1 + exp_term)
    a <- rbinom(N, 1, a_rwe_prob)
    y <- y1_star * a + y0_star * (1 - a)
    a_rwe <- a[delta0_ind]
    y_rwe <- y[delta0_ind]
  }

  list(x_rct = x_rct, x_rwe = x_rwe,
       a_rct = a_rct, a_rwe = a_rwe,
       y_rct = y_rct, y_rwe = y_rwe,
       newx_rlearner = newx_rlearner,
       tau = tau)
}
