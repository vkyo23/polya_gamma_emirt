library(Rcpp)
library(RcppArmadillo)

# Loading Cpp functions =======================
cpp_dir <- "../src/"
cpp_fn <- list.files(cpp_dir)
for (c in cpp_fn) {
  sourceCpp(paste0(cpp_dir, c))
}

# Loading wrapper function ====================
pg_emirt <- function(Y, init = NULL, prior = NULL, constraint = NULL, maxit = 500, tol = 1e-6, verbose = 10) {
  
  cat("===========================================================================\n")
  cat("2PL Item Response Model with EM Algorithm and Polya-Gamma Data Augmentation\n")
  cat("===========================================================================\n")
  
  # Measuring implementation times
  stime <- proc.time()[3]
  
  if (is.null(init)) {
    cat("Initial values is not supplied, so computing from the data.....")
    theta_init <- as.numeric(scale(rowMeans(Y, na.rm = TRUE)))
    alpha_init <- beta_init <- rep()
    for (j in 1:ncol(Y)) {
      c <- suppressWarnings(coef(glm(Y[, j] ~ theta_init, family = "binomial")))
      alpha_init[j] <- ifelse(abs(c[1]) > 10, 0.1, c[1])
      beta_init[j] <- ifelse(abs(c[2]) > 10, 0.1, c[2])
    }
    init <- list(theta = theta_init, alpha = alpha_init, beta = beta_init)
    cat("DONE!\n")
  }
  if (is.null(prior)) {
    cat("Prior is not supplied, so automatically setting:\n")
    cat("alpha ~ N(0, 100), beta ~ N(0, 100)\n")
    prior <- list(a0 = 0, A0 = 100, b0 = 0, B0 = 100)
  }
  
  # Get initial values
  theta <- init$theta
  alpha <- init$alpha
  beta <- init$beta
  
  # If theta constraint is not supplied
  if (is.null(constraint)) constraint <- which.max(theta)
  if (theta[constraint] < 0) theta <- -theta
  
  # Get priors
  a0 <- prior$a0
  A0 <- prior$A0
  b0 <- prior$b0
  B0 <- prior$B0
  
  # Convergence check statistics store (correlation between current and previous parameters)
  all_cor <- matrix(NA, maxit, 3)
  colnames(all_cor) <- c("alpha", "beta", "theta")
  
  # EM
  cat("Start EM:\n") 
  for (g in 1:maxit) {
    ## Previous parameters
    theta_old <- theta
    alpha_old <- alpha
    beta_old <- beta
    
    ## E-step (calculate E[omega])
    Omega <- draw_E_Omega(alpha_old, beta_old, theta_old)
    
    ## M-step (theta -> alpha -> beta)
    theta <- draw_theta(Y, alpha_old, beta_old, Omega, constraint)
    alpha <- draw_alpha(Y, beta_old, theta, Omega, a0, A0)
    beta <- draw_beta(Y, alpha, theta, Omega, b0, B0)
    
    ## Convergence stat
    all_cor[g, 1] <- cor(alpha_old, alpha)
    all_cor[g, 2] <- cor(beta_old, beta)
    all_cor[g, 3] <- cor(theta_old, theta)
    
    ## If converged
    if (min(1 - all_cor[g, ]) < tol) {
      ### Conv stat
      all_cor <- na.omit(all_cor)
      attr(all_cor, "na.action") <- NULL
      
      ### Parameters
      names(theta) <- rownames(Y)
      parameters <- list(alpha = alpha, beta = beta, theta = theta)
      
      ### Organize
      L <- list(parameters = parameters,
                converge = TRUE,
                conv_check = all_cor,
                iter = g,
                data = list(Y = Y,
                            init = init,
                            prior = prior,
                            constraint = constraint,
                            maxit = maxit,
                            tol = tol))
      
      ### Record time
      etime <- proc.time()[3]
      cat("Model converged at iteration", g, ": total time", round(etime - stime, 3), "sec\n")
      return(L)
    }
    
    ## If failed to converge
    if (g == maxit) {
      ### Parameters
      names(theta) <- rownames(Y)
      parameters <- list(alpha = alpha, beta = beta, theta = theta)
      
      ### Organize
      L <- list(parameters = parameters,
                converge = FALSE,
                conv_check = all_cor,
                iter = g,
                data = list(Y = Y,
                            init = init,
                            prior = prior,
                            constraint = constraint,
                            maxit = maxit,
                            tol = tol))
      
      ### Record time
      etime <- proc.time()[3]
      cat("Model failed to converge", ": total time", round(etime - stime, 3), "sec\n")
      warning("The result may not be reliable; return the result as it is.")
      return(L)
    }
    
    ### Print the status
    if (g %% verbose == 0) {
      mtime <- proc.time()[3]
      cat("Iteration", g, "eval =", min(1 - all_cor[g, ]), ": time", round(mtime - stime, 3), "sec\n")
    }
  }
} 

# Loading parametric bootstrap function ===================================
pg_em_boot <- function(model, boot = 100, verbose = 10) {
  
  cat("===========================================\n")
  cat("Parametric bootstrap for Polya-Gamma EM IRT\n")
  cat("===========================================\n")
  
  stime <- proc.time()[3]
  
  # For suppressing the messages
  quiet <- function(x) { 
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  } 
  
  # Extract the estimated parameters
  alpha <- model$parameters$alpha
  beta <- model$parameters$beta
  theta <- model$parameters$theta
  
  # Storage
  theta_boot <- matrix(NA, length(theta), boot)
  rownames(theta_boot) <- names(theta)
  alpha_boot <- beta_boot <- matrix(NA, length(alpha), boot)
  
  # Bootstrap
  for (b in 1:boot) {
    # Generating random binomial variates (prediction of Y)
    set.seed(b)
    p <- plogis(cbind(1, theta) %*% rbind(alpha, beta))
    y_pred <- matrix(rbinom(length(theta) * length(alpha), 1, p),
                     length(theta), length(alpha))
    
    # Estimation
    fit_boot <- quiet(pg_emirt(Y = y_pred,
                               init = model$data$init,
                               prior = model$data$prior,
                               constraint = model$data$constraint,
                               maxit = model$data$maxit,
                               tol = model$data$tol,
                               verbose = model$data$maxit + 1)) 
    
    # Saving
    theta_boot[, b] <- fit_boot$parameters$theta
    alpha_boot[, b] <- fit_boot$parameters$alpha
    beta_boot[, b] <- fit$parameters$beta
    
    if (b %% verbose == 0) {
      mtime <- proc.time()[3]
      cat("Bootstrap", b, "done : time", round(mtime - stime, 3), "sec\n")
    }
  }
  
  # Return the results
  etime <- proc.time()[3]
  L <- list(alpha = alpha_boot, beta = beta_boot, theta = theta_boot)
  cat("Bootstrap finished : total time", round(etime - stime, 3), "sec\n")
  return(L)
}
