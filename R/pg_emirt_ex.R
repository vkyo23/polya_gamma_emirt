# Load packages ============================================
library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
setwd("/mnt/c/Users/Kanetaka/Dropbox/R/polya_gamma_emirt/R")

# Load functions ===========================================
source("pg_emirt_fn.R")

# Load data ================================================
data("Senate", package = "MCMCpack")
Y <- as.matrix(Senate[, 6:ncol(Senate)])

# Generate initial values, constraint and priors ===========
# Calculate the standardized value of the probability of yea response 
# as initial value of theta
theta_init <- as.numeric(scale(rowMeans(Y, na.rm = TRUE)))

# With initial value of theta, get initial values of alpha and beta
# by regressing Y on theta
alpha_init <- beta_init <- rep()
for (j in 1:ncol(Y)) {
  c <- coef(glm(Y[, j] ~ theta_init, family = "binomial"))
  alpha_init[j] <- ifelse(abs(c[1]) > 10, 0.1, c[1])
  beta_init[j] <- ifelse(abs(c[2]) > 10, 0.1, c[2])
}

# Aggregate
init <- list(theta = theta_init,
             alpha = alpha_init,
             beta = beta_init)

# Set constraint such that theta of HELMS always takes a positive value
constraint <- which(rownames(Y) == "HELMS")
if (theta_init[constraint] < 0) theta_init <- -theta_init

# Priors
prior <- list(a0 = 0, A0 = 25, b0 = 0, B0 = 25)


# Run PG-EM IRT =============================================
model <- pg_emirt(Y = Y, init = init, prior = prior, constraint = constraint,
                  maxit = 500, tol = 1e-6, verbose = 10)

# Boostrap ==================================================
boot <- pg_em_boot(model)

# Result ====================================================
# Bias correction
## Calculate bootstrapm mean
boot_mean <- rowMeans(boot$theta)

## Bias (estimated mean - boot mean)
bias <- fit$parameters$theta - boot_mean

## 95% CI
lwr <- apply(boot$theta, 1, quantile, probs = .025) + bias
upr <- apply(boot$theta, 1, quantile, probs = .975) + bias

## Plot
pdf("../Output/result.pdf")
tibble(name = Senate$member, 
       party = Senate$party, 
       theta = fit$parameters$theta,
       lwr = lwr,
       upr = upr) %>% 
  mutate(party = if_else(party == 1, "Republican", "Democrat") %>% 
           factor(levels = c("Republican", "Democrat"))) %>% 
  ggplot(aes(x = theta, y = reorder(name, theta), 
             color = party)) +
  geom_pointrange(aes(xmin = lwr, xmax = upr), size = 0.3) +
  xlab("Ideal Point") + 
  theme_light() +
  theme(legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 5)) 
dev.off()
