# `polya_gamma_emirt`: Item Response Theory model with EM algorithm and Polya-Gamma data augmentation

Author: 

- Ukyo Kanetaka

## Description

In this repository, I provide `R` function for the implementation of 2PL Item Response Theory (IRT) model with EM Algrotihm and Polya-Gamma data augmentation (Polson *el al*. 2013). The algorithm here is based on the procedure proposed by Goplerud (2019). To estimate statistical uncertainty, I use the parametric bootstrap method proposed by Lewis and Poole (2004). For the fast estimation of the model, I use auxilary `C++` functions.

For those who are interested in the application, I provide an example at https://sunaninattahito.hatenablog.com/entry/2022/04/29/212901 (but only Japanese. I will write it in English soon). 

## R codes

R codes below are in [R](https://github.com/vkyo23/polya_gamma_emirt/tree/main/R) directory.

- `pg_emirt_ex.R`: R code example for using the function. I use roll call data for all Senators in the 106th Senate (`MCMCpack::Senate`).
- `pg_emirt_fn.R`: R code contatining 2 function.
  - `pg_emirt`: Function for implementation of Polya-Gamma EM IRT. For stopping rule, I employ the correlation between current and previous parameters. If the correation is much high (1 - cor < tolerance), the algorithm is diagnosed as converged.
  - `pg_em_boot`: Function for parametric bootstrap of the model.

## C++ codes

C++ codes below are in [src](https://github.com/vkyo23/polya_gamma_emirt/tree/main/src) directory.

- `draw_E_Omega.cpp`: Function for computing ω.
- `draw_alpha.cpp`: Function for update α (difficulty parameter or cut-point).
- `draw_beta.cpp`: Function for update β (discrimination parameter).
- `draw_theta.cpp`: Function for update Θ (ability or ideal point).

## Reference

- Goplerud, M. (2019). "A Multinomial Framework for Ideal Point Estimation". *Political Analysis*, 27(1), 69-89.
- Lewis, J. B., & Poole, K. T. (2004). "Measuring bias and uncertainty in ideal point estimates via the parametric bootstrap". *Political Analysis*, 12(2), 105-127.
- Polson, N. G., Scott, J. G., & Windle, J. (2013). "Bayesian inference for logistic models using Pólya–Gamma latent variables". *Journal of the American statistical Association*, 108(504), 1339-1349.
