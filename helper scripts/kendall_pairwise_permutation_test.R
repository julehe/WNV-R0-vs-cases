# ==============================================================================
# Author: Julian Heidecke        
# Email: julian.heidecke@iwr.uni-heidelberg.de or julian.heidecke@gmail.com
#
# Description: This provides a helper function to perform pairwise permutation 
# test to test statistical significance of differences of rank correlation coeffs.
#
# ==============================================================================

library(pcaPP)

kendall_pairwise_permutation_test <- function(df,
                                              target_col,
                                              predictor_cols,
                                              nsim = 1000,
                                              seed = 101,
                                              two_sided = FALSE) {
  set.seed(seed)
  # Extract response variable and predictor subset
  A <- df[[target_col]]
  pred_df <- df[, predictor_cols, drop = FALSE]
  predictor_names <- colnames(pred_df)
  
  comparison_results <- list()
  
  n <- length(A)
  
  # All unique pairs
  combn(seq_along(predictor_names), 2, simplify = FALSE) %>%
    walk(function(pair_idx) {
      i <- pair_idx[1]
      j <- pair_idx[2]
      
      p1 <- predictor_names[i]
      p2 <- predictor_names[j]
      
      B <- pred_df[[i]]
      C <- pred_df[[j]]
      
      B_rank <- rank(B)
      C_rank <- rank(C)
      
      tau_B <- cor.fk(A, B_rank)
      tau_C <- cor.fk(A, C_rank) # ---> cor.fk is much faster than cor.test
      
      d_obs <- tau_C - tau_B
      
      res <- numeric(nsim)
      for (i in 1:nsim) {
        mask <- sample.int(2, n, replace = TRUE) == 1L
        permuted_ranksB <- ifelse(mask, B_rank, C_rank)
        permuted_ranksC <- ifelse(mask, C_rank, B_rank)
        
        tau_permB <- cor.fk(A, permuted_ranksB)
        tau_permC <- cor.fk(A, permuted_ranksC)
        
        res[i] <- tau_permC - tau_permB
        
        print(i)
      }
      
      p_val <- if (two_sided) {
        mean(abs(res) >= abs(d_obs))
      } else {
        mean(res >= d_obs)
      }
      
      comparison_results[[paste(p1, p2, sep = "_vs_")]] <<- tibble(
        predictor1 = p1,
        predictor2 = p2,
        tau_B = tau_B,
        tau_C = tau_C,
        observed_diff = d_obs,
        p_value = p_val
      )
    })
  
  list(
    comparisons = bind_rows(comparison_results)
  )
}
