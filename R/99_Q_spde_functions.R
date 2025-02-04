# Q_spde functions

# Copy of sdmTMB's MakeH() in utils.h
# https://github.com/pbs-assess/sdmTMB/blob/main/src/utils.h#L232C1-L240C2
MakeH <- function(x) {
  H <- matrix(0, nrow = 2, ncol = 2)
  H[1, 1] <- exp(x[1])           
  H[2, 1] <- x[2]                
  H[1, 2] <- x[2]               
  H[2, 2] <- (1 + x[2]^2) / exp(x[1])  
  return(H)
}

# Copy of TMB's functions in R_inla.hpp
# https://kaskr.github.io/adcomp/R__inla_8hpp_source.html
Q_spde <- function(spde, kappa) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  kappa_pow4 * spde$M0 + 2 * kappa_pow2 * spde$M1 + spde$M2
}


Q_spde_aniso <- function(spde, kappa, H) {
  kappa_pow2 <- kappa * kappa
  kappa_pow4 <- kappa_pow2 * kappa_pow2
  
  n_s <- spde$n_s
  n_tri <- spde$n_tri
  Tri_Area <- spde$Tri_Area
  E0 <- spde$E0
  E1 <- spde$E1
  E2 <- spde$E2
  TV <- spde$TV
  G0 <- spde$G0
  G0_inv <- AD(spde$G0_inv)
  
  # G1_aniso <- matrix(0, nrow = n_s, ncol = n_s)
  # G2_aniso <-  matrix(0, nrow = n_s, ncol = n_s)
  G1_aniso <- AD(Matrix::Matrix(0, nrow = n_s, ncol = n_s))
  #G2_aniso <-  AD(Matrix::Matrix(0, nrow = n_s, ncol = n_s))
  
  # Calculate adjugate of H
  adj_H <- matrix(0, nrow = 2, ncol = 2)
  adj_H[1, 1] <- H[2, 2]
  adj_H[1, 2] <- -H[1, 2]
  adj_H[2, 1] <- -H[2, 1]
  adj_H[2, 2] <- H[1, 1]
  
  # Calculate G1 - pt. 1
  Gtmp <- AD(array(0, dim = c(n_tri, 3, 3)))
  for (i in 1:n_tri) {
    Gtmp[i, 1, 1] <- (E0[i, 1] * (E0[i, 1] * adj_H[1, 1] + E0[i, 2] * adj_H[2, 1]) +
                        E0[i, 2] * (E0[i, 1] * adj_H[1, 2] + E0[i, 2] * adj_H[2, 2])) / (4 * Tri_Area[i])
    Gtmp[i, 1, 2] <- (E1[i, 1] * (E0[i, 1] * adj_H[1, 1] + E0[i, 2] * adj_H[2, 1]) +
                        E1[i, 2] * (E0[i, 1] * adj_H[1, 2] + E0[i, 2] * adj_H[2, 2])) / (4 * Tri_Area[i])
    Gtmp[i, 1, 3] <- (E2[i, 1] * (E0[i, 1] * adj_H[1, 1] + E0[i, 2] * adj_H[2, 1]) +
                        E2[i, 2] * (E0[i, 1] * adj_H[1, 2] + E0[i, 2] * adj_H[2, 2])) / (4 * Tri_Area[i])
    Gtmp[i, 2, 2] <- (E1[i, 1] * (E1[i, 1] * adj_H[1, 1] + E1[i, 2] * adj_H[2, 1]) +
                        E1[i, 2] * (E1[i, 1] * adj_H[1, 2] + E1[i, 2] * adj_H[2, 2])) / (4 * Tri_Area[i])
    Gtmp[i, 2, 3] <- (E2[i, 1] * (E1[i, 1] * adj_H[1, 1] + E1[i, 2] * adj_H[2, 1]) +
                        E2[i, 2] * (E1[i, 1] * adj_H[1, 2] + E1[i, 2] * adj_H[2, 2])) / (4 * Tri_Area[i])
    Gtmp[i, 3, 3] <- (E2[i, 1] * (E2[i, 1] * adj_H[1, 1] + E2[i, 2] * adj_H[2, 1]) +
                        E2[i, 2] * (E2[i, 1] * adj_H[1, 2] + E2[i, 2] * adj_H[2, 2])) / (4 * Tri_Area[i])
  }
  
  # Calculate G1 - pt. 2
  for (i in 1:n_tri) {
    #  if (i == 18) browser()
    G1_aniso[TV[i, 2]  + 1 , TV[i, 1]  + 1 ] <- G1_aniso[TV[i, 2] + 1 , TV[i, 1] + 1 ] + Gtmp[i, 1, 2]
    G1_aniso[TV[i, 1]  + 1 , TV[i, 2]  + 1 ] <- G1_aniso[TV[i, 1] + 1 , TV[i, 2] + 1 ] + Gtmp[i, 1, 2]
    G1_aniso[TV[i, 3]  + 1 , TV[i, 2]  + 1 ] <- G1_aniso[TV[i, 3] + 1 , TV[i, 2] + 1 ] + Gtmp[i, 2, 3]
    G1_aniso[TV[i, 2]  + 1 , TV[i, 3]  + 1 ] <- G1_aniso[TV[i, 2] + 1 , TV[i, 3] + 1 ] + Gtmp[i, 2, 3]
    G1_aniso[TV[i, 3]  + 1 , TV[i, 1]  + 1 ] <- G1_aniso[TV[i, 3] + 1 , TV[i, 1] + 1 ] + Gtmp[i, 1, 3]
    G1_aniso[TV[i, 1]  + 1 , TV[i, 3]  + 1 ] <- G1_aniso[TV[i, 1] + 1 , TV[i, 3] + 1 ] + Gtmp[i, 1, 3]
    G1_aniso[TV[i, 1]  + 1 , TV[i, 1]  + 1 ] <- G1_aniso[TV[i, 1] + 1 , TV[i, 1] + 1 ] + Gtmp[i, 1, 1]
    G1_aniso[TV[i, 2]  + 1 , TV[i, 2]  + 1 ] <- G1_aniso[TV[i, 2] + 1 , TV[i, 2] + 1 ] + Gtmp[i, 2, 2]
    G1_aniso[TV[i, 3]  + 1 , TV[i, 3]  + 1 ] <- G1_aniso[TV[i, 3] + 1 , TV[i, 3] + 1 ] + Gtmp[i, 3, 3]
  }
  #browser()
  G2_aniso <- AD(G1_aniso %*% G0_inv %*% G1_aniso)
  
  kappa_pow4 * G0 + 2 * kappa_pow2 * G1_aniso + G2_aniso
}
