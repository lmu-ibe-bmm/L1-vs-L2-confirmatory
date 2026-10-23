###############################################################
# VISUALIZATION OF CORRELATION SUMMARIES FOR REAL-DATA ANALYSIS
###############################################################
# PURPOSE
# Check whether determinant or condition number of correlation matrix as 
# possible global collinearity measures represent the intended correlation 
# levels from the simulation study well.
#
#
# WHY THIs MATTERS
# For comparison of the results of the real-data analysis and the main 
# simulation, the same global correlation/collinearity measure has to be applied
# for all real datasets and the theoretical correlation matrix of the main 
# simulation. To ensure a meaningful comparison, the chosen measure should 
# reflect the intended correlation levels from the simulation study well.
#
# WHEN TO RUN
# This check is independent of all other pilots.
###############################################################
# PURPOSE
# Construct correlation matrix.
#
# WHY THIS FUNCTION EXISTS
# Need to construct correlation matrices for each combination of rho and p.
#
# INPUTS
# rho: pairwise correlation within each 10 x 10 block
# p: dimension of the matrix 
#
# RETURN VALUE
# pxp block-diagonal correlation matrix, where within each 10x10 block, pairwise
# correlations equal rho, and variables from different blocks are uncorrelated.
constructR <- function(rho, p) {
  stopifnot(p %% 10 == 0)
  G <- p / 10
  return(kronecker(diag(G), (matrix(rho, 10, 10) + diag(1 - rho, 10))))
}

# PURPOSE
# Calculate overall measures of correlation.
#
# WHY THIS FUNCTION EXISTS
# Calculate determinant and condition number of correlation matrices for each 
# combination of rho and p.
#
# INPUTS
# rho: pairwise correlation within each 10 x 10 block
# p: dimension of the matrix 
#
# RETURN VALUE
# Named numeric vector: 
#   d: determinant of block-diagonal correlation matrix.
#   k: condition number (kappa) of block-diagonal correlation matrix.
calcCorSum <- function(rho, p) {
  R <- constructR(rho, p)
  d <- det(R)
  k <- kappa(R)
  return(c(d = d, k = k))
}
###############################################################

# Calculate determinant and condition number for grid of possible combinations 
# of rho and p
vals <- expand.grid(rho = seq(0, 0.95, by = 0.05), p = c(40, 200, 1000, 5000))
res <- do.call(mapply, c(list(FUN = calcCorSum), vals))
res.cor <- data.frame(t(res), rho = vals$rho, p = vals$p)

# Plot determinant vs. rho for each p
pdf(file.path("outputs", "cor_summaries", "det_cor_mat.pdf"), height = 4, width = 6)
op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), las = 1, mar = c(4, 4, 2.5, 1), mgp = c(2, 0.75, 0))
lapply(split(res.cor, res.cor$p),
       function(x) {
         plot(x$rho, x$d, type = "b", pch = 16,
                        main = paste0("p = ", x$p[1]),
                        xlab = bquote(rho), ylab = expression(det(Sigma)))
       })
par(op)
dev.off()

# Plot condition number vs. rho for each p
pdf(file.path("outputs", "cor_summaries", "cond_cor_mat.pdf"), height = 4, width = 6)
op <- par(no.readonly = TRUE)
par(mfrow = c(2, 2), las = 1, mar = c(4, 4, 2.5, 1), mgp = c(2, 0.75, 0))
lapply(split(res.cor, res.cor$p), 
       function(x) {
         plot(x$rho, x$k, type = "b", pch = 16, 
              main = bquote(p == .(x$p[1])), 
              xlab = bquote(rho), ylab = expression(kappa(Sigma)), 
              panel.first = abline(h = c(10, 30), col = "darkgrey", lty = "dotted"))
       })
par(op)
dev.off()
