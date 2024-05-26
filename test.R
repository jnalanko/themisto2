log.bb.comp <- function(x, n, a, b) {
    lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(n + a, b)
}

log.norm.const <- function(x, n, a, b) {
    sum(log(a + n + x - 1:n) - log(a + b + 2*n - 1:n))
}

log.ll <- function(x, n, a, b) {
    log.bb.comp(x, n, a, b) + log.norm.const(x, n, a, b)
}

log.ll2 <- function(x, n, a, b) {
    lchoose(n, x) + lbeta(x + a, n - x + b) - lbeta(a, b)
}

c.sizes <- c(10, 20)
n.max <- c.sizes
e <- c.sizes*0.65
phi <- 1.0/(c.sizes - e + 0.01)
beta <- phi*(c.sizes - e)
alpha <- (e*beta)/(c.sizes - e)
joint <- matrix(0, nrow = n.max[1] + 1, ncol = n.max[2] + 1)
single <- list(numeric(n.max[1] + 1), numeric(n.max[2] + 1))
joint2 <- matrix(0, nrow = n.max[1] + 1, ncol = n.max[2] + 1)
single2 <- list(numeric(n.max[1] + 1), numeric(n.max[2] + 1))
for (i in 0:n.max[1]) {
    for (j in 0:n.max[2]) {
        single[[1]][i + 1] <- log.ll(i, n.max[1], alpha[1], beta[1])
        single[[2]][j + 1] <- log.ll(j, n.max[2], alpha[2], beta[2])
        single2[[1]][i + 1] <- log.ll2(i, n.max[1], alpha[1], beta[1])
        single2[[2]][j + 1] <- log.ll2(j, n.max[2], alpha[2], beta[2])
        joint[i + 1, j + 1] <- exp(single[[1]][i + 1] + single[[2]][j + 1])
        joint2[i + 1, j + 1] <- exp(single2[[1]][i + 1] + single2[[2]][j + 1])
    }
}

library("RColorBrewer")

gradient <- colorRampPalette(rev(c("#ca0020", "#f4a582", "white", "#92c5de", "#0571b0")))


##pdf(file = "img/mSWEEP_likelihood.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 0))
plot(0:10, single[[1]], type = 'p', col = "#ca0020", pch = 19, ylim = c(-30, 0), ylab = "Log-likelihood", xlab = "Number of successess", xaxt = 'n', yaxt = 'n', bty = 'n')
axis(side = 1, at = seq(0, 10, by = 2))
axis(side = 2, at = seq(0, -30, by = -5))
lines(0:10, single[[1]], lty = "dashed", col = "#ca0020")
lines(0:10, single2[[1]], type = 'p', col = "#0571b0", pch = 19)
lines(0:10, single2[[1]], lty = "dashed", col = "#0571b0")
title(main = "a) Log-likelihood for n = 10", font.main = 1, adj = 0.0)
legend("bottomright", col = c("#ca0020", "#0571b0"), fill = c(NA, NA), pch = c(19, 19), lty = c("dashed", "dashed"), legend = c("mSWEEP model", "Beta-binomial model"), bty = 'n', border = NA)
par(mar = c(4, 3, 1, 0))
plot(0:20, single[[2]], type = 'p', col = "#ca0020", pch = 19, ylim = c(-30, 0), ylab = "", xlab = "Number of successess", xaxt = 'n', bty = 'n')
axis(side = 1, at = seq(0, 20, by = 2))
axis(side = 2, at = seq(0, -30, by = -5))
lines(0:20, single[[2]], lty = "dashed", col = "#ca0020")
lines(0:20, single2[[2]], type = 'p', col = "#0571b0", pch = 19)
lines(0:20, single2[[2]], lty = "dashed", col = "#0571b0")
title(main = "b) Log-likelihood for n = 20", font.main = 1, adj = 0.0)
legend("bottomright", col = c("#ca0020", "#0571b0"), fill = c(NA, NA), pch = c(19, 19), lty = c("dashed", "dashed"), legend = c("mSWEEP model", "Beta-binomial model"), bty = 'n', border = NA)
##dev.off()
