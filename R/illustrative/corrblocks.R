

set.seed(1)
k <- 5
m <- 100
m1 <- 5
nn.pattern <- rep(0,m)
nn.pattern[sample(1:m, size = m1, replace=FALSE)] <- 1

mu <- 1.5 * nn.pattern

rho <- .995
z.block <- rnorm(m/k)
z <- mu + rho * rep(z.block, each=k) + sqrt(1-rho^2) * rnorm(m)


pdf("corrblocks.pdf", width=8, height=4)
par(mar=c(3.1, 4.5, 0.1, 0.1))
plot(z, xlab=expression(j), ylab=expression(hat(beta)[j]), type="n")
for(b in 0:(m/(2*k))) {
  polygon(x = .5+c(2*k*b,2*k*b,k*(2*b+1),k*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z)
dev.off()

pdf("corrblocks-red.pdf", width=8, height=4)
par(mar=c(3.1, 4.5, 0.1, 0.1))
plot(z, xlab=expression(j), ylab=expression(hat(beta)[j]), type="n")
for(b in 0:(m/(2*k))) {
  polygon(x = .5+c(2*k*b,2*k*b,k*(2*b+1),k*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z, col=1+nn.pattern, pch=1+15*nn.pattern, cex=1+0.5*nn.pattern)
dev.off()



set.seed(1)
block.size <- 10
n.blocks <- 10
m <- block.size * n.blocks
m1 <- 5
nn.indices <- sample(1:m, size = m1, replace=FALSE)
nn.pattern <- rep(0,m)
nn.pattern[nn.indices] <- 1

mu <- 1.5 * nn.pattern

rho <- .8
z.block <- rnorm(n.blocks)
z <- mu + rho * rep(z.block, each=block.size) + sqrt(1-rho^2) * rnorm(m)


wd <- 8
ht <- 5
mar.vec <- c(4.1, 4.7, 0.1, 0.1)


pdf("exercise1.pdf", width=wd, height=ht)
par(mar=mar.vec)
plot(rep(1:m,2), c(z,z+.1), xlab=expression(j), ylab=expression(hat(beta)[j]), type="n")
for(b in 0:(m/(2*block.size))) {
  polygon(x = .5+c(2*block.size*b,2*block.size*b,block.size*(2*b+1),block.size*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z)
dev.off()


pdf("exercise1-red.pdf", width=wd, height=ht)
par(mar=mar.vec)
plot(rep(1:m,2), c(z,z+.1), xlab=expression(j), ylab=expression(hat(beta)[j]), type="n")
for(b in 0:(m/(2*block.size))) {
  polygon(x = .5+c(2*block.size*b,2*block.size*b,block.size*(2*b+1),block.size*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z, col=1+nn.pattern, pch=1+15*nn.pattern, cex=1+0.5*nn.pattern)
dev.off()


pdf("exercise1-rank1.pdf", width=wd, height=ht)
par(mar=mar.vec)
plot(rep(1:m,2), c(z,z+.1), xlab=expression(j), ylab=expression(hat(beta)[j]), type="n")
for(b in 0:(m/(2*block.size))) {
  polygon(x = .5+c(2*block.size*b,2*block.size*b,block.size*(2*b+1),block.size*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z, col=1+nn.pattern, pch=1+15*nn.pattern, cex=1+0.5*nn.pattern)

ord <- order(abs(z), decreasing = TRUE)
rnk <- rank(-abs(z))
text(ord[1:5], z[ord[1:5]], labels=1:5, pos = 3)
dev.off()

pdf("exercise1-rank2.pdf", width=wd, height=ht)
par(mar=mar.vec)
plot(rep(1:m,2), c(z,z+.1), xlab=expression(j), ylab=expression(hat(beta)[j]), type="n")
for(b in 0:(m/(2*block.size))) {
  polygon(x = .5+c(2*block.size*b,2*block.size*b,block.size*(2*b+1),block.size*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z, col=1+nn.pattern, pch=1+15*nn.pattern, cex=1+0.5*nn.pattern)

ord <- order(abs(z), decreasing = TRUE)
rnk <- rank(-abs(z))
text(ord[1:5], z[ord[1:5]], labels=1:5, pos = 3)
text(nn.indices, z[nn.indices], labels = ifelse(rnk[nn.indices] <= 5, "", rnk[nn.indices]), pos=3)
dev.off()




grp.meds <- apply(matrix(z, nrow=block.size), MARGIN=2, median)
z2 <- z - rep(grp.meds, each=block.size)

pdf("exercise1-subtracted.pdf", width=wd, height=ht)
par(mar=mar.vec)
plot(rep(1:m,2), c(z2,z2+.1), xlab=expression(j), ylab=expression(hat(beta)[j] - "block median"), type="n")
for(b in 0:(m/(2*block.size))) {
  polygon(x = .5+c(2*block.size*b,2*block.size*b,block.size*(2*b+1),block.size*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z2, col=1+nn.pattern, pch=1+15*nn.pattern, cex=1+0.5*nn.pattern)

ord <- order(abs(z2), decreasing = TRUE)
rnk <- rank(-abs(z2))
text(ord[1:5], z2[ord[1:5]], labels=1:5, pos = 3)
text(nn.indices, z2[nn.indices], labels = ifelse(rnk[nn.indices] <= 5, "", rnk[nn.indices]), pos=3)
dev.off()




par(mar=c(3.1, 4.5, 0.1, 0.1))
plot(z2, xlab=expression(j), ylab=expression(hat(beta)[j]), type="n")
for(b in 0:(m/(2*block.size))) {
  polygon(x = .5+c(2*block.size*b,2*block.size*b,block.size*(2*b+1),block.size*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z2)
points(z2, col=1+nn.pattern, pch=1+15*nn.pattern, cex=1+0.5*nn.pattern)






p <- length(z)
n <- 5*p
blockSigma <- matrix(rho, block.size, block.size)
diag(blockSigma) <- 1
beta_cov <- as.matrix(diag(n.blocks) %x% blockSigma)
cov_mat <- solve(beta_cov)

R <- chol(cov_mat)
basis <- qr.Q(qr(matrix(rnorm(n*p), n)))
X <- basis %*% R

proj_y <- X %*% z

y <- rnorm(n, sd = 1)
y <- proj_y + y - lm(formula = y ~ X + 0)$fitted.values

library(glmnet)
lasso.cv <- cv.glmnet(X, y, nfolds = 10, intercept = F, standardize = F,
                      standardize.response = F, family = "gaussian")

lasso.fit <- glmnet::glmnet(X, y, lambda = 1*lasso.cv$lambda.min,
                            intercept = F, standardize = F,
                            standardize.response = F, family = "gaussian")

z2 <- as.vector(lasso.fit$beta)

pdf("exercise1-lasso.pdf", width=wd, height=ht)
par(mar=mar.vec)
plot(rep(1:m,2), c(z2,z2+.1), xlab=expression(j), ylab=expression(hat(beta)[j]("from lasso")), type="n")
for(b in 0:(m/(2*block.size))) {
  polygon(x = .5+c(2*block.size*b,2*block.size*b,block.size*(2*b+1),block.size*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z2, col=1+nn.pattern, pch=1+15*nn.pattern, cex=1+0.5*nn.pattern)

ord <- order(abs(z2), decreasing = TRUE)
rnk <- rank(-abs(z2))
text(ord[1:5], z2[ord[1:5]], labels=1:5, pos = 3)
text(nn.indices, z2[nn.indices], labels = ifelse(rnk[nn.indices] <= 5, "", rnk[nn.indices]), pos=3)
dev.off()



library(cknockoff)
result <- cknockoff(X, y, alpha = 0.2)
ckn.selected <- rep(NA, p)
ckn.selected[result$selected] <- z2[result$selected]

pdf("exercise1-ckn.pdf", width=wd, height=ht)
par(mar=mar.vec)
plot(rep(1:m,2), c(z2,z2+.1), xlab=expression(j), ylab=expression(hat(beta)[j]("from lasso")), type="n")
for(b in 0:(m/(2*block.size))) {
  polygon(x = .5+c(2*block.size*b,2*block.size*b,block.size*(2*b+1),block.size*(2*b+1)), y = c(-10,10,10,-10), col=rgb(.95,.95,.95), lty=3)
}
points(z2, col=1+nn.pattern, pch=1+15*nn.pattern, cex=1+0.5*nn.pattern)
points(ckn.selected, col=4, pch=1, cex=2, lwd=2)
dev.off()


