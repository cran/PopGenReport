#calculate mantel tests
wassermann <- function(gen.mat, cost.mat, eucl.mat, plot=TRUE)
{
mant1 <- mantel.partial(gen.mat, cost.mat, eucl.mat,permutations=9999)
mant2 <- mantel.partial(gen.mat, eucl.mat, cost.mat,permutations=9999)

mantel.mat <- data.frame(model=NA, r=NA, p=NA)
mantel.mat[1,] <- c("Gen ~ Friction | Eucl.Distance", round(mant1$statistic,4), round(mant1$signif,4))
mantel.mat[2,] <- c("Gen ~ Eucl.Distance | Friction", round(mant2$statistic,4), round(mant2$signif,4))

if (plot==TRUE)
{
plot(density(mant1), main="Gen ~ Friction | Eucl.Distance", cex.main=0.6)
plot(density(mant2), main="Gen ~ Eucl.Distance | Friction", cex.main=0.6)
}

return(list(mantel.mat=mantel.mat))
}
