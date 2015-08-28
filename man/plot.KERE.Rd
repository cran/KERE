\name{plot.KERE}
\alias{plot.KERE}

\title{Plot coefficients from a "KERE" object
}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{KERE}} object.
}
\usage{
\method{plot}{KERE}(x, ...)
}
\arguments{
\item{x}{fitted \code{\link{KERE}} model.}
\item{\dots}{other graphical parameters to plot.}
}
\details{
A coefficient profile plot is produced. The x-axis is \eqn{log(\lambda)}. The y-axis is the value of fitted \eqn{\alpha}'s.
}

\author{
Yi Yang, Teng Zhang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>
}

\references{
Y. Yang, T. Zhang, and H. Zou. "Flexible Expectile Regression in Reproducing Kernel Hilbert Space." ArXiv e-prints: stat.ME/1508.05987, August 2015.
}


\examples{
# create data
N <- 200
X1 <- runif(N)
X2 <- 2*runif(N)
X3 <- 3*runif(N)
SNR <- 10 # signal-to-noise ratio
Y <- X1**1.5 + 2 * (X2**.5) + X1*X3
sigma <- sqrt(var(Y)/SNR)
Y <- Y + X2*rnorm(N,0,sigma)
X <- cbind(X1,X2,X3)

# set gaussian kernel 
kern <- rbfdot(sigma=0.1)

# define lambda sequence
lambda <- exp(seq(log(0.5),log(0.01),len=10))

# run KERE
m1 <- KERE(x=X, y=Y, kern=kern, lambda = lambda, omega = 0.5) 

# plot the solution paths
plot(m1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
