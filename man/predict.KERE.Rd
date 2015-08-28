\name{predict.KERE}
\alias{predict.KERE}
\title{
make predictions from a "KERE" object.}
\description{
Similar to other predict methods, this functions predicts fitted values and class labels from a fitted \code{\link{KERE}} object.
}
\usage{
\method{predict}{KERE}(object, kern, x, newx,...)
}
\arguments{
	\item{object}{fitted \code{\link{KERE}} model object.}
  	\item{kern}{the built-in kernel classes in \pkg{KERE}. Objects can be created by calling the rbfdot, polydot, tanhdot, vanilladot, anovadot, besseldot, laplacedot, splinedot functions etc. (see example.)
	}
  \item{x}{the original design matrix for training \code{KERE}.
}
\item{newx}{matrix of new values for \code{x} at which predictions are
to be made. NOTE: \code{newx} must be a matrix with each row as an observation. \code{predict} function does not accept a vector or other formats of \code{newx}.}
\item{\dots}{other parameters to \code{predict} function.}
}
\details{
The fitted \eqn{\alpha_0 + K * \alpha} at newx is returned as a size \code{nrow(newx)*length(lambda)} matrix for various lambda values where the \code{KERE} model was fitted.
}
\value{
The fitted \eqn{\alpha_0 + K * \alpha} is returned as a size \code{nrow(newx)*length(lambda)} matrix. The row represents the index for observations of newx. The column represents the index for the lambda sequence.
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
N <- 100
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

# create newx for prediction
N1 <- 5
X1 <- runif(N1)
X2 <- 2*runif(N1)
X3 <- 3*runif(N1)
newx <- cbind(X1,X2,X3)

# make prediction
p1 <- predict.KERE(m1, kern, X, newx)
p1
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
