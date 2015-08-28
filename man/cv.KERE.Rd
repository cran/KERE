\name{cv.KERE}
\alias{cv.KERE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cross-validation for KERE
}
\description{
Does k-fold cross-validation for \code{KERE}, produces a plot, and returns a value for \code{lambda}.
}
\usage{
\method{cv}{KERE}(x, y, kern, lambda = NULL, nfolds = 5, foldid, omega = 0.5, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}
  \item{y}{response variable.
}
\item{kern}{the built-in kernel classes in \pkg{KERE}. 
The \code{kern} parameter can be set to any function, of class
kernel, which computes the inner product in feature space between two
vector arguments. \pkg{KERE} provides the most popular kernel functions
which can be initialized by using the following
functions:
\itemize{
\item\code{rbfdot} Radial Basis kernel function,
\item\code{polydot} Polynomial kernel function,
\item\code{vanilladot} Linear kernel function,
\item\code{tanhdot} Hyperbolic tangent kernel function,
\item\code{laplacedot} Laplacian kernel function,
\item\code{besseldot} Bessel kernel function,
\item\code{anovadot} ANOVA RBF kernel function,
\item\code{splinedot} the Spline kernel. 
}
Objects can be created by calling the rbfdot, polydot, tanhdot, vanilladot, anovadot, besseldot, laplacedot, splinedot functions etc. (see example.)
}
\item{lambda}{a user supplied \code{lambda} sequence. It is better to supply a decreasing sequence of \code{lambda} values, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.
}
  \item{nfolds}{
	number of folds - default is 5. Although \code{nfolds}
	can be as large as the sample size (leave-one-out CV), it is not
	recommended for large datasets. Smallest value allowable is \code{nfolds=3}.
	}
  \item{foldid}{	an optional vector of values between 1 and \code{nfold}
	identifying what fold each observation is in. If supplied,
	\code{nfold} can be missing.}
  \item{omega}{the parameter \eqn{\omega}{omega} in the expectile regression model. The value must be in (0,1). Default is 0.5.
}
  \item{\dots}{other arguments that can be passed to \code{KERE}.}
}
\details{
The function runs \code{\link{KERE}} \code{nfolds}+1 times; the
first to get the \code{lambda} sequence, and then the remainder to
compute the fit with each of the folds omitted. The average error and standard deviation over the
folds are computed.
}
\value{
an object of class \code{\link{cv.KERE}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length
		\code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlo}{lower curve = \code{cvm-cvsd}.}
		\item{name}{a character string "Expectile Loss"}
		\item{lambda.min}{the optimal value of \code{lambda} that gives minimum
		cross validation error \code{cvm}.}
		\item{cvm.min}{the minimum
		cross validation error \code{cvm}.}
}
\author{
Yi Yang, Teng Zhang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>
}

\references{
Y. Yang, T. Zhang, and H. Zou. "Flexible Expectile Regression in Reproducing Kernel Hilbert Space." ArXiv e-prints: stat.ME/1508.05987, August 2015.
}
\examples{

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

cv.KERE(x=X, y=Y, kern, lambda = lambda, nfolds = 5, omega = 0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
