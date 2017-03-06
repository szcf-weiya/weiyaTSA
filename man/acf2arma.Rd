\name{acf2arma}
\alias{acf2arma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine ARMA model by ACF
%%  ~~function to do ... ~~
}
\description{
Given Auto-Covariance Function to determine the parameter of ARMA model
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
acf2arma(x,p,q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the sequence of auto-covariance function values}
  \item{p}{the degree of AR model}
  \item{q}{the degree of MA model}
%%     ~~Describe \code{x} here~~
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ARcoeff}{the coefficients of AR model}
\item{MAcoeff}{the coefficients of MA model}
\item{VarNoise}{the variance of noise}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
weiya
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x = c(4.61, -1.06, 0.29, 0.69, -0.12)
acf2arma(x,2,2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Time series}% use one of  RShowDoc("KEYWORDS")
%\keyword{}% __ONLY ONE__ keyword per line
