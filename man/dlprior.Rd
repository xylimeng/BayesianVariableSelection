\name{dlprior}
\alias{dlprior}
\title{Function to simulate the dirichlet laplace shrinkage prior}
\usage{
dlprior(hyper=1/2,p=10000000,plt=TRUE,min=-5,max=5,sigma=1)
}

\arguments{

\item{hyper}{important hyperparameter that related to posterior shrinkage scales and prior distribution.}
\item{plt}{whether to plot the dirichlet laplace prior. default TRUE means plot the distribution.}
\item{min}{left point of the plot graph.}
\item{max}{right point of the plot graph.}
\item{sigma}{the value equals to normal noises' standard deviations}
}

\value{
\item{beta}{A p*1 vector. p observations from the distribution}
}

\description{
This function generates random deviates from dirichlet laplace shrinkage prior and can plot the distribution function.
}
\examples{
theta=dlprior(hyper=1/2,p=10000000,plt=TRUE,min=-5,max=5,sigma=1)
}
