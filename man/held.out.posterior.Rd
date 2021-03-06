% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/held-out.r
\name{held.out.posterior}
\alias{held.out.posterior}
\title{Calculate posterior covariance and estimate parameters for
held out genes given pseudotimes estimated by DeLorean model.}
\usage{
held.out.posterior(dl, held.out, posterior.sample = dl$best.sample)
}
\arguments{
\item{dl}{de.lorean object}

\item{held.out}{Held out gene expression levels}

\item{posterior.sample}{Posterior sample to use}
}
\description{
Calculate posterior covariance and estimate parameters for
held out genes given pseudotimes estimated by DeLorean model.
}
