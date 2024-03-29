% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetworkUtils.R
\name{BuildMinConnectedGraphs}
\alias{BuildMinConnectedGraphs}
\title{Compute minimum connected network composed of seed nodes using shortest path based approach}
\usage{
BuildMinConnectedGraphs(dataSetObj = NA, max.len = 200)
}
\arguments{
\item{dataSetObj}{Input the name of the created dataSetObj (see Init.Data)}

\item{max.len}{Maximum number of seeds, if more than this number, take top 100 seed nodes based on degrees}
}
\value{
check overall.graph object for result
}
\description{
Compute minimum connected network composed of seed nodes using shortest path based approach
}
