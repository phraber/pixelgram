% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{color.edges}
\alias{color.edges}
\title{Recursively propagate tip colors internally to color interior branches}
\usage{
color.edges(T, tip_colors)
}
\arguments{
\item{T}{ape::phylo() tree object.}

\item{tip_colors}{vector of colors, whether rgb hex or integers, whose order matches the names in the tree tip.label slot.}
}
\value{
A vector of colors of length equal to the number of tree branches.
}
\description{
Branches with descendents that have different colors will be represented with the default color.
}
\details{
NB: the default color is currently hard-coded, which seems unwise.
Some trees may not be amenable to the recursion scheme used here.
Alternatively, starting to label on the root node may be a bad idea.
One tree was reported to have caused a problem, so it is worth keeping in mind.
}
