% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_bfcif.R
\name{get_args}
\alias{get_args}
\title{get_args}
\usage{
get_args(env = parent.frame(), to_ignore = character(), ...)
}
\arguments{
\item{env}{Environment to collect variable from. By default this is the
calling environment, i.e. the function containing a call to get_args. This
may also be a function that returns an environment.}

\item{to_ignore}{character names of variables to ignore.}

\item{...}{Additional variables that should be considered as part of
environment.}
}
\value{
parameters of calling function as a named list.
}
\description{
returns parameters of calling function as a named list.
}
\examples{
#The most common usage is to simply collect all local variables in a function
test_fun = function(x = 1, y = 2){
  get_args()
}
test_fun()

#Specified variables may be ignored
test_fun2 = function(x = 1, y = 2){
  get_args(to_ignore = "x")
}
test_fun2()

#Additional variables can also be added from higher environments
global_z = 3
test_fun3 = function(x = 1, y = 2){
  get_args(env = parent.frame, to_ignore = character(), z = global_z)
}
test_fun3()
}
