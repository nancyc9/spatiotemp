#!/usr/bin/R


get_trend <- function(fname)
{
  pos <- regexpr('\\.t', fname)[[1]]
  tt <- substr(fname, pos+2, pos+2)
  return (as.numeric(tt))
}

get_pls <- function(fname)
{
  pos <- regexpr('\\.c', fname)[[1]]
  fname <- substr(fname, pos+2, nchar(fname))
  pos <- regexpr('\\.', fname)[[1]]
  comps <- substr(fname, 1, pos-1)
  return(as.numeric(strsplit(comps, split='')[[1]]))
}

get_beta <- function(fname)
{
  betas <- get_trend(fname) + 1
  my.cov.beta <- list(covf=c(rep(NA,betas)), nugget=rep(NA, betas))
  full.beta <- strsplit(fname,'\\.')[[1]][3]
  i=1
  while (nchar(full.beta))
  {
    b.structure <- regexpr('(exp|iid|matern|spherical)', full.beta)
    my.length <- attr(b.structure, 'match.length')
    my.cov.beta$covf[i] <- substr(full.beta, 1, my.length)
    full.beta <- substr(full.beta, (my.length+1), nchar(full.beta))
    b.nugget <- regexpr('(T|F)', full.beta)
    my.length <- attr(b.nugget, 'match.length')
    my.cov.beta$nugget[i] <- as.logical(substr(full.beta, 1, my.length))
    full.beta <- substr(full.beta, (my.length+1), nchar(full.beta))
    i <- i + 1
  }
  return(my.cov.beta)
}

get_nu <- function(fname)
{
  my.cov.nu <- list(covf=NA, nugget=NA, random.effect=NA)
  full.nu <- strsplit(fname,'\\.')[[1]][5]
  n.structure <- regexpr('(exp|iid|matern|spherical)', full.nu)
  my.length <- attr(n.structure, 'match.length')
  my.cov.nu$covf <- substr(full.nu, 1, my.length)
  full.nu <- substr(full.nu, (my.length+1), nchar(full.nu))
  n.re <- regexpr('(T|F)', full.nu)
  my.length <- attr(n.re, 'match.length')
  my.cov.nu$nugget <- as.logical(substr(full.nu, n.re[[1]], my.length+n.re[[1]]-1))
  my.cov.nu$random.effect <- as.logical(substr(full.nu, my.length+n.re[[1]], nchar(full.nu)))
  return(my.cov.nu)
}

