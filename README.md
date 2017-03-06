# WEST's Evidence of Absence Package

These are routines to estimate the number of missed entities after a series of
searches. Often these will be
used to estimate true number of carcassess on a wind power facility's grounds
after field searches. Implements the **Evidence of Absence (EoA)** model of Huso et al. (2015) and the
**Informed Evidence of Absence (IEoA)** approaches.

## How to git it

Assuming you have access to GitLab (i.e., you are inside the WEST network), issue the following 
commands in your git shell: 

`cd <directory you want>`
`git clone 'https://lar-git.west-inc.com/tmcdonald/evoab.git'`

The above commands will download all source from GitLab to your computer.  Among other things, 
you should see a `DESCRIPTION` file and `R` directory.  

#### Using `devtools`

Open R and `setwd()` to the directory containing the `DESCRIPTION` file. In R issure the following:

`library(devtools)`
`document()`
`install()`

#### Manual install

Open a command window, change directory to the folder containing `DESCRIPTION` and issue 
the following command: 

`r CMD INSTALL evoab`


## To Contribute

If you change something, and it's useful, issue a [*merge request* here.](https://lar-git.west-inc.com/tmcdonald/evoab/merge_requests)

## Usage Example

At this time, the main routine is `estimateL.EoA.MultiYear`.  Here is an example of how 
it is run: 


`# --- Fake data from a three year study`
`# The alpha and beta parameters for g-value Beta distributions, one per year. `
`g <- data.frame(
  alpha = c( 69.9299, 63.5035,  84.6997),
  beta = c(  736.4795,  318.3179, 759.9333 )
)`

`# The number of carcasses found each year
X <- c( 0, 1, 3)`

`# The regular un-informed eoa estimator
eoa <- estimateL.EoA.MultiYear( X, g, LMax=500 )  `

`# The informed eoa estimator
ieoa <- estimateL.EoA.MultiYear( X, g, Lprior="normal", Lprior.mean=20, Lprior.sd=4) `

`# To check convergence of the latter, run traceplot and Gelman stats
plot(ieoa$out) # tracePlot
gelman.diag(ieoa$out) # gelmanStats
gelman.plot(ieoa$out) # gelmanPlot
`
