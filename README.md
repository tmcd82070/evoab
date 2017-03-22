# WEST's Evidence of Absence Package

These are routines to estimate the number of missed entities after a series of
searches. Often these will be
used to estimate true number of carcassess on a wind power facility's grounds
after field searches. Implements the **Evidence of Absence (EoA)** model of Huso et al. (2015) and the
**Informed Evidence of Absence (IEoA)** approaches.

## How to git it:

Assuming you have access to GitLab (i.e., you are inside the WEST network), issue the following 
commands in your git shell: 

`cd <directory you want>`  
`git clone 'https://lar-git.west-inc.com/tmcdonald/evoab.git'`  

The above commands will download all source from GitLab to your computer.  Among other things, 
you should see a `DESCRIPTION` file and `R` directory.  

## Intalling:

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

At this time, the main routine is `eoa`.  It takes a count vector, model for lambda, and g-values, 
Here is an example of how 
it is run: 

This is fake data from a three year study on seven sites.  First, the 
alpha and beta parameters for g-value distributions, one per year.   
`ns <- 3  
ny <- 7  
g <- data.frame(  
  alpha = rnorm(ns*ny,70,2),  
  beta = rnorm(ns*ny,700,25)  
)`  

This is the carcasses count vector, one count per site per year.  

`Y <- rbinom(ns*ny, c(rep(20,ny), rep(40,ny), rep(60,ny)), g$alpha/(g$alpha+g$beta))`

This is the covariate data frame.  This data frame contains Year as a linear 
effect (1,2,3,etc) and Year as a factor (2015, 2016, 2017, etc).  

`df <- data.frame(year=factor(c(rep("2015",ny),rep("2016",ny),rep("2017",ny))),  
    Year=c(rep(1,ny),rep(2,ny),rep(3,ny)))`  

The following computes un-informed EoA (vague priors for coefficients):     

`eoa.1 <- eoa(Y~year, g, df )`

This computes IEoA:


`# Assume prior mean is 10 and prior sd is 3`  
`# Fit intercept-only model to get one mean lambda
prior <- data.frame(mean=log(10), sd=log(3), row.names="(Intercept)")
eoa.1 <- eoa(Y~1, g, df, priors=prior )`


After either run, you should check convergence.  
To do so, run a traceplot and Gelman stats.  Any Rhats > 1.1 indicate suspect 
convergence. 

`library(lattice)  
xyplot(ieoa.1$out[,labels(ieoa.1)])  
acfplot(ieoa.1$out[,labels(ieoa.1)])  
densityplot(ieoa.1$out[,labels(ieoa.1)])  
gelman.diag(ieoa.1$out) # gelmanStats  
gelman.plot(ieoa.1$out) # gelmanPlot`
