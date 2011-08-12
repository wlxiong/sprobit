******** generate data for the simulation study
version 11
cap log close
set more off
clear

******** global variables
global M = 3 // number of choice alternatives for each individual
global N = 4 // number of individuals in each household
global NM = $N*$M // number of observations in each household

// turn on log
clear
cd ~/Workspace/Stata/sprobit
log using simulate_data.log, replace

******** prepare the simulation data

// set number of observations in the simulation
local gobs 200
disp (`gobs'*$NM)
set obs 2400	 // `gobs'*$NM
set matsize 2400 // `gobs'*$NM

// true parameters
scalar _rho = 0.6 // 0.2
scalar _b1 = 1 // two independent variables
scalar _b2 = 2

// weight matrix
matrix Wk = ( 0, .5, .5, .5 \ .5,  0,  0,  0 \ .5,  0,  0,  0 \ .5,  0,  0,  0)
// matrix Wk = ( 0, .5, .5, .5 \ .5,  0,  0,  0 \ .5,  0,  0,  0 \ .5,  0,  0,  0)
// matrix Wk = ( 0, .5, .5, .5 \ .5,  0,  0,  0 \ .5,  0,  0,  0 \ .5,  0,  0,  0)
matrix W  = Wk # I($M)

// covariance matrix
matrix _A    = I($NM) - _rho*W
matrix _invA = invsym(_A)
matrix _covU = invsym(_invA'*_invA)

// generate household identity
gen sampno = int((_n-1)/$NM)+1
sort sampno

// generate correlated random errors
local randutils " "
forvalues j = 1/$NM {
	local randutils "`randutils' u_`j'"
}
drawnorm `randutils', cov(_covU)
forvalues j = 1/$NM {
	by sampno: replace  u_`j' = u_`j'[1]
}

// generate observed attributes
gen x1 = 5*runiform() - 3
gen x2 = 3*runiform() - 2

// generate latent utility
gen utils = _b1*x1 + _b2*x2
forvalues j = 1/$NM{
	by sampno: replace utils = utils + u_`j'
}

// generate observed choices
gen choice = utils > 0
tab choice

// save simulated data
save simulated_data, replace

// turn off log
log off
