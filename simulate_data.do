******** generate data for the simulation study
version 11
cap log close
set more off
clear

******** global variables
global M = 3 // number of choice alternatives for each individual
global N = 2 // number of individuals in each household
global NM = $N*$M // number of observations in each household

// turn on log
clear
cd ~/Workspace/Stata/sprobit
log using log/simulate_data.log, replace

******** set the simulation parameters

// set number of observations in the simulation
local gobs = 400
disp (`gobs'*$NM)
local nobs = `gobs'*$NM
set obs `nobs'	 // `gobs'*$NM
set matsize `nobs' // `gobs'*$NM

// true parameters
matrix _R = ( .7,  0,  0  \  0, -.5,  0  \  0,  0, 0.05)
matrix _R = I($N)#_R
scalar _b1 = -0.05
scalar _b2 = 0.3
scalar _b3 = 0.2
scalar _b4 = -0.1
scalar _cons1 = 5.0
scalar _cons2 = 2.0
scalar _cons3 = 6.0

// weight matrix
global Wk_name Wk1
matrix Wk1 = ( 0, 1 \ 1,  0)
// matrix Wk1 = ( 0, .5, .5, .5 \ .5,  0,  0,  0 \ .5,  0,  0,  0 \ .5,  0,  0,  0)
matrix Wk2 = ( 0, .5,  0, .5 \ .5,  0, .5,  0 \  0, .5,  0, .5 \ .5,  0, .5,  0)
matrix Wk3 = ( 0, .5, .5, .5 \ .5,  0, .5, .5 \ .5, .5,  0, .5 \ .5, .5, .5,  0)
matrix W = $Wk_name # I($M)

// covariance matrix
matrix _A    = I($NM) - _R*W
matrix _invA = invsym(_A)
matrix _covU = invsym(_invA'*_invA)

// generate household identity
gen sampno = int((_n-1)/$NM)+1
sort sampno

// generate person identity
by sampno: gen persno = int((_n-1)/$M)+1
sort sampno persno

// generate activity code
by sampno persno: gen actcode = _n
sort sampno persno actcode

******** generate correlated random errors
local randutils " "
forvalues j = 1/$NM {
	local randutils "`randutils' u_`j'"
}
// draw from multivariate normal distribution
drawnorm `randutils', cov(_covU)
matrix list _covU
corr u_*, covariance
// reformat random errors
gen uj = .
forvalues j = 1/$NM{
	by sampno: replace uj = u_`j'[1] if _n == `j'
}

******** generate latent utility
// generate socioeconomic attibutes
gen age    = 48*runiform() + 12
gen gender = cond(persno==1, 1, 0)
gen employ = runiform() > 0.1
// generate travel time
gen ttime = 55*runiform() + 5
sum age gender employ ttime
// calculate latent utility
gen _xb = _b1*age + _b2*gender + _b3*employ + _b4*ttime
forvalues k = 1/$M {
	by sampno persno: replace _xb = _xb + _cons`k' if _n == `k'
}
forvalues i = 1/$NM {
	by sampno: gen _xb`i' = _xb[`i']
}
// generate xb array
local _xbs " "
forvalues i = 1/$NM {
	local _xbs  "`_xbs'  _xb`i'"
}
// multiply _XB with (I-rho*W)
mkmat `_xbs', matrix(_XB)
matrix _AXB = _XB*(_invA)'
capture drop _axb*
svmat _AXB, name(_axb)
// add random errors
gen utils = .
forvalues j = 1/$NM{
	by sampno: replace utils = _axb`j' + uj if _n == `j'
}

******** generate observed choices
gen choice = utils > 0
tab actcode choice
tab gender choice
tab choice

// save simulated data
save dat/simulated_data_$Wk_name, replace

// turn off log
log off
