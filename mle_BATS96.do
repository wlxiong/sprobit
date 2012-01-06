******** spatial probit model of intra-household interactions
version 11
cap log close
set more off
set matsize 6000

******** read data
clear
cd ~/Workspace/Stata/sprobit
use dat/merged_actv_pers
// turn on log
log using log/mle_BATS96.log, replace

******** global varibles
qui tab relate
global N = r(r)
qui tab actcode
global M = r(r)
global NM = $N*$M
disp $NM
matrix Wk = (0.0, 1.0 \ 1.0, 0.0)
matrix W  = Wk # I($M)
/*
	TODO Define economic distance matrix: 
		 Differences of total working hours between husband and wife
		 Or the non-overlapping working hours
	NOTE The spatial weight matrix W is highly sparse. 
	NOTE The resulting weights are meaningful, finite and non-negative.
	NOTE It is also important to maintain the weights matrix as exogenous.
	NOTE When a scaling factor, such as a spatial autoregressive coefficient, 
		is included together with parameterized weights, both sets of parameters 
		are not necessarily identified.
*/

******** define household identity
global hid sampno
sort sampno persno actcode

/*
	TODO extend to alternative specified parameters
	TODO identification of alternative specified parameters and 
		 covariance matrix of multivariate probit (not MNP)
*/
// generate alternative specific constants and AR coefficients
global ascons " "
global asrho " "
global rhoeq " "
local j = 0
foreach actv of numlist 14 15 16 21 { // omit activity 21 to avoid collinearity
	gen _cons_`actv' = cond(`actv'==actcode, 1, 0)
	global ascons "$ascons _cons_`actv'"
	local j = `j' + 1
	global asrho "$asrho rho_`j'"
	global rhoeq "$rhoeq /rho_`actv'"
}
disp "$ascons"
disp "$asrho"
disp "$rhoeq"

******** deifne dependent and independent variables 
global y choice
global X age i.gender i.employ $ascons
/*
	CHANGED remove i.student from independent variables
	CHANGED and also consider significance of _cons after removal of i.student
	CHANGED add drive license as independent variables (or in weight matrix?)
	CHANGED remove drive license (not siginificant)
	CHANGED Heteroscedastic covariance matrix: 
			The full spatial model is premultiplied by the variance-normalizing 
			transformation diagonal matrix.
	TODO number of cars within the household
	TODO number of children under 6-year old 6 to 12-year
*/

******** get initial b0 from probit
set rmsg on
probit $y $X, nocons
set rmsg off
matrix b0 = e(b)
matrix list b0

******** estimation procedure: increase `drnum' gradually
local drlist 2 10 20 50
foreach drnum of local drlist {
	// create `drnum' Halton draws
	set rmsg on
	mdraws, dr(`drnum') neq($NM) prefix(z) burn(15) antithetics replace
	set rmsg off
	return list
	global dr = r(n_draws)

	// call simulation-based ML
	ml model d0 sprobit_d0 (choice: $y = $X, nocons) $rhoeq, tech(nr 4 dfp 8) ///
	title(Spatial Probit Model, $dr Random Draws)

	disp "run simulated maximum likelihood"
	// set initial values
	ml init b0
	// start the numerical method
	set rmsg on
	ml maximize, difficult
	set rmsg off
	// update initial b0
	matrix b0 = e(b)
}

// turn off log
log off
