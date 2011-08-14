******** spatial probit model of intra-household interactions
/*
	TODO extend to alternative specified parameters
	TODO identification of alternative specified parameters and 
		covariance matrix of multivariate probit (not MNP)
*/
version 11
cap log close
set more off
set matsize 10000

******** read data
clear
cd ~/Workspace/Stata/sprobit
use merged_actv_pers
// turn on log
log using mle_BATS96.log, replace

******** global varibles
qui tab actcode
global M = r(r)
qui tab relate
global N = r(r)
global NM = $N*$M
disp $NM
matrix Wk = (0.0, 0.5 \ 0.5, 0.0)
matrix W  = Wk # I($M)
/*
	TODO Define economic distance matrix: 
		 Differences of total working hours between husband and wife
	NOTE The spatial weight matrix W is highly sparse. 
	NOTE The resulting weights are meaningful, finite and non-negative.
	NOTE It is also important to maintain the weights matrix as exogenous.
	NOTE When a scaling factor, such as a spatial autoregressive coefficient, 
		is included together with parameterized weights, both sets of parameters 
		are not necessarily identified.
*/

******** define household identity
global hid sampno
sort sampno persno

******** deifne dependent and independent variables 
global y choice
global X age i.gender i.employ
/*
	CHANGED remove i.student from independent variables
	CHANGED and also consider significance of _cons after removal of i.student
	TODO add drive license as independent variables (or in weight matrix?)
	TODO Heteroscedastic covariance matrix: 
			The full spatial model is premultiplied by the variance-normalizing 
			transformation diagonal matrix.
	TODO number of cars within the household
	TODO number of children under 6-year old 6 to 12-year
*/

******** get initial b0 from probit
set rmsg on
probit $y $X
set rmsg off
matrix b0 = e(b)

******** estimation procedure: increase `drnum' gradually
local drlist 10 20 50 100
foreach drnum of local drlist {
	// create `drnum' Halton draws
	set rmsg on
	mdraws, dr(`drnum') neq($NM) prefix(z) burn(15) antithetics replace
	set rmsg off
	return list
	global dr = r(n_draws)

	// call simulation-based ML
	ml model d0 sprobit_d0 (choice: $y = $X) /r, tech(nr) ///
	title(Spatial Probit Model, $dr Random Draws)
	ml init b0
	set rmsg on
	ml maximize, difficult
	set rmsg off

	// update initial b0
	matrix b0 = e(b)
}

// turn off log
log off
