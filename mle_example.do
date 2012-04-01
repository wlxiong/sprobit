******** spatial probit estimation for simulation studies
version 11
cap log close
set more off
set matsize 5000
global Wk_name Wk1

******** read data
clear
use dat/simulated_data_$Wk_name
// turn on log
log using log/mle_example_$Wk_name.log, replace

******** define household identity
global hid sampno
sort sampno

// generate alternative specific constants and AR coefficients
global ascons " "
global asrho " "
global rhoeq " "
local j = 0
foreach actv of numlist 1 2 3 { // omit activity 21 to avoid collinearity
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
global X age gender employ ttime $ascons

******** get initial b0 from probit
set rmsg on
probit $y $X, nocons
set rmsg off
matrix b0 = e(b)
matrix list b0

// set trace on
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
// set trace off

// turn off log
log off
