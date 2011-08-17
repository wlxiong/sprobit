******** spatial probit estimation for simulation studies
version 11
cap log close
set more off
set matsize 2400

******** read data
clear
cd ~/Workspace/Stata/sprobit
use simulated_data_$Wk_name
// turn on log
log using mle_example_$Wk_name.log, replace

******** define household identity
global hid sampno
sort sampno

******** deifne dependent and independent variables 
global y choice
global X x1 x2

******** get initial b0 from probit
set rmsg on
probit $y $X
set rmsg off
matrix b0 = e(b)

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
	ml model d0 sprobit_d0 (choice: $y = $X) /rho /lnsigma, tech(nr) ///
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
