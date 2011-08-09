******** spatial probit model of intra-household interactions
version 11
cap log close
set more off
set matsize 10000

******** likelihood function: derivative 0 form
capture program drop sprobit_d0
program define sprobit_d0
	version 11
	args todo b lnf
	
	tempvar theta
	tempname rho
	mleval `theta' = `b', eq(1)
	mleval `rho'   = `b', eq(2) scalar

	tempname A invA covU L
	matrix `A'    = I($NM) - `rho'*W
	matrix `invA' = inv(`A')
	matrix `covU' = inv(`A'*`A')
	matrix `L'    = cholesky(`covU')
	/*
		TODO consider unrestricted covariance matrix
	*/
	// declare xb*
	forvalues i = 1/$NM {
		tempvar  xb`i'
	}
   	// generate k* and xb*
	capture drop k*
	quietly{
		forvalues i = 1/$NM {
			by $hid: gen double k`i'  = (2*$ML_y1[`i']) - 1
			by $hid: gen double `xb`i'' = `theta'[`i']
		}
		// generate xb array
		local xbs " "
		forvalues i = 1/$NM {
			local xbs  "`xbs'  `xb`i''"
		}
		// multiply XB with (I-rho*W)
		tempname XB AXB
		mkmat `xbs', matrix(XB)
		matrix AXB = XB*(`invA')'
		// convert matrix into varlist
		capture drop axb*
		svmat AXB, name(axb)
		tempvar last fi
		by $hid: gen byte `last' = (_n==$NM)
		egen `fi'   = mvnp(axb*), chol(`L') dr($dr) prefix(z) signs(k*)
		mlsum `lnf' = ln(`fi') if `last'
	}
end

// turn on log
log using sprobit_ml.log, replace

******** read data
clear
cd ~/Workspace/Stata/sprobit
use merged_actv_pers
sort sampno persno

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
	TODO define economic distance matrix
*/

******** define household identity
global hid sampno

******** deifne dependent and independent variables 
global y choice
global X age i.gender i.employ i.student
/*
	TODO remove i.student from independent variables
	TODO and also consider significance of _cons
	TODO extend to alternative specified parameters
*/

******** get initial b0 from probit
set rmsg on
probit $y $X
set rmsg off
matrix b0 = e(b)

******** estimation procedure: increase `drnum' gradually
local drlist 10 20 50
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
