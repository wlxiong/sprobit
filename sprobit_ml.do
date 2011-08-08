// spatial probit model
version 11
cap log close
set more off
set matsize 10000

// likelihood function

// linear form
capture program drop sprobit_lf
program define sprobit_lf

end

// deravative 0 form
capture program drop sprobit_d0
program define sprobit_d0
	/// ml general form form 0
	version 11
	args todo b lnf
	
	tempvar theta
	tempname rho
	mleval `theta' = `b', eq(1)
	mleval `rho'   = `b', eq(2) scalar

	tempname A invA covU L
	// display "matrix oeprations 0"
	// matrix list W
	matrix `A'    = I($NM) - `rho'*W
	// matrix list `A'
	matrix `invA' = inv(`A')
	// matrix list `invA'
	matrix `covU' = inv(`A'*`A')
	// matrix list `covU'
	matrix `L'    = cholesky(`covU')
	// matrix list `L'
	
	// display "matrix oeprations 1"
	
	forvalues i = 1/$NM {
		tempvar  xb`i'
		tempname xbvec`i'
	}
	// display "tempvar over"
	// list `theta' in 1/10
	quietly{
		forvalues i = 1/$NM {
			by $hid: gen double k`i'  = (2*$ML_y1[`i']) - 1
			by $hid: gen double `xb`i'' = `theta'[`i']
			mkmat `xb`i'', matrix(`xbvec`i'')
//			by $hid: list `theta'[`i']
		}
		// disp `" xbs: "$xbs" "'
		// list $xbs in 1/20
		// sum $xbs
		// mkmat $xbs, matrix(XB)
		tempname XB AXB
		matrix XB = `xbvec1'
		forvalues i = 2/$NM {
			matrix XB = XB, `xbvec`i''
		}
		// matrix dir 
		// matrix list XB
		matrix AXB = XB*(`invA')'
		// matrix list AXB
		svmat AXB, name(axb)
		tempvar last fi
		by $hid: gen byte `last' = (_n==$NM)
		su axb* k*
		egen `fi'   = mvnp(axb*), chol(`L') dr($dr) prefix(z) signs(k*)
		mlsum `lnf' = ln(`fi') if `last'
		// clear variables
		drop axb* k*
	}	
end

// turn on log
log using sprobit_ml.log, replace

// read data
clear
use merged_actv_pers
sort sampno persno

// global varibles
qui tab actcode
global M = r(r)
qui tab relate
global N = r(r)
global NM = $N*$M
disp $NM
matrix Wk = (0.0, 0.5 \ 0.5, 0.0)
matrix W  = Wk # I($M)

global hid sampno
global y choice
global X age i.gender i.employ i.student

mdraws, dr(50) neq($NM) prefix(z) replace
global dr = r(n_draws)

// get initial value from probit
probit $y $X
matrix b0 = e(b)

// generate xb array
// global xbs " "
// forvalues i = 1/$NM {
// 	global xbs "$xbs \`xb`i''"
// }

// call simulation-based ML
ml model d0 sprobit_d0 (choice: $y = $X) /r
ml init b0
ml maximize, trace

// turn off log
log off
