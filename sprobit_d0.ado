******** likelihood function: derivative 0 form
capture program drop sprobit_d0
program define sprobit_d0
	version 11
	args func b lnf
	
	tempvar theta
	tempname $asrho sigma
	mleval `theta' = `b', eq(1)
	// get activity-specific auto-regressive coefficients
	local j = 1
	forvalues i = 1/$M {
		local j = `j' + 1
		mleval `rho_`i'' = `b', eq(`j') scalar
	}
	// for identification purpose, we assume sigma = 1
	scalar `sigma' = 1

	tempname A R invA covU L sqrtU Z corrU
	// define auto-regressive coefficient matrix
	matrix `R' = I($M)
	forvalues i = 1/$M {
		// the entries in matrix R are within the interval (0,1)
		matrix `R'[`i',`i'] = (`rho_`i'')
	}
	// calculate (I-R*W)^-1 and the covariance matrix
	matrix `A'    = I($NM) - (I($N)#`R')*W
	matrix `invA' = invsym(`A')
	matrix `covU' = invsym((`A')'*`A')*(`sigma'*`sigma')
	// square root of the diag of the covariance matrix
	matrix `sqrtU' = cholesky(diag(vecdiag(`covU')))
	// variance-normalizing diagonal matrix
	matrix `Z'     = invsym(`sqrtU')
	// normalized A^-1 and the correlation matrix
	matrix `invA'  = `Z'*`invA'
	matrix `corrU' = `Z'*`covU'*`Z'
	matrix `L'     = cholesky(`corrU')
	/*
		CHANGED consider i.i.d. normally distributed errors
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
		mkmat `xbs', matrix(`XB')
		matrix `AXB' = `XB'*(`invA')'
		// convert matrix into varlist
		capture drop axb*
		svmat `AXB', name(axb)
		tempvar last fi
		by $hid: gen double `last' = (_n==$NM)
		egen `fi'   = mvnp(axb*), chol(`L') dr($dr) prefix(z) signs(k*)
		mlsum `lnf' = ln(`fi') if `last'
	}
end
