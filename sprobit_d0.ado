******** likelihood function: derivative 0 form
capture program drop sprobit_d0
program define sprobit_d0
	version 11
	args func b lnf
	
	tempvar theta
	tempname lnsig sigma
	mleval `theta' = `b', eq(1)
	mleval `lnsig' = `b', eq(2) scalar
	scalar `sigma' = exp(`lnsig')
	// get activity-specific auto-regressive coefficients
	local num_eq = $M+2
	forvalues i = 3/`num_eq' {
		local ri = `i' - 2
		tempname rho_`ri'
		mleval `rho_`ri'' = `b', eq(`i') scalar
	}

	tempname A R invA covU L sqrtU Z corrU
	// define auto-regressive coefficient matrix
	matrix `R' = I($M)
	forvalues i = 1/$M {
		matrix `R'[`i',`i'] = (`rho_`i'')
	}
	// calculate (I-R*W)^-1 and the covariance matrix
	matrix `A'    = I($NM) - (I($N)#`R')*W
	matrix `invA' = invsym(`A')
	matrix `covU' = invsym((`A')'*`A')*`sigma'*`sigma'
	// square root of the diag of the covariance matrix
	matrix `sqrtU' = cholesky(diag(vecdiag(`covU')))
	// variance-normalizing diagonal matrix
	matrix `Z'     = invsym(`sqrtU')
	// normalized A^-1 and the correlation matrix
	matrix `invA'  = `Z'*`invA'
	matrix `corrU' = `Z'*`covU'*`Z'
	matrix `L'     = cholesky(`corrU')
	/*
		CHANGED consider i.i.d. normal distributed errors
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
		by $hid: gen byte `last' = (_n==$NM)
		egen `fi'   = mvnp(axb*), chol(`L') dr($dr) prefix(z) signs(k*)
		mlsum `lnf' = ln(`fi') if `last'
	}
end
