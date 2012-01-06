******** likelihood function: derivative 0 form
capture program drop sprobit_d0
program define sprobit_d0
	version 11
	args func b lnf
	
	tempvar theta
	tempname rho lnsig sigma
	mleval `theta' = `b', eq(1)
	mleval `rho'   = `b', eq(2) scalar
	mleval `lnsig' = `b', eq(3) scalar
	scalar `sigma' = exp(`lnsig')

	// generate household working time span
	sort sampno persno
	by sampno: gen persid = _n
	forvalues j = 1/$N {
		by sampno: gen dist_`j' = max(ending[_n], ending[`j']) - ///
								  min(begin[_n], begin[`j'])
	}

	// generate weight matrix
	by sampno: gen persid = _n
	forvalues j = 1/$N {
		by sampno: gen weight_`i'`j' = exp(1/span_`i'`j')
	}

	tempname A invA covU L sqrtU Z corrU
	matrix `A'    = I($NM) - `rho'*W
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
