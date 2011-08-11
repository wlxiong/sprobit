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
	matrix `covU' = inv((`A')'*`A')
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
