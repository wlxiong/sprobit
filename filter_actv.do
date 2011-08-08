// import activity data
use ../BATS96/compACT, clear

// drop fields
keep sampno persno dayno actcode
// only consider a subset of activities and the first day
keep if dayno == 1
keep if actcode == 14 | actcode == 15 | actcode == 16 | actcode == 21

// remove duplicates
duplicates drop sampno persno dayno actcode, force
tab actcode

// prepare binary data strcuture for estimation function
gen choice = 1
reshape wide choice, i(sampno persno dayno) j(actcode)
mvencode choice*, mv(0)
reshape long

// save refined data
save filtered_actv, replace
