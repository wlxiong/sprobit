// load activity data
use dat/filtered_actv, clear
// merge activity data with individual's socio-economic attributes
merge m:1 sampno persno using dat/filtered_pers, ///
	keepusing(persid relate gender age employ student license numpers has_spouse dist_*)
tab _merge
keep if _merge == 3
drop _merge

// if hasing missing values, drop the household
sort sampno
egen has_miss = rowmiss(_all)
by sampno: replace has_miss = sum(has_miss)
drop if has_miss
drop has_miss
// if only one person, also drop the household
by sampno: gen numactv = sum(actcode < .)
by sampno: replace numactv = numactv[_N]
keep if numactv == 8
drop numactv
// save merged data
save dat/merged_actv_pers, replace
