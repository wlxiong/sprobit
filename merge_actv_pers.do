// load activity data
use filtered_actv, clear
// merge activity data with individual's socio-economic attributes
merge m:1 sampno persno using filtered_pers, ///
	keepusing(relate gender age employ student license numpers has_spouse)
tab _merge
keep if _merge == 3
drop _merge
sort sampno
// if hasing missing values, drop the household
egen has_miss = rowmiss(_all)
by sampno: replace has_miss = sum(has_miss)
drop if has_miss
// if only one person, also drop the household
by sampno: gen numactv = sum(actcode < .)
by sampno: replace numactv = numactv[_N]
keep if numactv == 8
// save merged data
save merged_actv_pers, replace
