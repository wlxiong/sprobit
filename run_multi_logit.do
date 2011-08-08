/// import activity data
use ../BATS96/compACT, clear
/// drop fields
keep sampno persno dayno actcode
keep if actcode == 14 | actcode == 15 | actcode == 16 | actcode == 21
/// remove duplicates
duplicates drop sampno persno dayno, force
/// prepare binary data strcuture for estimation function
merge m:1 sampno persno using ../BATS96/compPER, ///
	keepusing(relate gender age employ student income)
/// remove missing observations
drop if employ == 3 | gender == 3 | student == 3
/// only consider the household heads
drop if relate > 1
/// replace missing variable 
replace age = . if age == 999
replace income = . if income == 98 | income == 99
/// the multinomial logit model
mlogit actcode age i.gender i.employ i.student
/// the multinomial probit model
mprobit actcode age i.gender i.employ i.student
