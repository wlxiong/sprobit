// import individual socio-economic data
use ../BATS96/compPER, clear
sort sampno

// merge work activity
merge 1:1 sampno persno using filtered_work
tab _merge
keep if _merge == 3
drop _merge

******** filter household heads
// replace missing variables
// replace age = . if age == 999
// replace income = . if income == 98 | income == 99

// remove missing observations
drop if age == 999
drop if income == 98 | income == 99
drop if employ == 3 | gender == 3 | student == 3
drop if license == 3

// only consider the household heads
keep if relate < 2
by sampno: gen numpers = _N
keep if numpers == 2
by sampno: gen has_spouse = sum(relate==1)
by sampno: replace has_spouse = has_spouse[_N]
keep if has_spouse == 1

******** calculate economic distance
// number of persons in each household
qui tab relate
global N = r(r)

// generate household working time span
sort sampno persno
by sampno: gen persid = _n
forvalues j = 1/$N {
	by sampno: gen dist_`j' = max(ending[_n], ending[`j']) - ///
							  min(begin[_n], begin[`j'])
}

// save refined data
save filtered_pers, replace
