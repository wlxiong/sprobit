// import individual socio-economic data
use ../BATS96/compPER, clear
sort sampno

// replace missing variables
// replace age = . if age == 999
// replace income = . if income == 98 | income == 99
// remove missing observations
drop if age == 999
drop if income == 98 | income == 99
drop if employ == 3 | gender == 3 | student == 3

// only consider the household heads
keep if relate < 2
by sampno: gen numpers = _N
keep if numpers == 2
by sampno: gen has_spouse = sum(relate==1)
by sampno: replace has_spouse = has_spouse[_N]
keep if has_spouse == 1

// save refined data
save filtered_pers, replace
