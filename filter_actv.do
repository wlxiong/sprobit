// import activity data
use ../BATS96/compACT, clear

******** filter recreational and maintenance activities
preserve
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
save dat/filtered_actv, replace

******** filter work-related activities
restore
// keep work activity in the first day
keep if dayno == 1
keep if actcode == 1

// get starting and ending time for each activity
sort sampno persno actno
gen tail = mod(actendhr,12)*60 + actendmn + cond(actendam==1,0,720)
gen leng = int(actdur/100)*60 + mod(actdur,100)
gen head = tail - leng

// get the working time span of individual
by sampno persno: egen begin  = min(head)
by sampno persno: egen ending = max(tail)

// save refined data
by sampno persno: keep if _n == _N
gen span = ending - begin
save dat/filtered_work, replace
