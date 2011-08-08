// read data
clear
do read_data.do

// calculate choice frequencies
sort sampno persno dayno
by sampno persno dayno: gen chosum = sum( choice )
by sampno persno dayno: replace chosum = . if _n < _N
tab chosum
gen chosen = actcode*choice
replace chosen = . if chosen == 0
tab chosen
// create an indentity
egen spdid = concat(sampno persno dayno)
destring spdid, replace

// conditional logit model
clogit choice age i.gender i.employ i.student, group(id)
// alternative specified conditional logit model
xi: asclogit choice, case(id) //
	casevars(age i.gender i.employ i.student) ///
	alternatives(actcode) basealternative(14)
// alternative specified multinomial probit model
// xi: asmprobit choice, case(id) //
// 	casevars(age i.gender i.employ i.student) ///
// 	alternatives(actcode) basealternative(14) intmethod(halton) correlation(independent)
