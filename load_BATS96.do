/// make a subdir for the data
mkdir ../BATS96
cd ../BATS96
/// download BATS96 and unzip it
shell wget "ftp://ftp.abag.ca.gov/pub/mtc/planning/BATS/BATS96/BATS96G ASCII.zip"
shell unzip "BATS96G ASCII.zip"
/// import travel survey data
insheet using compACT.txt, clear
save compACT
insheet using compPER.txt, clear
save compPER
cd ..
*use BATS96/compACT.dta, clear
*merge m:1 sampno persno using BATS96/compPER.dta
