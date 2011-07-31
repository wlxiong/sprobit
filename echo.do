// capture program drop echo
program define echo
	version 9.1
	if `"`0'"' != "" display as text `"`0'"'
end
