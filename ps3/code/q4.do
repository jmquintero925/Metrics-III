**PSET 3 QUESTION 4
*** Isaac Norwich

**If you are reading this, don't worry! I will be updating the formatting when I don't have to deal with the annoying remote server.


net install scheme-modern, from("https://raw.githubusercontent.com/mdroste/stata-scheme-modern/master/")
set scheme modern

ssc install estout
ssc install tabout

cap net install sg97_5.pkg
*For frmttable

cd "C:\Users\inorwich\OneDrive - The University of Chicago\ps3\code"
use "..\data\lottery.dta", clear

* B) Assess instrument relevance
/* To test instrument relevance we can regress D on Z and test if the 
coefficient on Z is 0 or not */
eststo: reg d z
esttab using "..\tables\q4_partB.tex", se nostar label nonum tex replace
/* Ok we get a coefficient of .5203044 with SE of .02 */

* C) Estimate the return to medical school on earnings in 2007 using IV.
*Two ways, 2SLS or IVREGRESS
eststo clear
eststo: ivregress 2sls lnw (d = z)
* Do we need robust standard errors? Slide 21 of week3.pdf has robust in the command

reg d z
predict d_hat
la var d_hat "$\widehat{D}$ from D on Z"
eststo: reg lnw d_hat
* Weeee exact same coefficient and std error is basically the same

* Gives coefficient on d_hat of .187
esttab using "..\tables\q4_partC.tex", se nostar label nonum mtitles("2SLS manually" "2SLS using ivregress") tex replace


* D) Count the number of compliers. Compare to number of applicants by gender.
tab d z

**MHE pg 167 E[D|Z=1]-E[D|Z=0]
summ d if z == 1
local d_z1 = `r(mean)'
summ d if z == 0
local d_z0 = `r(mean)'
di "The proportion of compliers in the sample is " %04.3f `d_z1'-`d_z0'
di "The number of compliers in the sample is " %7.3f (`d_z1'-`d_z0')*_N

cap matrix drop partD
foreach type in 0 1 2  {
	** So outcome is # of female, number of male, and total compliers
	qui summ d if z==1 & female != `type'
	local stuff1 = `r(mean)'
	qui summ d if z==0 & female != `type'
	local stuff2 = `r(mean)'
	qui count if female != `type'
	local count_`type' = `r(N)'
	local prob_complier =  `stuff1' - `stuff2'
	local total_compliers_`type' = `prob_complier' * `r(N)'
	di `total_compliers_`type''
	matrix partD = nullmat(partD) \ `prob_complier', `total_compliers_`type''
}
matrix rownames partD = Female Male Total
matrix colnames partD = $\%$ N
esttab matrix(partD, fmt(%04.3f %7.3f)) using "..\tables\q4_partD.tex", nomtitles replace tex


* F) Estimate mean and distribution of Y_0 and Y_1 for compliers

*To estimate E[Y_0|C=c], see slide 20 of week3.pdf
gen d2 = 1-d
gen outcome2 = lnw * d2
ivregress 2sls outcome2 (d2=z)
matrix partF = e(b)[1,1] \ .

*To estimate E[Y_1|C=c], see slide 20 of week3.pdf
gen outcome1 = lnw * d
ivregress 2sls outcome1 (d=z)
matrix partF[2,1] = e(b)[1,1]

*Save table with results
matrix rownames partF = "$ Y\_0 $ for compliers" "$ Y\_1 $ for compliers"
matrix colnames partF = Mean
esttab matrix(partF, fmt(%04.3f)) using "..\tables\q4_partF.tex", substitute(\_ _) nomtitles replace tex

*Ok but we want full distribution
tab z d, matcell(table)

local p_a = table[1,2] / (table[1,2]+table[1,1])
local p_n = table[2,1] / (table[2,1]+table[2,2])
local p_c = 1 - `p_a' - `p_n'

gen a = z==0 & d==1
la var a "Always Taker - Observed"
gen n = z==1 & d==0
la var n "Never Taker - Observed"

kdensity lnw if a==1, at(lnw) generate(den_a)
la var den_a "Y{sub:1}, Always-takers"
kdensity lnw if n==1, at(lnw) generate(den_n)
la var den_n "Y{sub:0}, Never-takers"

kdensity lnw if d==0 & z==0, at(lnw) generate(den_00)
kdensity lnw if d==1 & z==1, at(lnw) generate(den_11)

gen den_c1 = den_11 * ((`p_c'+`p_a')/`p_c') - den_a * (`p_a'/`p_c')
la var den_c1 "Y{sub:1}, Compliers"
line den_c1 lnw, sort xtitle("Log Wage") ytitle("Density") title("Counterfactual Density of Y{sub:1} for Compliers")
graph export "..\figures\q4_partF_compliers0.pdf", replace

gen den_c0 = den_00 * ((`p_c'+`p_n')/`p_c') - den_n * (`p_n'/`p_c')
la var den_c0 "Y{sub:0}, Compliers"
twoway line den_c0 lnw, sort  xtitle("Log Wage") ytitle("Density") title("Counterfactual Density of Y{sub:0} for Compliers")
graph export "..\figures\q4_partF_compliers1.pdf", replace

twoway line den_c0 den_c1 lnw, sort  xtitle("Log Wage") ytitle("Density") title("Counterfactual Density for Compliers") legend(on  position(6) col(2))
graph export "..\figures\q4_partF_compliers.pdf", replace

* G) Y_0 and Y_1 for always- and never-takers

*Estimate means
qui summ lnw if a==1
matrix partG = `r(mean)' \ .
qui summ lnw if n==1
matrix partG[2,1] = `r(mean)'

*Save table with results
matrix rownames partG = "$ Y\_0 $ for never-takers" "$ Y\_1 $ for always_takers"
matrix colnames partG = Mean
esttab matrix(partG, fmt(%04.3f)) using "..\tables\q4_partG.tex", substitute(\_ _) nomtitles replace tex

twoway line den_a den_n lnw, sort  xtitle("Log Wage") ytitle("Density") title("Observed Density for Always- and Never-takers") legend(on  position(6) col(2))
graph export "..\figures\q4_partG_observed.pdf", replace


* H) Estimation two ways
eststo clear

*First way
egen dummy=group(lotcateg year)
levelsof dummy

levelsof year, local(yearlist)
levelsof lotcateg, local(lotlist)

foreach y in `yearlist' {
	foreach l in `lotlist' {
		qui ivregress 2sls lnw (d = z) if year == `y' & lotcateg == `l'
		matrix partH = nullmat(partH) \ `y', `l', e(b)[1,1]
	}
}
frmttable using "..\tables\q4_partH.tex", statmat(partH) sfmt(g,g,f) tex norowtitl  replace 

*Second way
ivregress 2sls lnw lotcateg#year (lotcateg#year)*z (d = z)
