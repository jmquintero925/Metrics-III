****************************************************************************************************
****************************************************************************************************
***********************************    Empirical Analysis III    ***********************************
***********************************        Problem Set 3         ***********************************
***********************************          Question 4          ***********************************
***********************************   Written By: Isaac Norwich  ***********************************
***********************************    Date: April 18th, 2022    ***********************************
****************************************************************************************************
****************************************************************************************************

*Setup
clear all
mat drop _all
set more off
cd "C:\Users\inorwich\OneDrive - The University of Chicago\ps3\code"

/*
*Install useful/necessary packages
net install scheme-modern, ///
  from("https://raw.githubusercontent.com/mdroste/stata-scheme-modern/master/") replace
set scheme modern

ssc install estout, replace
ssc install tabout, replace

*Sometimes this one doesn't work right - it's for frmttable
cap net install sg97_5.pkg, replace
*/

*Load the data
use "..\data\lottery.dta", clear

****************************************************************************************************
* B) Assess instrument relevance
****************************************************************************************************
*Regress D on Z and test if the coefficient on Z is 0 or not
eststo: reg d z
esttab using "..\tables\q4_partB.tex", se nostar label nonum tex replace

****************************************************************************************************
* C) Estimate the return to medical school on earnings in 2007 using IV.
****************************************************************************************************
*TSLS
eststo clear
eststo: ivregress 2sls lnw (d = z)

*"By hand"
reg d z
predict d_hat
la var d_hat "$\widehat{D}$ from D on Z"
eststo: reg lnw d_hat

*Output both results
esttab using "..\tables\q4_partC.tex", se nostar label nonum ///
  mtitles("2SLS manually" "2SLS using ivregress") tex replace

****************************************************************************************************
* D) Count the number of compliers. Compare to number of applicants by gender.
****************************************************************************************************
*First way: MHE pg 167 E[D|Z=1]-E[D|Z=0]
qui summ d if z == 1
local d_z1 = `r(mean)'
qui summ d if z == 0
local d_z0 = `r(mean)'
di "The proportion of compliers in the sample is " %04.3f `d_z1'-`d_z0'
di "The number of compliers in the sample is " %7.3f (`d_z1'-`d_z0')*_N

*Second way: 
foreach type in 0 1 2  {
	*Note that female != 2 is the total population! Tricky tricky
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

****************************************************************************************************
* F) Estimate mean and distribution of Y_0 and Y_1 for compliers
****************************************************************************************************
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
esttab matrix(partF, fmt(%04.3f)) using "..\tables\q4_partF.tex", ///
  substitute(\_ _) nomtitles replace tex

*Ok but we want full distribution, not just the means. Follow the appendix from week3.pdf.
tab z d, matcell(table)

*Generate probabilities for each group
local p_a = table[1,2] / (table[1,2]+table[1,1])
local p_n = table[2,1] / (table[2,1]+table[2,2])
local p_c = 1 - `p_a' - `p_n'

gen a = z==0 & d==1
la var a "Always Taker - Observed"
gen n = z==1 & d==0
la var n "Never Taker - Observed"

*Generate densities for always- and never-takers
kdensity lnw if a==1, at(lnw) generate(den_a)
la var den_a "Y{sub:1}, Always-takers"
kdensity lnw if n==1, at(lnw) generate(den_n)
la var den_n "Y{sub:0}, Never-takers"

*Generate other observed densities
kdensity lnw if d==0 & z==0, at(lnw) generate(den_00)
kdensity lnw if d==1 & z==1, at(lnw) generate(den_11)

*Generate counterfactual distribution for c1 compliers and plot
gen den_c1 = den_11 * ((`p_c'+`p_a')/`p_c') - den_a * (`p_a'/`p_c')
la var den_c1 "Y{sub:1}, Compliers"
tw line den_c1 lnw, sort ///
  xtitle("Log Wage") ytitle("Density") title("Counterfactual Density of Y{sub:1} for Compliers")
graph export "..\figures\q4_partF_compliers0.pdf", replace

*Generate counterfactual distribution for c0 compliers and plot
gen den_c0 = den_00 * ((`p_c'+`p_n')/`p_c') - den_n * (`p_n'/`p_c')
la var den_c0 "Y{sub:0}, Compliers"
tw line den_c0 lnw, sort ///
  xtitle("Log Wage") ytitle("Density") title("Counterfactual Density of Y{sub:0} for Compliers")
graph export "..\figures\q4_partF_compliers1.pdf", replace

*Plot both of them together
two line den_c0 den_c1 lnw, sort                                                     ///
  xtitle("Log Wage") ytitle("Density") title("Counterfactual Density for Compliers") ///
  legend(on  position(6) col(2))
graph export "..\figures\q4_partF_compliers.pdf", replace

****************************************************************************************************
* G) Calculate mean and distribution of Y_0 and Y_1 for always- and never-takers
****************************************************************************************************
*Estimate means
qui summ lnw if a==1
matrix partG = `r(mean)' \ .
qui summ lnw if n==1
matrix partG[2,1] = `r(mean)'

*Save table with results
matrix rownames partG = "$ Y\_0 $ for never-takers" "$ Y\_1 $ for always_takers"
matrix colnames partG = Mean
esttab matrix(partG, fmt(%04.3f)) using "..\tables\q4_partG.tex", ///
  substitute(\_ _) nomtitles replace tex

*Plot the distributions calcuted above
twoway line den_a den_n lnw, sort                                                             ///
  xtitle("Log Wage") ytitle("Density") title("Observed Density for Always- and Never-takers") ///
  legend(on  position(6) col(2))
graph export "..\figures\q4_partG_observed.pdf", replace

****************************************************************************************************
* H) LATE estimation two ways
****************************************************************************************************
*First way: X-specific betas and weights
levelsof year, local(yearlist)
levelsof lotcateg, local(lotlist)

foreach y in `yearlist' {
	foreach l in `lotlist' {
		*Estimate betas
		qui ivregress 2sls lnw (d = z) if year == `y' & lotcateg == `l'
		local b = e(b)[1,1]

		*Weights
		qui reg d z#lotcateg#year lotcateg#year
		local pi = _b[1.z#`l'.lotcateg#`y'.year]
		qui summ z if year == `y' & lotcateg == `l'
		local num = `pi'^2 * `r(Var)'
		local share = `r(N)' / _N
		matrix partH = nullmat(partH) \ `y', `l', `b', `num', `share'
	}
}
matrix colnames partH = "year" "lotcateg" "b_x" "num_x" "share_x"
tempfile data
save `data', replace

clear
svmat partH, names(col)

*Construct denominator for weight, check that shares sum to 1
qui summ share_x
assert float(`r(sum)') == 1
gen denom_x = num_x*share_x
qui summ denom_x
local denom = `r(sum)'

*Construct weight for each x
gen wt_x = num_x / `denom'

*Note: beta=E[beta(x)*omega(x)] is the sum over x of beta(x)*omega(x) TIMES Pr(X=x)
gen late_x = share_x*b_x*wt_x
summ late_x
local late_wt: di %06.5f `r(sum)'
di "The weighted LATE is estimated as " `late_wt'

*Export results to table
keep year lotcateg b_x wt_x share_x
mkmat _all, matrix(partH)
matrix colnames partH = "Year" "Lottery Category" "$\beta(x)$" "Proportion" "Weight"
frmttable using "..\tables\q4_partH_tab1.tex",                  ///
  note("$\beta=E[\beta(x)\omega(x)]$ is calculated as `late_wt'") ///
  statmat(partH) sfmt(g,g,f,f,f) norowtitl fragment tex replace

*Second way: control for lotcateg x year dummies and interact Z with these dummies
use `data', clear
eststo clear
eststo: ivregress 2sls lnw lotcateg#year (d = z#lotcateg#year)

la var year "Year"
esttab using "..\tables\q4_partH_tab2.tex", b(%6.5f) se nostar nobaselevels label ///
  mtitle("Log Wage") nonum tex replace


