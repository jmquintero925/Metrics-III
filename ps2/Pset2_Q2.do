********************************************************************
********************************************************************
*********	 Empirical Analysis III 		************************
*********	 Problem Set 2 					************************
********* 	 Question 2 					************************
********************************************************************
********************************************************************

* Housekeeping
clear all

* Read data frane
use "Data/lalonde2.dta", clear 

* Covariate set used throughout
global controls age educ black married nodegree hisp kids18 kidmiss re74

* (a)
* Create mean difference table for randomization
eststo D1: quietly estpost summarize $controls if (treated == 1 & sample==1)
eststo D0: quietly estpost summarize $controls if (treated == 0 & sample==1)
eststo mdiff: quietly estpost ttest $controls if sample==1, by(treated) 
esttab D1 D0 mdiff using Tables/ps2q2atab1.tex , cells("mean(pattern(1 1 0) fmt(3)) sd(pattern(1 1 0)) b(star pattern(0 0 1) fmt(3)) t(pattern(0 0 1) par fmt(3))") label  mtitles("Treated" "Untreated" "Mean Difference") tex replace

* (b)
* Calculate the mean difference in the randomize sample
eststo Y1s1: quietly estpost summarize re78 if (treated == 1 & sample==1)
eststo Y0s1: quietly estpost summarize re78 if (treated == 0 & sample==1)
eststo ateS1: quietly estpost ttest re78 if sample==1, by(treated) 

*Export results
esttab Y1s1 Y0s1 ateS1 , cells("mean(pattern(1 1 0) fmt(3)) sd(pattern(1 1 0)) b(star pattern(0 0 1) fmt(3)) t(pattern(0 0 1) par fmt(3))") label  mtitles("Treated" "Untreated" "ATE")
esttab Y1s1 Y0s1 ateS1 using Tables/ps2q2btab1.tex , cells("mean(pattern(1 1 0) fmt(3)) sd(pattern(1 1 0)) b(star pattern(0 0 1) fmt(3)) t(pattern(0 0 1) par fmt(3))") label  mtitles("Treated" "Untreated" "ATE") tex replace 

* (c)
* Create new sample for OLS
gen samplec = (treated == 1 & sample==1) | (sample==2)
* Create treated variable for new sample 
gen D = .
replace D=1 if treated==1
replace D=0 if (sample==2) 
* Run OLS 
reg re78 D 
outreg2 using Tables/ps2q2ctab1, nocons dec(3) addtext("Controls", "No") replace 

reg re78 D $controls
outreg2 using Tables/ps2q2ctab1, nocons dec(3) addtext("Controls","Yes") append drop($controls)

* (d)
* Balancing in covariates
eststo D1: quietly estpost summarize $controls if (D == 1)
eststo D0: quietly estpost summarize $controls if (D == 0)
eststo mdiff: quietly estpost ttest $controls if D!=., by(D) 
esttab D1 D0 mdiff using Tables/ps2q2dtab1.tex , cells("mean(pattern(1 1 0) fmt(3)) sd(pattern(1 1 0)) b(star pattern(0 0 1) fmt(3)) t(pattern(0 0 1) par fmt(3))") label  mtitles("Treated" "Untreated" "Mean Difference") tex replace

* Histogram for age
twoway (histogram age if D==0, fcolor(cranberry%70) lcolor(maroon) width(1.8) start(16)) /// 
	   (histogram age if D==1, fcolor(ebblue%70) lcolor(navy) width(1.8) start(16)),  graphregion(fcolor(white)) ///
	   legend(order(1 "Untreated" 2 "Treated") cols(1) region(lcolor(white)) ring(0) pos(2)) ylabel(,nogrid format(%3.2f)) ///
	   xti("Age")
graph export "Figures/ps2_q2d_age.pdf", replace

* Histogram for educ
twoway (histogram educ if D==0, fcolor(cranberry%70) lcolor(maroon) width(1) start(0)) /// 
	   (histogram educ if D==1, fcolor(ebblue%70) lcolor(navy) width(1) start(0)),  graphregion(fcolor(white)) ///
	   legend(order(1 "Untreated" 2 "Treated") cols(1) region(lcolor(white)) ring(0) pos(2)) ylabel(,nogrid format(%3.2f)) ///
	   xti("Years of Schooling") 
graph export "Figures/ps2_q2d_educ.pdf", replace

* Histogram for kids18
twoway (histogram kids18 if D==0, fcolor(cranberry%70) lcolor(maroon) width(1) start(0)) /// 
	   (histogram kids18 if D==1, fcolor(ebblue%70) lcolor(navy) width(1) start(0)),  graphregion(fcolor(white)) ///
	   legend(order(1 "Untreated" 2 "Treated") cols(1) region(lcolor(white)) ring(0) pos(2)) ylabel(,nogrid format(%3.2f)) ///
	   xti("# Kids younger than 18")
graph export "Figures/ps2_q2d_kids18.pdf", replace

* Histogram for re74
twoway (histogram re74 if D==0, fcolor(cranberry%70) lcolor(maroon) width(2000) start(0) frac) /// 
	   (histogram re74 if D==1, fcolor(ebblue%70) lcolor(navy) width(2000) start(0) frac),  graphregion(fcolor(white)) ///
	   legend(order(1 "Untreated" 2 "Treated") cols(1) region(lcolor(white)) ring(0) pos(2)) ylabel(,nogrid format(%3.2f)) ///
	   xti("Real Earnings 1974")
graph export "Figures/ps2_q2d_re74.pdf", replace


* Run Psmatch
probit D $controls
predict _pscore
gen _pscore2 = _pscore 
replace _pscore2 = . if _pscore < 0.05 | _pscore > 0.95
* Histogram for re74
twoway (histogram _pscore2 if D==0, fcolor(cranberry%70) lcolor(maroon) width(0.05) start(0.05)) /// 
	   (histogram _pscore2 if D==1, fcolor(ebblue%70) lcolor(navy) width(0.05) start(0.05)),  graphregion(fcolor(white)) ///
	   legend(order(1 "Untreated" 2 "Treated") cols(1) region(lcolor(white)) ring(0) pos(2)) ylabel(,nogrid format(%3.2f)) ///
	   xti("Propensity Score")
graph export "Figures/ps2_q2d_pscore.pdf", replace

* Calculate 
psmatch2 D, out(re78) n(1) pscore(_pscore2)
matrix A = (r(att))

reg re78 D $controls [pweight = _weight]
matrix A = (A, _b[D])

esttab matrix(A)


qui lpoly D _pscore2, kernel(rec) deg(1) nodraw
* estimate matching using LLR and this bandwidth
psmatch2 D, outcome(re78) llr qui bw(`r(bwidth)') pscore(_pscore2) 






