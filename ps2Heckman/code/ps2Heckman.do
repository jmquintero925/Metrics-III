****************************************************************************************************
****************************************************************************************************
***********************************    Empirical Analysis III    ***********************************
***********************************    Part B Problem Set 2      ***********************************
***********************************          Question 8          ***********************************
***********************************   Written By: Isaac Norwich  ***********************************
***********************************     Date: May 17th, 2022     ***********************************
****************************************************************************************************
****************************************************************************************************

****************************************************************************************************
* Import, Clean, and Merge
*  NB: Open Q8.zip, copy the CSVs from inside the Q8 folder, and put them in the data folder
****************************************************************************************************
*Setup
clear all
mat drop _all
set more off
cd "C:\Users\inorwich\OneDrive - The University of Chicago\Metrics III\ps2Heckman\code"

/*
*Install useful/necessary packages
net install scheme-modern, ///
  from("https://raw.githubusercontent.com/mdroste/stata-scheme-modern/master/") replace
set scheme modern, perm

ssc install estout, replace
ssc install tabout, replace

*Sometimes this one doesn't work right - it's for frmttable
cap net install sg97_5.pkg, replace
*/

*Load the data
forv n=1/5 {
    import delimited "..\data\econ31200-ps2-dataset-`n'.csv", varnames(1) clear
    
    rename (y d) (y_`n' d_`n')
    
    tempfile temp`n'
    save `temp`n''
}

*Merge the data
use `temp1', clear
forv n = 2/5 {
    merge 1:1 v1 x z using `temp`n'', assert(3) nogen
}
assert v1 == _n
drop v1

save "..\data\merged.dta", replace

****************************************************************************************************
* A) Estimate, for a given X, Pr(D = 1) for each data set and graph the estimate as a function of Z. 
*  What is the subjective treatment effects for each data set? Define the graph for each data set.
****************************************************************************************************
use "..\data\merged.dta", clear

forv n = 1/5 {
    probit d_`n' x, nocons
    predict prob_d_`n'

    tw scatter prob_d_`n' z  || lfit prob_d_`n' z, ylabel(0.45(.025).55, format(%4.3f)) title("Results from Probit of D on X for Dataset `n'") legend(off)
    graph export "..\figures\q8_parta_d`n'.png", replace
}

****************************************************************************************************
* B) Heckman two-step procedure
****************************************************************************************************
use "..\data\merged.dta", clear

rename z w

forv n = 1/5 {
    qui probit d_`n' x w, nocon
    local gamma_x = e(b)[1,1]
    local gamma_w = e(b)[1,2]
    predict p_score
    la var p_score "Propensity Score"
    gen mu_d_z = invnormal(p_score)
    gen inv_mill_tilde = - normalden(mu_d_z) /      normal(mu_d_z)
    gen inv_mill       =   normalden(mu_d_z) / (1 - normal(mu_d_z))

    qui reg y_`n' x inv_mill_tilde if d_`n' == 1, nocons
    local beta_1 = e(b)[1,1]
    local rho_1  = e(b)[1,2]

    qui reg y_`n' x inv_mill       if d_`n' == 0, nocons
    local beta_0 = e(b)[1,1]
    local rho_0  = e(b)[1,2]

    matrix row = `gamma_x', `gamma_w', `beta_1', `beta_0', `rho_1', `rho_0'
    matrix ests = nullmat(ests) \ row

    gen TE = (`beta_1'-`beta_0')*x
    qui summ TE
    local ATE = `r(mean)'
    tw scatter TE p_score || lfit TE p_score, ytitle("TE") title("Dataset `n'") legend(off)
    graph export "..\figures\q8_partb_d`n'_te.png", replace

    gen TT = TE + (`rho_1'-`rho_0')*inv_mill_tilde
    qui summ TT
    local ATT = `r(mean)'
    tw scatter TT p_score || lfit TT p_score, ytitle("TT") title("Dataset `n'") legend(off)
    graph export "..\figures\q8_partb_d`n'_tt.png", replace

    gen TU = TE + (`rho_1'-`rho_0')*inv_mill
    qui summ TU
    local ATUT = `r(mean)'
    tw scatter TU p_score || lfit TU p_score, ytitle("TU") title("Dataset `n'") legend(off)
    graph export "..\figures\q8_partb_d`n'_tu.png", replace

    *PRTE with z=0.5 to zprime=1, -0.5
    gen mu_z0 = TE - (`gamma_x'*x + `gamma_w'*0.5)
    local count = 1
    foreach z in 1 -0.5 {
        gen mu_z1 = TE - (`gamma_x'*x + `gamma_w'*`z')
        gen PRTE`count' = TE*(normal(mu_z1)-normal(mu_z0)) - (`rho_1'-`rho_0')*(normalden(mu_z1)-normalden(mu_z0))

        qui summ PRTE`count'
        local PRTE`count' = `r(mean)'
        
        local count = `count' + 1
        drop mu_z1 
    }
    tw scatter PRTE1 PRTE2 p_score || lfit PRTE1 p_score || lfit PRTE2 p_score, ytitle("PRTE") title("Dataset `n'") legend(off)
    graph export "..\figures\q8_partb_d`n'_prte.png", replace

    *MTE
    *Note this v is equal to mu_d_z from above
    gen v = (`gamma_x'*x + `gamma_w'*w)
    la var v "(Gamma_x*X + Gamma_w*W)"
    *gen v = -`gamma_w'*w
    *la var v "-Gamma_w*W"
    gen MTE = TE + (`rho_1'-`rho_0')*v
    la var MTE "MTE"
    tw scatter MTE p_score || lfit MTE p_score, ytitle("MTE") title("Dataset `n'") legend(off)
    graph export "..\figures\q8_partb_d`n'_mte.png", replace

    *Graph All
    tw scatter TE TT TU PRTE1 PRTE2 MTE p_score, title("Dataset `n'") xlabel(0.35(.05).65, format(%4.3f))
    graph export "..\figures\q8_partb_d`n'_all.png", replace

    *Save results
    matrix results = nullmat(results) \ `ATE', `ATT', `ATUT', `PRTE1', `PRTE2'

    rename (mu_d_z p_score v TE TT TU PRTE1 PRTE2 MTE) (mu_d_z_`n' p_score_`n' v_`n' TE_`n' TT_`n' TU_`n' PRTE1_`n' PRTE2_`n' MTE_`n')

    *cap drop p* 
    cap drop mu_z0
    cap drop inv_mill*
    *cap drop TE TT TU PRTE* MTE

}
mat li ests
mat li results

frmttable using "..\tables\q8_partb.tex", statmat(results) ctitles("Dataset","ATE","ATT","ATUT","PRTE1","PRTE2","MTE") rtitles("1"\"2"\"3"\"4"\"5") fragment tex replace

gen n = _n
reshape long y_ d_ mu_d_z_ p_score_ v_ TE_ TT_ TU_ PRTE1_ PRTE2_ MTE_, i(n) j(dataset)
rename *_ *
save "..\data\final.dta", replace