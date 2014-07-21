#!/bin/bash

if [ $1 == 'weibull' ] 
then
    echo 'Running All Weibull Fit Code'
    Rscript prot_ex_weibull.R cens                                                                                                     
    Rscript prot_ex_weibull.R int
    Rscript prot_ex_weibull.R uncens
    Rscript prot_ex_weibull.R compl   
fi

”if [ $1 == 'piecewise' ]
then
    echo 'Running All Piecewise Fit Code'
    Rscript prot_ex_pwfit.R cens
    Rscript prot_ex_pwfit.R int
    Rscript prot_ex_pwfit.R uncens
    Rscript prot_ex_pwfit.R compl
fi





