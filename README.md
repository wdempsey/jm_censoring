Revival Model : Estimation under Censoring
============

The code fits a specified revival model according to the joint likelihood of the health process and survival time.  The revival model is defined [here](http://www.stat.uchicago.edu/~pmcc/reports/revival.pdf).

# Data Format

We require two separate tables for the data.  **Table 1** contains all the Survival Data.  The first column is an 'id' associated to each patient, followed by their survival or censoring time, and finally an indicator of whether the patient was censored or not.  Additional columns for covariates may be included.  A typical subset of Table 1 would be of the form:


| Id        | Survival/Censoring Time | Censored  |
| ------------- |:-------------:| -----:|
| 1     | 2.25  | 0 |
| 2     | 0.06      |   1 |
| 3     | 7.51      |    0 |


**Table 2** contains all the observations of the health process for each patient, as well as any covariates that may be used in modeling the revival process.  The first column is an 'id' which must match Table 1 'id' column.  Column 2 are the observed values of the health process, and column 3 the observation times as measured since recruitment.  A typical subset of Table 2 would be:

| Id        | Y | Observation Time  |
| ------------- |:-------------:| -----:|
| 1     | 10.20 |  0      |  
| 1     | 6.10  |  1      |  
| 2     |  5.22 |  0      |
| 3     | 13.10 | 0      |
| 3     |11.44  |   1      |  
| 3     |  10.21 | 2      |

# User Determined Functions

# Sample Fit
