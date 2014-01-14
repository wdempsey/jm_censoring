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

The User is required to define several functions in order to fit the model.

* Mean Model :  We require two functions **X1** and **X2**.  Each takes a subset of Table_2 for a specific patient, and produces the covariates to be used in the mean model.  **X2** will also takes the survival time, as an input, and is used to estimate covariates which are dependent on the survival time. A simple example would be where **X2** is just the revival function:

<pre><code> X2 = function(t, table) {
return(revival = table2$obs_times - as.numeric(lapply(table$id, survival)))
}
</code></pre>

where 'survival' is a function that returns the survival of the patient in question.

* Covariance Model : We require a function **SigmaCalc** which computes the covariance matrix given a set of parameters, **covparams**.  An example would be iid patients where there is a exponential covariance term depending on the observation times for a particular patient.  Such an example would give:

<pre><code> SigmaCalc = function(covparams, table) {
sigmasq0 = covparams[1]
sigmasq1 = covparams[2]
lambda = covparams[3]
Cov = sigmasq0 * diag(length(table)) + sigmasq1* exp(-abs(outer(table$obstimes, table$obstime)))
return(Cov)
}
</code></pre>

* Initial Values : The User is required to provide initial values for the parameter estimates.  These typically are the result of analysis done on only the uncensored individuals.


# Sample Fit

We provide a set of simulated [examples](./examples.R) to test the code.  In each, the model is fit by first sourcing the [fitting code](./MLE_censoring.R), and then running:

<pre><code> revival_model(Table1, Table2, X1, X2, Sigmacalc, meanparams, covparams, theta)
</code></pre>

where **theta** denotes the parameters associated with the survival model.  We also provide a [case study](./prothrombin_example.R).
