TODO:
1. set hyper-parameters for
    1. change `X[100]` to `X[200]`
        - use `X[200]` gives unusually large negative log-likelihood with table4's parameter
    2. drop `ln(deriv)`
2. organize MLE results into table
    1. for each country: `[batch_num]_[algo_name]`: estimated parameters; `X[200]`: T/F; `ln(derive)`: T/F; `time_lapse`: run_time
3. Plot neg log vs. X[i], i:0-199 using table 4 parameter
4. tweal tolerance level for both `SLSQP` and `trust-constr`

Data:
1. 1Y,3Y,5Y CDS data are not of same # of observations for added countries

1. discount fxn a bit different
2. truncate to have same start date as in PSAEJ for CDS data?
    - tweak code to allow variable length CDS data
    - careful matching D fxn
3. find macro variables for risk recomposition

<!-- 3. numba the code -->

### date log:
1. AEJ
    1. 2021_8_4_19_20: ALGOs: 'SLSQP':dict(maxiter=5000, disp=True, ftol=1e-16), 'trust-constr': dict(maxiter=5000, verbose=2, xtol=1e-16, gtol=1e-16)
    2. 2021_8_10_23_30: ALGOs: 'SLSQP':dict(maxiter=5000, disp=True, ftol=1e-8), 'trust-constr': dict(maxiter=5000, verbose=2, xtol=1e-8, gtol=1e-8)
        - SLSQP terminates prematurely
2. EXT
    1. "2021_8_20": starting value being XGUESSES[-1]



### Batch Hyperparameter
See input/data_v/batch_*.csv





XNYU901
PICKLE