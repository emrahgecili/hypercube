# hypercube
This repository is for code of hypercube procedure described in paper: 'Hypercubes to identify key environmental exposures predictive of rapid cystic fibrosis lung disease progression'.

The hypercube procedure includes three steps:

1. Reduction.Phase

2. Exploratory.Phase

3. ModelSelection.Phase

Check R code file in each subfolder for details.

The example input data 'input_data.RData' in folder 'Reduction.Phase' is generated from simulation study, with settings:

1. Includes 152 patients, with patients age from 6 to 20 years.
2. The number of FEV1 measurements for each patient varies between 1 and 87.
3. The number of total predictors is 50, with number of signal variables is 5.
4. The 50 predictors are generated from a p-dimensional multivariate normal distribution with mean 0 and compound symmetric covariance matrix, where the diagonal
elements are 1 and off-diagonal elements are 0.1.
5. The mean effect is fixed at 90, time effect is fixed at 0.15. Coefficient for signal variables is 5, and for noise variables is 0. So the coefficient setting for all predictors is {90, 0.15, 5, 5, 5, 5, 5, 0, 0, 0, ..., 0}, with first 5 variables are signal variables ('P1' ~ 'P5').
6. Variance terms sigma^2 = 6, omega^2 = 125, tau^2 = 77.
