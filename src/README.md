## Package requirement

To run the codes,  please ensure you have installed the following R packages.

| Package Name | Description               |
| ------------ | ------------------------- |
| quantreg     | for quantile regression   |
| lpSolve      | for linear programming    |
| parallel     | for detectCores()         |
| snowfall     | for parallel programming  |
| LowRankQP    | for quadratic programming |
| expectreg    | to find expectile         |
| purrr        | for cross3()              |
| quadprog     | for quadratic programming |



## Main code

These codes conduct the simulations and two empirical applications.

| File Name              | Description                                                  |
| ---------------------- | ------------------------------------------------------------ |
| functions.R            | Functions used below                                         |
| DGP1.R                 | Perform the simulations of DGP 1 in Simulation Design I      |
| DGP2.R                 | Perform the simulations of DGP 2 in Simulation Design I      |
| DGP3.R                 | Perform the simulations of DGP 3 in Simulation Design I      |
| design2.R              | Perform Simulation Design II                                 |
| data_prep_collection.R | Prepare the data for the prediction of the collection volume of shipping carrier stores |
| main_collection.R      | Perform the prediction of the collection volume of shipping carrier stores |
| data_prep_stock.R      | Prepare the data for the forecast of excess stock returns    |
| main_stock_rolling.R   | Perform the forecast of excess stock returns                 |