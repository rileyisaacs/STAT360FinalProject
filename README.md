# STAT360FinalProject
Group members: Darcie Reiter, Riley Isaacs, Carter Phillips 

Files in project folder

#### folder: data  
Contains data used in manual tests and examples.  
- marstestdata.rda  
  - contains data for test 1
- laptop_price.csv  
  - contains data for test 2
- Real estate.csv  
  - contains data for test 3

#### folder: data-raw  
Contains data used in other tests.  
- DATASET
- testbwd_stepwise (R workspace)  
  - tests the backward stepwise portion of the MARS algorithm
- testfwd_stepwise (R workspace)  
  - tests the forward stepwise portion of the MARS algorithm
- testmc (R workspace)  
  - tests the mars.control object

#### folder: R
Defines methods used for MARS.  
- mars.R
  - defines the mars method, including its helper functions and the mars.control object
- anova.mars.R
  - defines the anova.mars method, which returns a one-way ANOVA table with F-statistics and tests for $H_0 : \beta_i = 0$
- plot.mars.R
  - defines the plot.mars method, which implements the plot.lm method for MARS objects  
- predict.mars.R  
  - defines the predict.mars method, which returns the predicted values of the model on either a new dataset or the data contained within the MARS object
- print.mars.R  
  - defines the print.mars method, which prints the initial mars method call and the coefficients of the model
- summary.mars.R
  - defines the summary.mars method, which returns a detailed summary of the MARS object, including:
    - initial mars method call
    - information on data contained within the MARS object
    - residuals of prediction
    - prediction coefficients and t-tests
    - standard deviation of fitted model
    - degrees of freedom associated with each coefficion
    - $R^2$ and adjusted $R^2$ of the model
    - information on F-statistic of model
    - covariance matrix

#### folder: tests
Allows for testing of MARS package.  
- testthat.R 
  - sets up automatic testing of MARS package
- test.R
  - contains manual tests and examples

#### folder: testthat
Contains files used for automatic testing of the mars algorithm.
- testbwd_stepwise (R workspace)
- testbwd_stepwise.R
- testfwd_stepwise (R workspace)
- testfwd_stepwise.R
- testmc (R workspace)
- testmc.R
- testmars (R workspace) 
- testmars.R
- testpredict (R workspace)
- testpredict.R

DESCRIPTION  
- contains description of the package
LICENSE  
- contains license of the package
NAMESPACE

Contributions

Darcie: helped set up main mars.R folder and its functions, set up all the tests in the testthat folder, completed print.mars.R, plot.mars.R. and anova.mars.R, found datasets and completed test.R file, did part of documentation, set up README.md file

Riley: helped set up the main mars.R folder and its functions, set up the github, completed summary.mars, helped write documentation 

Carter: wrote documentation and description, altered anova.mars and plot.mars, wrote README.md file
