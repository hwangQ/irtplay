
<!-- README.md is generated from README.Rmd. Please edit that file -->

# irtplay\_1.6.4 (2022-03-29)

o added a beta version of new function, ‘catsib()’, for computing CATSIB
statistic (Nandakumar & Roussos, 2004) to detect DIF on items in CAT.

o added a new function, ‘est\_mg()’, for the multiple group (MG) item
calibration (Bock & Zimowski, 1997). The new function also supports the
MG fixed item calibration (MG-FIPC) method (e.g., Kim & Kolen, 2016).

o added a new argument, ‘se’, in the ‘est\_irt()’ function. If ‘se =
FALSE’, the standard errors of the item parameter estimates are not
computed.

o added a new argument, ‘fix.id’, in the ‘est\_irt()’ function. When
implementing the fixed item parameter calibration (FIPC), the fixed
items can be specified by either of the ‘fix.loc’ or ‘fix.id’ arguments.

o updated the `write.flexmirt()` function so that it can also create a
“-prm.txt” file for flexMIRT software with multiple groups. The previous
version only worked for a single group.

o updated the `write.flexmirt()` function so that it can also create a
“-prm.txt” file for flexMIRT software with multiple groups. The previous
version only worked for a single group.

o added a new simulated data sets called `simMG` with three multiple
groups.

# irtplay\_1.6.3 (2021-11-01)

o resolved the issue occurred when fixing the item guessing parameters
to a specific (e.g., 0.1) in the `est_irt()` function.

o updated the `est_irt()` function to estimate the population latent
ability distribution only when all item parameters are fixed in a test
using the fixed item parameter calibration (FIPC).

o fixed the `est_irt()` function so that the log-likelihood, AIC, and
BIC can be computed based on both the fixed- and freely estimated items
when the FIPC is implemented. In the previous version, those valused
were computed based on only freely estimated items.

o added a new argument of ‘item.id’ in the ‘est\_irt()’ and
‘est\_item()’ functions where a user can provide item IDs.

o added a new ‘rdif()’ function which computes RDIF statistics (Lim,
Choe, & Han, 2022; Lim, Choe, Han, Lee, & Hong, 2021) for analyzing DIF.

o updated ‘plot.test.info()’ function so that (a) multiple item
information functions can be displayed in one panel by setting ‘overlap
= TRUE’ and (b) a plot of conditional standard error of estimation at a
test level can be shown by setting ‘csee = TRUE’.

o updated ‘est\_score()’ function to make it return NA values for
examinees who have all missing responses.

o fixed an error of ‘est\_score()’ function which occurs when an
examinee has missing data for all polytomous items or dichotomous items
(thanks to Craig Wells).

o fixed a few minor issues of ‘irtfit()’ function (thanks to Dimitrios
Zacharatos).

o Updated ‘bring.flexmirt()’ function to read the empirical histogram of
population distribution from “-prm.txt” file.

o Updated ‘run\_flexmirt()’ function to run flexMIRT in which version is
&gt;= 3.6.

o Updated ‘est\_item()’ function to produce a variance-covariance matrix
for item parameter estimates.

o Added ‘vcov()’ method for ‘est\_item’ function.

# irtplay\_1.6.2 (2020-12-14)

o Added ‘coef’, ‘logLik’ methods for ‘est\_item’ function, and added
‘coef’, ‘logLik’, and ‘vcov’ methods for ‘est\_irt’ function.

o Updated ‘est\_irt’ and ‘est\_item’ functions so that the argument of
‘verbose’ can suppress all messages of the parameter estimation
progress.

o Updated ‘bring.flexmirt’ function by including a new argument of
“rePrm.gpc”. In the previous version, the nominal model parameters in
the flexMIRT parameter output file are reparameterized into the (G)PCM
parameters when (G)PCM is fit to data by setting ‘rePrm = TRUE’. In the
new version, however, the nominal model parameters are reparameterized
into the (G)PCM slope/difficulty parameters only when ‘rePrm.gpc =
TRUE’.

# irtplay\_1.6.1 (2020-08-13)

o Included a new function of ‘post\_den’ to compute updated prior
(a.k.a. posterior) densities of the latent ability distribution given a
prior ability distribution, item parameters, and item response data.

o Updated ‘est\_score’ function so that the standard errors of ability
estimates can be computed using eigther the observed item information
function or the expected item information (a.k.a. Fisher information)
function when MLE, MLE with fence (MLEF), or MAP scoring method is used.

o Updated ‘est\_irt’ function to compute the loglikelihood-based fit
statistics of Akaike information criterion (AIC) and Bayesian
information criterion (BIC).

o Updated ‘est\_irt’ function to tally the number of freely estimated
parameters taking the mean and variance parameters of the latent ability
distribution into consideration when ‘fipc = TRUE’.

o Updated ‘est\_irt’ function to suppress printing the observed data
log-likelihood after each EM cycle using the argument of ‘verbose’.

o Fixed an error of the ‘est\_irt’ function when only dichotomous items
are used with ‘fipc = TRUE’. In that condition, an error message of
“subscript out of bounds” was returned in the previous version. No error
message is shown in the updated version. (thanks to Ahmet GUVEN)

o Fixed the ‘lwrc’ function so that it can return the probability
results even when only a single theta value is used.

# irtplay\_1.6.0 (2020-07-14)

The package has been updated significantly in this verstion. In this
version, I have:

o Updated ‘est\_score’ function to estimate ability parameters much
faster than the previous version of the function.

o Updated ‘est\_irt’ and ‘est\_item’ functions to estimate item
parameters much faster than the previous version of the functions.

o Updated ‘test.info’ function to compute items infomation and test
information much faster than the previous version of the function.

o Added an option to use a prior distribution of the item difficulty (or
threshold) parameters in ‘est\_irt’, ‘est\_item’, and ‘llike\_item’
functions.

o Solved unstable item parameter estimation of ‘est\_irt’ and
‘est\_item’ functions which occured when the scaling factor of ‘D’ is
other than 1.0 and ‘use.aprior = TRUE’.

o Fixed an error which occured in the function ‘est\_irt’ when the data
set contains missing values and ‘fix.a.1pl = FALSE’.

# irtplay\_1.5.1 (2020-06-16)

o Included ‘summary’ method to summarize the IRT calibration results
from ‘est\_irt’ or ‘est\_item’ objects.

o Included a new function of ‘getirt’ to extract various estimates
results from ‘est\_irt’ or ‘est\_item’ objects.

o Fixed an error which happens when “DRM” is specified in the model name
in the function ‘est\_irt’.

o Included total computation time in the function ‘est\_irt’.

# irtplay\_1.5.0 (2020-04-12)

o Changed the title of ‘irtplay’ package to “Unidimensional Item
Response Theory Modeling”.

o Included a new function of ‘est\_irt’ to fit unidimensional IRT models
to mixture of dichotomous and polytomous item data using the marginal
maximum likelihood estimation with expectation-maximization (MMLE-EM;
Bock & Aitkin, 1981) algorithm.

o Included the fixed item parameter calibration (FIPC; Kim, 2006)
approach, which is one of useful online calibration methods, in the
function ‘est\_irt’.

o Updated the documentation to explain how to implement the new function
‘est\_irt’.

o Included well-known LSAT6 dichotomous response data set from Thissen
(1982).

o Fixed a problem of inaccurate item parameter estimation in the
function ‘est\_item’ when a prior distribution of the slope parameter is
used with a scaling factor other than D = 1.

o Updated the function ‘bring.flexmirt’ to read the item parameters of
the generalized partial credit model when the number of score categories
are two.

o Updated the function ‘est\_score’ to find a smart starting value when
MLE is used. More specifically, the smart starting value is a theta
value where the log-likelihood is the maximum at the highest peak.

# irtplay\_1.4.1 (2020-02-21)

o Included the function ‘run\_flexmirt’ to implement flexMIRT software
(Cai, 2017) through R.

o Applied a prior distribution to the slope parameters of the IRT 1PL
model when the slope parameters are constrained to be equal in the
function of ‘est\_item’.

o Fixed a problem of using staring values to estimate item parameters in
the function of ‘est\_item’.

# irtplay\_1.4.0 (2020-01-23)

o Fixed a non-convergence problem of the maximum likelihood estimation
with fences (MLEF) in the function of ‘est\_score’.

o Updated the description and introduction of the package.

o Updated the documentation to explain how to implement the function
“est\_item” in more detail.

o Updated the README.md file to explain how to implement the function
“est\_item” in more detail.

# irtplay\_1.3.0 (2019-11-17)

o Included the function ‘llike\_score’ to compute the loglikelihood
function of ability for an examinee.

o Updated the function ‘est\_item’ to find better starting values for
item parameter calibration.

o Updated the function ‘est\_item’ to exclude items that contains no
item response data during the item parameter estimation.

o Updated the function ‘est\_item’ to count the number of item responses
for each item used to estimate the item parameters.

o Updated the function ‘est\_score’ to find better starting values when
MLE is used.

o Updated the function ‘est\_score’ to address NaNs of gradient values
and NaNs of hessian values when MLE, MLEF, or MAP is used.

o Fixed a problem of the function ‘est\_score’, which returned an error
message when a vector of an examinee’s response data was used in the
argument of ‘x’.

# irtplay\_1.2.0 (2019-10-16)

o Fixed a problem of the function ‘est\_score’, which returned an error
message when only one dichotomous item or one polytomous item was
included in the item meta data set.

o Fixed a problem of the function ‘est\_item’, which returned an error
message when the inverse of hessian matrix is not obtainable.

o Included the ‘maximum likelihood estimation with fences scoring method
(Han, 2016) in the function ’est\_score’.

o Included the ‘inverse test characteristic curve (TCC)’ scoring method
(e.g., Stocking, 1996) in the function ‘est\_score’.

o Included the function ‘llike\_item’ to compute the loglikelihood
values of items.

# irtplay\_1.1.0 (2019-09-15)

o For the function ‘est\_item’, default parameters of a-parameter prior
distribution were revised

o Updated the function ‘est\_item’ to find better starting values for
item parameter calibration.

o Updated the function ‘est\_score’ to estimate an ability in a brute
force way when MLE or MAP fails to find the solution.

o Updated the function ‘irtfit’ to compute the likelihood ratio
chi-square fit statistic (G2; Mckinley & Mills, 1985).

# irtplay\_1.0.0 (2019-08-21)

o initial release on CRAN
