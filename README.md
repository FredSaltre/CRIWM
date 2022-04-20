# CRIWM

Code associated to the paper "Estimating extinction time using radiocarbon dates" by Salvador Herrando Perez and Frederik Saltre

## Abstract:

1.	The extinction of a species is a key demographic event often signalling major climatic, ecological and/or evolutionary shifts that can be investigated using the fossil record. In that context, radiocarbon (14C) dating has become a popular tool to time and test for scenarios of extinction that can inform on how species respond to past and future ecosystem changes.
2.	We develop CRIWM, a method for estimating extinction (and arrival) time from time series of 14C dates while accounting for probability density functions (PDF) deriving from 14C calibration. The sister method GRIWM assumes normal chronometric errors and is inappropriate for 14C chronologies as PDFs are often non-Gaussian and multi-modal. Compared to GRIWM, CRIWM reduces by 4-fold the gap between true and estimated extinction times and the width of the confidence intervals, and is less sensitive to sample size, dating errors and the temporal distribution of fossils in a species’ fossil record.
3.	We build the R package Rextinct with three user-friendly functions for computing CRIWM and GRIWM, and the PDF moments, modes and intercepts of 14C calibrations. CRIWM and GRIWM accept time series comprising only 14C dates or observations (fossils, sightings) dated by multiple chronometric methods, and calibration curves for fossils sampled in the Northern and Southern Hemispheres and marine environments. For both methods, we implement two different estimators of extinction time whether they are reliant on a p-values (original conceptualization) or not (novel).
4.	CRIWM and Rextinct are robust tools to infer time extinction and arrival events using 14C chronologies spanning the entire Holocene to 50,000 calendar years Before Present, and can be used to investigate demographic phenomena queried through the fossil record such as migrations, domestications, extirpations and rewilding.


## Folders:

- Rextinct: contains the R package <em>Rextinct</em> (with three user-friendly functions) designed  for computing CRIWM (error-distribution free), and its sister method GRIWM (normal dating errors), and for calibrating 14C dates, respectively. The package accepts observations (fossils, sightings) dated by any combination of other chronometric methods. For 14C chronologies, Rextinct handles fossil dates from the Northern and Southern Hemispheres and marine environments as they require specific 14C calibrations. 

- Sensitivity: R-script to examine the robustness of CRIWM and GRIWM. The script runs a sensitivity analysis in three stages by: <em>(i)</em> simulating times series of 14C dates for different sample sizes and time periods, <em>(ii)</em> comparing accuracy and precision of estimating extinction time by both methods, and <em>(iii)</em> evaluating the sensitivity (based on Hilbert-Schmidt Independence Criterion (HSIC, Gretton et al. 2005) of such accuracy and precision to the statistical properties of the simulated time series. Accuracy = the difference between the simulated true extinction time and CRIWM’s (or GRIWM’s) estimate of extinction time (the lower the difference, the higher the accuracy). Precision = the width of the confidence interval around CRIWM’s (or GRIWM’s) estimate (the narrower the confidence interval, the higher the precision).


## Data:

<em>moa210_holdaway.txt</em>: moa 14C dates spanning 564 to 5,503 14C years BP obtained from 270 fossils collected on the South Island of New Zealand (Holdaway et al., 2014). 

Reference:  Holdaway, R.N., Allentoft, M.E., Jacomb, C., Oskam, C.L., Beavan, N.R. & Bunce, M. (2014) An extremely low-density human population exterminated New Zealand moa. Nature Communications, 5, 5436. https://doi.org/10.1038/ncomms6436


*************************
Frederik Saltre, Flinders University, frederik.saltre@flinders.edu.au April 2021

![Screen Shot 2022-04-19 at 4 08 06 pm](https://user-images.githubusercontent.com/46954120/163941223-e32cc3ce-5562-4de3-af6a-165c5adf8a15.png)
