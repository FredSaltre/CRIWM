#######################################################################################
## package: Rextinct v.1
## function: Rextinct::griwm by Salvador Herrando-Perez (salherra@gmail.com)
## Estimation of extinction time using times series of fossil dates
## Date = 27 March 2021
##
## Based on original code by Corey Bradshaw (cjabradshaw@gmail.com) and Frederik Saltre (frederik.saltre@gmail.com)
## Bradshaw et al 2012, Quaternary Science Reviews 33: 14-19 / http://doi.org/10.1016/j.quascirev.2011.11.021
## Saltre et al 2015, Quaternary Science Reviews 112: 128-37 / https://doi.org/10.1016/j.quascirev.2015.01.022
##
## INPUT DATA = text file with 2 columns (file extension and separators detected automatically)
## column 1 of 2 = 14C dates
## column 2 of 2 = standard deviations of 14C dates
## OUTPUT TO CONSOLE = list 'out' with two objects:
## object 1 of 2: out$calibration = non-calibrated and calibrated 14C dates with their errors
## object 2 of 2: out$griwm = extinction or arrival time (2.5%, 50% and 97.5% quartiles)
## OUTPUT TO WORKING DIRECTORY (OPTIONAL) = up to three text files
## file 1 of 3: non-calibrated and calibrated 14C dates with their errors
## file 2 of 3: re-sampled extinction or arrival times
## file 3 of 3: median extinction or arrival time with 95% Confidence Interval
######################################################################################
#' @name griwm
#'
#' @title Gaussian Re-sampled Inverse-Weighted McInerney (GRIWM)
#'
#' @description Estimate extinction time using a time series of dated observations with normal errors following Bradshaw et al (2012).
#'
#' @details The fossil record is highly fragmented such that the youngest fossil of an extinct species is unlikely to belong to the last individual of that species (so-called Signor-Lipps effect by Raup 1986). Among the methods available for estimating extinction time (see Wang and Marshall 2016), GRIWM controls for the Signor-Lipps effect based on the frequency of dated fossils prior to the extinction event.
#' Analogously, GRIWM can also estimate arrival time by conceptualizing that the oldest fossil of a species recorded in a region is unlikely to belong to the first individual colonizing that region.
#' GRIWM accounts for dating uncertainty by means of a re-sampling procedure assuming a normal distribution of every fossil age in calendar years (Bradshaw et al 2012; Saltre et al 2015). Re-sampling implemented using using \code{\link[stats]{rnorm}}.
#'
#' GRIWM defines extinction time as the calendar year at which the probability \code{alpha} of finding a fossil younger than the youngest known fossil is low (alpha = 0.05).
#' Analogously, GRIWM defines arrival time as the calendar year at which the probability \code{alpha} of finding a fossil older than the oldest known fossil is low (alpha = 0.05).
#' The method can be applied not only to fossil data but to any time series of historic observations each consisting of a date and an associated error.
#'
#' \code{griwm} can compute biased and unbiased estimates of extinction (or arrival) time.
#' The biased estimator follows Bradshaw et al (2012) and generates extinction (or arrival) time given a pre-defined probability \code{alpha}.
#' The biased estimator generates an absolute extinction (or arrival) time independent of a pre-defined probability \code{alpha}.
#'
#' \code{griwm} accepts any combination of chronological data (e.g., radiocarbon, ESR, OSL, TL, Uranium series).
#' Non-calibrated radiocarbon dates are calibrated in calendar years using \code{\link[rcarbon]{calibrate}} and expressed as a weighted average plus/minus a weighted standard deviation of the calibrated probability density function (Telford et al 2004).
#' Calibration curves include \code{"intcal20"}, \code{"intcal13nhpine16"} or \code{"intcal13"} for fossils from the Northern Hemisphere, \code{"shcal20"}, \code{"shcal13shkauri16"} or \code{"shcal13"} for fossils from the Southern Hemisphere, and \code{"marine20"} or \code{"marine13"} for marine fossils.
#' \code{"intcal20"}, \code{"shcal20"} and \code{"marine20"} allow calibrations up to 55000 calendar years Before Present, and the other calibration curves allow calibrations up to 50000 calendar years Before Present.
#' Dates from 1951 to present (post-bomb) not supported.
#' Default \code{cal_curve = "intcal20"}.
#'
#' Input data read from working directory using \code{\link[data.table]{fread}}; file extension and separators detected automatically.
#' Output data saved to working directory using \code{\link[data.table]{fwrite}}.
#' \code{cal_save = TRUE} saves the original data with dates in years (non-calibrated) and calendar years (calibrated). Non-calibrated dates only included if \code{calibra = TRUE}. Name of output file = outputGRIWM_calibration.
#' \code{resample_save = TRUE} saves re-sampled extinction (or arrival) values used to estimate median extinction (or arrival) time with 95 percent quartile ranges. Number of values set by \code{resample}. Name of output file = outputGRIWM_resampled_Extinction (\code{signor_lipps = "ext"}) or outputGRIWM_resampled_Arrival (\code{signor_lipps = "arr"}).
#' \code{griwm_save = TRUE} saves median extinction (or arrival) time with 95 percent quartile ranges. Name of output file = outputGRIWM_medianCI_Extinction (\code{signor_lipps = "ext"}) or outputGRIWM_medianCI_Arrival (\code{signor_lipps = "arr"}).
#'
#' @param chrono_data Text file read from working directory.
#' File must contain two columns without row names or column headings, column 1 includes dates, column 2 includes errors.
#' Unlimited number of dates (rows).
#' All dates must be Before Present, where Present = 1950.
#'
#' @param signor_lipps Character. GRIWM accounts for the Signor-Lipps effect at the youngest (extinction) or oldest (arrival) end of a time series of dated observations.
#' Default \code{signor_lipps = "ext"} for extinction time, and \code{signor_lipps = "arr"} for arrival time.
#'
#' @param biased Logical. Two estimators of extinction or arrival time.
#' Default is the original GRIWM estimator \code{biased = TRUE} dependent on parameter \code{alpha}.
#' Option \code{biased = FALSE} calculates an absolute time of extinction or arrival so is independent of parameter \code{alpha}. See details.
#'
#' @param alpha Numeric. Sets the probability of finding a new observation younger than the youngest known observation (extinction) or older than the oldest known observation (arrival) of population or species.
#' Default \code{alpha = 0.05}. Lacks functionality for GRIWM unbiased estimator \code{biased = FALSE}. See details.
#'
#' @param resample Numeric. Sets the number of samples taken from the normal distribution of each dated observation in the time series.
#' Default \code{resample = 10000}.
#'
#' @param radiocarbon Character or numeric. Specifies if the time series contains radiocarbon dates.
#' Default \code{radiocarbon = "all"} assumes that all records in the time series are non-calibrated-radiocarbon dated observations.
#' \code{radiocarbon = 0} assumes that all records are in calendar years (whether dated by radiocarbon or not).
#' For data with mixture of non-calibrated radiocarbon dates and/or calibrated radiocarbon dates and/or non-radiocarbon dates, \code{radiocarbon} must be a vector of integers indicating the position of non-calibrated radiocarbon dates by row number. For instance \code{radiocarbon = c(2)} for one non-calibrated radiocarbon date in row 2, or \code{radiocarbon = c(2:4,7,11)} for five non-calibrated radiocarbon dates in rows 2 to 4, 7 and 11.
#'
#' @param calibra Logical. Calibrates non-calibrated radiocarbon dates in calendar years.
#' Default \code{calibra = TRUE} assumes presence of non-calibrated radiocarbon dates if \code{radiocarbon = "all"} or \code{radiocarbon = c(row numbers)}.
#' See details.
#'
#' @param cal_curve Character. Sets calibration curve for calibrating non-calibrated radiocarbon dates if present. See details.
#'
#' @param upper14C Numeric. Sets the upper temporal boundary for radiocarbon calibrations.
#' Default \code{upper14C = 55000} for calibration curves \code{"intcal20"}, \code{"shcal20"} or \code{"marine20"}.
#' Use \code{upper14C = 50000} for calibration curves \code{"intcal13"}, \code{"shcal13"}, \code{"marine13"}, \code{"intcal13nhpine16"} or \code{"shcal13shkauri16"}
#'
#' @param cal_save Logical. Determines if calibrated records are saved as text file to working directory.
#' Default \code{cal_save = "FALSE"}.
#'
#' @param resample_save Logical. Determines if raw extinction (or arrival) times from each re-sampling iteration are saved as text file to working directory.
#' Default \code{resample_save = "FALSE"}.
#'
#' @param griwm_save Logical. Determines if median extinction (or arrival) time and 95 percent CI are saved as text file to working directory.
#' Default \code{griwm_save = "FALSE"}.
#'
#' @return Returns a list containing two elements:
#' \itemize{
#' \item {\code{griwm_arguments}} {Data.frame with 2 columns.
#'  Column 1 = arguments in function \code{griwm}, Column 2 = values selected by user.}\cr
#'  \item {\code{calibration_BP}} {Data.table with 5 columns and as many rows as dates in the input data.
#'  Column 1 = non-calibrated dates, Column 2 = non-calibrated errors, Column 3 = calibrated dates, Column 4 = calibrated errors, Column 5 = type of date (radiocarbon, non-radiocarbon).
#'  Rows ordered from younger to older dates.}\cr
#'  \item {\code{griwm}} {Data.frame with 3 columns and 3 rows.
#'  Columns 1, 2 and 3 for 2.5 percent quartile CI, 50 percent (median) and upper 97.5 percent quartile CI of extinction (or arrival) time, respectively.
#'  Time in calendar years BCE (row 1: BCE = Before Current Era where CE starts at year 0), BP (row 2: BP = Before Present where Present = year 1950), BCY (row 3: BCY = Before Current Year).}
#' }
#'
#' @examples
#' #None
#'
#' @references
#' Bradshaw, C.J.A., A. Cooper, C.S.M. Turney & B.W. Brook (2012) Robust estimates of extinction time in the geological record. Quaternary Science Reviews, 33, 14-19.\cr
#'
#' Raup, D.M. (1986) Biological extinction in earth history. Science, 231, 1528-1533.
#'
#' Saltre, F., B.W. Brook, M. Rodriguez-Rey, A. Cooper, C.N. Johnson, C.S.M. Turney & C.J.A Bradshaw (2015) Uncertainties in dating constrain model choice for inferring extinction time from fossil records. Quaternary Science Reviews, 112, 128-137.\cr
#'
#' Telford, R.J., Heegaard, E. & Birks, H.J.B. (2004) The intercept is a poor estimate of a calibrated radiocarbon age. Holocene, 14, 296-298.
#'
#' Wang, S.C. & C.R. Marshall (2016) Estimating times of extinction in the fossil record. Biology Letters, 12, 20150989.\cr
#'
#' @seealso \code{\link[rcarbon]{calibrate}} (package \bold{rcarbon}),
#' \code{\link[data.table]{data.table}} (package \bold{data.table})
#' \code{\link[data.table]{fread}} (package \bold{data.table})
#' \code{\link[data.table]{fwrite}} (package \bold{data.table})
#'
#' @importFrom stats na.omit
#' @importFrom stats quantile
#' @importFrom stats rnorm
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#'
#' @export
utils::globalVariables(c("i_date", "i_griwm", "i_iter", "i_resample")) #assign non-binding global variables (used as index in loops)
griwm <- function(chrono_data, #input data
                  signor_lipps = "ext", biased = TRUE, alpha = 0.05, resample = 10000, #griwm parameters
                  radiocarbon = "all", calibra = TRUE, cal_curve = "intcal20", upper14C = 55000, #calibration parameters
                  cal_save = FALSE, resample_save = FALSE, griwm_save = FALSE) #saving parameters
{

  ##############################################
  ## ABORT SESSION IF PARAMETERS ENTERED INCORRECTLY
  if (!signor_lipps %in% c("ext", "arr")) {
    stop("Parameter signor_lipps must be ext or arr, please check spelling\nComputation aborted\n")
  }
  if (alpha == 0) {
    stop("Probability alpha must be > 0 and <= 1, please check values\nComputation aborted\n")
  }
  if (alpha > 1) {
    stop("Probability alpha must be > 0 and <= 1, please check values\nComputation aborted\n")
  }
  if (is.character(radiocarbon[1]) & length(radiocarbon) > 1) {
    stop("Parameter radiocarbon must be all, 0 or a vector of integers, please check values\nComputation aborted\n")
  }
  if (is.character(radiocarbon[1]) & length(radiocarbon) == 1 & radiocarbon[1] != "all") {
    stop("Parameter radiocarbon must be all, 0, or a vector of integers, please check values\nComputation aborted\n")
  }
  if (length(radiocarbon) > 1 & !is.numeric(radiocarbon)) {
    stop("Parameter radiocarbon must be all, 0 or a vector of integers concanated with c() and indicating row numbers, please check values\nComputation aborted\n")
  }
  if (length(radiocarbon) > 1 & length(which(radiocarbon %in% 0)) > 0) {
    stop("Parameter radiocarbon must be all, 0, or a vector of integers concanated with c() and indicating row numbers (row 0 does not exist), please check values\nComputation aborted\n")
  }
  if (!cal_curve %in% c("intcal20", "shcal20", "marine20", "intcal13", "shcal13", "marine13", "intcal13nhpine16", "shcal13shkauri16")) {
    stop("Calibration curves accepted comprise intcal20, shcal20, marine20, intcal13, shcal13, marine13, intcal13nhpine16 or shcal13shkauri16, please modify or check spelling of parameter cal_curve\nComputation aborted\n")
  }
  if (cal_curve %in% c("intcal13", "shcal13", "marine13", "intcal13nhpine16", "shcal13shkauri16") & upper14C > 50000) {
    stop("Calibration window older than 50000 years only possible with calibration curves intcal20, shcal20 or marine20, please modify parameters cal_curve and/or upper14C\nComputation aborted\n")
  }

  ##############################################
  ## SAVE selected arguments
  griwm_args <- match.call()
  argsNames <- names(as.list(griwm_args))[-1] #get names of arguments
  argsValues <- as.character(griwm_args)[-1] #get values selected per argument
  griwm_args <- data.frame(argsNames, argsValues) #dataframe arguments and their values
  colnames(griwm_args) <- c("argument", "value") #rename columns
  rm(argsNames, argsValues) #remove objects

  ##############################################
  ## PREPARE DATA

  # LOAD data
  message("\n")
  if (signor_lipps == "ext") {message(paste("########GRIWM-extinction analysis started / Time = ", date(), "\n"))}
  if (signor_lipps == "arr") {message(paste("########GRIWM-arrival analysis started / Time = ", date(), "\n"))}
  message("Loading data\n")
  message("Input data assumed to be a text file, without row names or column headings, column 1 = dates in years Before Present (integers), column 2 = errors in years (integers)\n")
  chrono_data <- data.frame(data.table::fread(file = chrono_data))
  chrono_data <- chrono_data[, c(1:2)] #avoid reading column with unnoticed values beyond column 1 and 2

  ##################################
  ###BENCHMARKING MODULE
  # if (isTRUE(doBenchmark == 1)) {
  #   age <- sample(c(1:9707), size = 1)
  #   age <- sample(c(age:(age+4999)), size = datasetN, replace = FALSE)
  #   chrono_data <- data.frame(age)
  #   chrono_data$sd <- i_error
  #   rm(age)
  # }
  ##################################

  colnames(chrono_data) <- c("age", "sd") #rename columns
  dateN <- dim(chrono_data)[1] #number of dates in data
  dat_uncal <- chrono_data #create object for output file
  message(paste(dateN, " dates loaded\n", sep = ""))

  # MESSAGE type of data to console
  if(radiocarbon[1] == "all" & isTRUE(calibra)) {message("All dates are non-calibrated (calibra = TRUE) radiocarbon (radiocarbon = 'all')\n")}
  if(radiocarbon[1] == "all" & isFALSE(calibra)) {message("All dates are radiocarbon (radiocarbon = 'all') in calendar years (calibra = FALSE)\n")}
  if(radiocarbon[1] == 0) {message("All dates assumed to be non-radiocarbon (radiocarbon = 0) in calendar years\n")}
  if(is.numeric(radiocarbon[1]) & radiocarbon[1] != 0) {
    date14C <- length(radiocarbon)
    if (isTRUE(calibra)) {
      message(paste("Data contains ", date14C, " non-calibrated (calibra = TRUE) radiocarbon dates (radiocarbon = c(row numbers)) and ", dateN - date14C, " non-radiocarbon dates\n", sep=""))
    }
    if (isFALSE(calibra)) {
      message(paste("Data contains ", date14C, " calibrated (calibra = FALSE) radiocarbon dates (radiocarbon = c(row numbers)) and ", dateN - date14C, " non-radiocarbon dates\n", sep=""))
    }
  }

  # CALIBRATE radiocarbon data and convert to data.table
  #prevent error if users request calibration even if all input dates are in calendar years
  if (isTRUE(calibra) & radiocarbon[1] == 0) { #only if radiocarbon dates are NOT present
    chrono_data <- data.frame(as.numeric(chrono_data[, 1]), as.numeric(chrono_data[, 2])) #dataframe repository and coerce mean and sd to numeric
    colnames(chrono_data) <- c("age", "sd") #rename columns
    dat_cal <- data.table::data.table(chrono_data)
  }
  #subset non-calibrated radiocarbon dates
  if (isTRUE(calibra) & radiocarbon[1] != 0) { #radiocarbon dates present
    #subset radiocarbon dates if list includes dates obtained through chronological methods
    if (is.character(radiocarbon[1]) == TRUE) { #all dates are non-calibrated radiocarbon
      chrono_data <- chrono_data
    }
    if (is.numeric(radiocarbon[1]) == TRUE) { #mixed radiocarbon and non-radiocarbon dates prsent
      dat1 <- chrono_data[-radiocarbon, ] #object contains calibrated dates
      chrono_data <- chrono_data[radiocarbon, ] #object contains non-calibrated dates
    }

    #calibrate
    chrono_data$curve <- cal_curve
    calRaw <- rcarbon::calibrate(chrono_data[, 1], errors = chrono_data[, 2], calCurves = chrono_data[, 3], timeRange = c(upper14C, 0), F14C = TRUE, eps = 0)
    chrono_data <- matrix(NA, nrow = nrow(chrono_data), ncol = 2) #repository of calibrated point estimates of radiocarbon dates
    for (i_date in 1:length(calRaw)) {
      cal_date <- as.numeric(data.frame(calRaw$grids[i_date])[, 1]) #get calibrated values
      cal_probs <- as.numeric(data.frame(calRaw$grids[i_date])[, 2]) #get PDF
      cal_w_mean <- round(sum(cal_date * cal_probs), 0) #estimate calibrated weighted mean
      cal_w_sd <- round(sqrt(sum(((cal_date - cal_w_mean) ^ 2) * cal_probs)), 0) #estimate calibrated weighted sd
      chrono_data[i_date, ] <- c(cal_w_mean, cal_w_sd) #store cal weighted mean and sd
      # message(paste("date", i_date))
    }
    chrono_data <- data.frame(as.numeric(chrono_data[, 1]), as.numeric(chrono_data[, 2])) #dataframe repository and coerce mean and sd to numeric
    colnames(chrono_data) <- c("age", "sd") #rename columns

    #rbind calibrated radiocarbon dates with non-radiocarbon dates, re-order by row number and assign data to data.table object
    if (is.numeric(radiocarbon[1]) == TRUE & radiocarbon[1] != 0) {
      rownames(chrono_data) <- radiocarbon #row number for radiocarbon dates
      chrono_data <- data.frame(rbind(chrono_data, dat1)) #rbind radiocarbon and non-radiocarbon dates
      chrono_data <- data.table::data.table(chrono_data[order(as.numeric(rownames(chrono_data))), ]) #get original row order
      rm(dat1) #remove object
    }
    #create object for output file when at least one date is radiocarbon
    dat_cal <- data.table::data.table(chrono_data)
    rm(calRaw, cal_date, cal_probs, cal_w_mean, cal_w_sd) #clean space from R environment
  }

  #assign data to data.table object when all input dates are calibrated
  if (isFALSE(calibra)) {
    chrono_data <- data.table::data.table(chrono_data)
    colnames(chrono_data) <- c("age", "sd") #rename columns
    dat_cal <- chrono_data #create object for output file
  }

  #assign data to data.table object when no input dates are calibrated or radiocarbon
  if (is.character(radiocarbon[1]) == TRUE | radiocarbon[1] == 0) {
    chrono_data <- data.table::data.table(chrono_data)
  }

  # SORT data
  if (signor_lipps == "ext") {chrono_data <- chrono_data[order(chrono_data[, 1], decreasing = FALSE), c(1:2)]} #extinction
  if (signor_lipps == "arr") {chrono_data <- chrono_data[order(chrono_data[, 1], decreasing = TRUE), c(1:2)]} #arrival

  # MESSAGE calibration options to console
  if (isTRUE(calibra) & is.character(radiocarbon[1])) {
    message(paste("Non-calibrated radiocarbon dates calibrated in calendar years Before Present with ", cal_curve, sep = ""), "\n")
  }
  if (isTRUE(calibra) & is.numeric(radiocarbon[1]) & radiocarbon[1] != 0) {
    message(paste("Non-calibrated radiocarbon dates calibrated in calendar years Before Present with ", cal_curve, sep = ""), "\n")
  }
  if (isTRUE(calibra) & radiocarbon[1] == 0) {
    message("Radiocarbon calibration requested (calibra = TRUE) but all records assumed to be in calendar years Before Present (radiocarbon = 0), so calibration not implemented \n")
  }
  if (isFALSE(calibra) & is.numeric(radiocarbon[1]) & radiocarbon[1] != 0) {
    message("Radiocarbon calibration not requested (calibra = FALSE). Since input data contains radiocarbon and non-radiocarbon dates, all dates assumed to be in calendar years Before Present\n")
  }

  ##############################################
  ## MESSAGE GRIWM PARAMETERS TO CONSOLE
  message("########GRIWM re-sampling started\n")

  if (signor_lipps == "ext") {
    message(paste("Extinction time (signor_lipps == 'ext') from ", resample, " re-sampled estimates (resample = ", resample, ")", sep = ""), "\n")
  }
  if (signor_lipps == "arr") {
    message(paste("Arrival time (signor_lipps == 'arr') from ", resample, " re-sampled estimates (resample = ", resample, ")", sep = ""), "\n")
  }
  if (isTRUE(biased)) {
    message(paste("Biased estimator applied (biased = TRUE) for alpha = ", alpha, " (alpha = ", alpha, ")" , sep = ""), "\n")
  } else {
    message("Unbiased estimator applied (biased = FALSE)\n")
  }

  ## COMPUTE GRIWM

  # VECTORIZE data, set griwm parameters, define output repositories
  # datF <- c(chrono_data[, 1])$age #get dates
  # sd.vec <- c(chrono_data[, 2])$sd #get errors
  datF <- data.frame(chrono_data)$age #get dates
  sd.vec <- data.frame(chrono_data)$sd #get errors
  rm(chrono_data) #remove data object
  k <- length(datF) #get number of dates
  out <- matrix(0, 1, 3) #create repository of final griwm results
  T.up.vec <- T.mci.vec <- w.T.mci.vec <- rep(0, resample) #create repository of temporal griwm results

  # LOOP re-sampling
  #set progress bar
  pb <- txtProgressBar(min = 1, max = resample, style = 3)

  #run loop
  for (i_iter in 1:resample) { #loop over i_iter iterations (default = 10000)

    #run progress bar
    setTxtProgressBar(pb = pb, i_iter)

    #fix random generator
    set.seed(i_iter)

    #re-sample one value from each of the k dates and store in vector date.samp = re-sampled time series
    date.samp <- rep(0,k) #repository of re-sampled dates
    for (i_resample in 1:k) { #loop over k dates
      date.samp[i_resample] <- round(rnorm(1, datF[i_resample], sd.vec[i_resample]))
    }

    #set estimate: extinction versus arrival
    if (signor_lipps == "ext") {date.samp <- sort(date.samp)}
    if (signor_lipps == "arr") {date.samp <- sort(date.samp, decreasing = TRUE)}

    #estimate weights per date
    last.diff <- 1 / (date.samp - date.samp[1])[-1]
    weight <- last.diff / last.diff[1]
    if (last.diff[1] == Inf) {
      weight <- last.diff / last.diff[2]
      weight <- weight[-1]
    }

    #estimate extinction or arrival parameters
    ldate <- length(date.samp)
    T.mci.lst.vec <- rep(0, ldate - 1)
    for (i_griwm in 1:(ldate-1)) { #loop over ldate = k - 1 dates
      date.it <- date.samp[1:(1 + i_griwm)]
      date.age.it <- date.samp[1:(1 + i_griwm)]
      date.mci.it <- rev(max(date.it) + 1 - date.it)
      k <- length(date.it)
      t.n <- date.mci.it[k]
      n <- k
      T.rng <- t.n - date.mci.it[1]
      #biased versus unbiased estimate
      if (isTRUE(biased)) {
        i = t.n - t.n * log(alpha) / n #GRIWM biased estimate (original)
      }
      if (isFALSE(biased)) {
        i = t.n * (n + 1) / n #GRIWM unbiased estimate
      }
      T.mci.lst.vec[i_griwm] <- max(date.it) + 1 - i
    }

    #inverse-weight dates
    if (last.diff[1] == Inf) {
      w.T.mci.vec[i_iter] <- round((sum(weight * T.mci.lst.vec[-1])) / sum(weight), 0)
    }
    if (last.diff[1] != Inf) {
      w.T.mci.vec[i_iter] <- round((sum(weight * T.mci.lst.vec)) / sum(weight), 0)
    }
  }
  #close progress bar
  close(pb)

  # CALCULATE median extinction (or arrival) time and 95% CI
  T.wmci.vec.CI <- quantile(na.omit(w.T.mci.vec), probs = c(0.025, 0.5, 0.975));
  out[1] <- round(T.wmci.vec.CI[1], 0) #2.5% quartile (lower CI bound)
  out[2] <- round(T.wmci.vec.CI[2], 0) #50% quartile (median)
  out[3] <- round(T.wmci.vec.CI[3], 0) #97.5% quartile (upper CI bound)
  out <- data.frame(out) #convert to dataframe

  #convert BP = Before Present
  outBCE <- out - 1949 #to BCE = Before Current Era
  outBCY <- out  + as.numeric(format(Sys.Date(), "%Y")) - 1950 #to BCY = Before Current Year
  out <- data.frame(rbind(outBCE, out, outBCY)) #data.frame the three time scales

  if (signor_lipps == "ext") {rownames(out) <- c("Extinction_Time_BCE", "Extinction_Time_BP", "Extinction_Time_BCY")}
  if (signor_lipps == "arr") {rownames(out) <- c("Arrival_Time_BCE", "Arrival_Time_BP", "Arrival_Time_BCY")}
  colnames(out) <- c("2.5CI", "median", "97.5CI")

  # REMOVE GRIWM loop inputs from console
  rm(datF, sd.vec, k, date.samp, last.diff, weight, ldate, date.it, date.age.it, date.mci.it, t.n, n, T.rng, outBCE, outBCY)

  # MESSAGE end of computations
  message("\n")
  message("########GRIWM re-sampling ended\n")
  if (signor_lipps == "ext") {message("GRIWM extinction time reported in calendar years BP, BCE and BCY\n")}
  if (signor_lipps == "arr") {message("GRIWM arrival time reported in calendar years BP, BCE and BCY\n")}
  message("BP = Before Present (Present = year 1950)\n")
  message("BCE = Before Current Era (Current Era starts at year 0)\n")
  message(paste("BCY = Before Current Year (Current Year = ", as.numeric(format(Sys.Date(), "%Y")), ")", sep = ""), "\n")
  if (signor_lipps == "ext") {message("Note = extinction time in BP or BCE will be negative if earlier than 1950 or 0, respectively\n")}
  if (signor_lipps == "arr") {message("Note = arrival time in BP or BCE will be negative if earlier than 1950 or 0, respectively\n")}

  ##############################################
  ## SAVE GRIWM MEDIAN EXTINCTION (OR ARRIVAL) TIME AND 95% CI
  if (griwm_save == TRUE & signor_lipps == "ext") {
    data.table::fwrite(x = out, file = "outputGRIWM_medianCI_Extinction", col.names = TRUE)
    message("Median extinction time and 95% confidence interval saved to working directory (griwm_save = TRUE) with file name = outputGRIWM_medianCI_Extinction\n")
  }
  if (griwm_save == TRUE & signor_lipps == "arr") {
    data.table::fwrite(x = out, file = "outputGRIWM_medianCI_Arrival", col.names = TRUE)
    message("Median arrival time and 95% confidence interval saved to working directory (griwm_save = TRUE) with file name = outputGRIWM_medianCI_Arrival\n")
  }

  ##############################################
  ## SAVE RE-SAMPLED EXTINCTION (OR ARRIVAL) TIMES
  if (signor_lipps == "ext" & isTRUE(resample_save)) {
    data.table::fwrite(x = data.frame(w.T.mci.vec), file = "outputGRIWM_resampled_Extinction", col.names = FALSE)
    message("Re-sampled extinction times saved to working directory (resample_save = TRUE) with file name = outputGRIWM_resampled_Extinction\n")
  }
  if (signor_lipps == "arr" & isTRUE(resample_save)) {
    data.table::fwrite(x = data.frame(w.T.mci.vec), file = "outputGRIWM_resampled_Arrival", col.names = FALSE)
    message("Re-sampled arrival times saved to working directory (resample_save = TRUE) with file name = outputGRIWM_resampled_Arrival\n")
  }
  # REMOVE GRIWM loop outputs from console
  rm(T.up.vec, T.mci.vec, w.T.mci.vec, T.wmci.vec.CI, T.mci.lst.vec) #remove GRIWM outputs from console

  ##############################################
  ## SAVE CALIBRATION

  # CREATE vector assigning dates to radiocarbon or non-radiocarbon dating
  if (is.character(radiocarbon) == TRUE) {chrono_type <- rep("radiocarbon", dateN)}
  if (is.numeric(radiocarbon) == TRUE & radiocarbon[1] != 0) {
    chrono_type <- rep("non-radiocarbon", dateN)
    chrono_type[radiocarbon] <- "radiocarbon"
  }
  if (radiocarbon[1] == 0) {chrono_type <- rep("non-radiocarbon", dateN)}

  # ASSIGN calibrated results to object calibration, and (if requested) save object to working directory
  if (isTRUE(calibra)) { #calibration requested
    #coerce input data to NA calibrated dates present
    if (radiocarbon[1] == 0) {dat_uncal[c(1:dateN),] <- NA} #all calibrated dates coerced to NA in dat_uncal
    if (is.numeric(radiocarbon[1]) & radiocarbon[1] != 0) {non_radiocarbon <- setdiff(c(1:dateN), radiocarbon); dat_uncal[non_radiocarbon, ] <- NA} #only calibrated dates coerced to NA  in dat_uncal
    if (is.character(radiocarbon[1])) {dat_uncal <- dat_uncal} #no dates coerced to NA  in dat_uncal
    #save calibration object
    calibration <- data.table::data.table(uncal.age = dat_uncal[, 1], uncal.sd = dat_uncal[, 2], cal.age = dat_cal[, 1], cal.sd = dat_cal[, 2], type = chrono_type)
    calibration <- calibration[order(calibration$cal.age), ] #order from younger to older dates
    if (isTRUE(cal_save) & radiocarbon[1] != 0) { #save non-calibrated and calibrated dates
      data.table::fwrite(x = data.frame(calibration), file = "outputGRIWM_calibration", col.names = TRUE)
      message("Dates in years and calendar years saved to working directory (cal_save = TRUE) with file name = outputGRIWM_calibration\n")
    }
    if (isTRUE(cal_save) & radiocarbon[1] == 0) { #save calibrated dates
      data.table::fwrite(x = data.frame(calibration)[, c(3:5)], file = "outputGRIWM_calibration", col.names = TRUE)
      message("All input dates were in calendar years (radiocarbon = 0) and have been saved to working directory (cal_save = TRUE) with file name = outputGRIWM_calibration\n")
    }
  }
  if (isFALSE(calibra)) { #calibration not requested
    calibration <- data.table::data.table(cal.age = dat_cal[, 1], cal.sd = dat_cal[, 2], type = chrono_type)
    calibration <- calibration[order(calibration$cal.age), ] #order from younger to older dates
    if (isTRUE(cal_save)) { #save calibrated dates
      data.table::fwrite(x = data.frame(calibration), file = "outputGRIWM_calibration", col.names = TRUE)
      message("All input dates are in calendar years (calibra = FALSE) and have been saved to working directory (cal_save = TRUE) with file name = outputGRIWM_calibration\n")
    }
  }

  ##############################################
  ## PRINT RESULTS TO CONSOLE AND CLEAN R ENVIRONMENT
  if (signor_lipps == "ext") {message(paste("########GRIWM-extinction analysis ended / Time = ", date(), "\n"))}
  if (signor_lipps == "arr") {message(paste("########GRIWM-arrival analysis ended / Time = ", date(), "\n"))}
  out <- list(griwm_arguments = griwm_args, calibration_BP = calibration, griwm = out)
  return(out)
  rm(list = ls(all.names = TRUE)) #clean R environment
}

## END
##############################################
