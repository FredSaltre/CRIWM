#######################################################################################
## package: Rextinct v.1
## function: Rextinct::calendar by Salvador Herrando-Perez (salherra@gmail.com)
## Calibration of radiocarbon dates
## Date = 01 April 2021
##
##
## INPUT DATA = text file with 3 columns or 4 columns (file extension and separators detected automatically)
## column 1 of 4 = 14C dates
## column 2 of 4 = standard deviations of 14C dates
## column 3 of 4 = date alpha-numeric code
## column 4 of 4 = names of calibration curves
## OUTPUT TO CONSOLE = list with calibration results
## OUTPUT TO WORKING DIRECTORY (OPTIONAL) = text file with calibration results
######################################################################################
######################################################################################
#' @name calendar
#'
#' @title Calibration of radiocarbon dates
#'
#' @description Calibrate radiocarbon dates into point-estimate intercepts (mean, median, mode) in calendar years and describe the shape of calibrated probability density functions (PDF).
#'
#' @details The concentration of the radioactive isotope carbon 14 (14C) in an organism starts decaying at a near-constant rate (Libby half-life = 5568 years) after death. The date of death can be then estimated by measuring the amount of 14C in the tissue of the organism in question as radioactivity (beta counting) or number of 14C atoms (accelerator mass spectrometry). Organic materials older some 48 thousand years retain negligible 14C content today so, currently, radiocarbon cannot be used for dating samples beyond that time interval.
#'
#' Critically, 14C decay rates are non-linear so raw 14C dates in years (as obtained from dating facilities) must be adjusted to 14C decay trends based on calibration curves (as obtained from natural archives including corals, sediments, speleothems and tree rings). Calibration curves differ depending on the region (Northern versus Southern Hemisphere) or environment (terrestrial versus marine) of provenance of the sample. Overall, 14C calibration aims to express a raw 14C date plus/minus error (in years) in a probability density function (PDF) that quantifies the probability that the sample has a range of possible ages (in calendar years) given the amount of 14C present in layers of natural archives of known age (Keenan 2012).
#'
#' \code{calendar} emulates the calibration protocol of the field-standard software \code{oxcal} (Bronk Ramsey 1995, 2007) using \code{\link[rcarbon]{calibrate}} (Crema & Bevan 2021) in F14C (radiocarbon concentration) space (Heaton et al. 2020a), and PDFs are normalized to a maximum probability of 1. The overarching goal is to provide a friendly tool for calibrating 14C dates, saving simple and tidy calibration outputs to the working directory, and computing statistical descriptors for examining PDF shape. The function employs one single computer core.
#'
#' Calibration curves include \code{"intcal20"}, \code{"intcal13nhpine16"} or \code{"intcal13"} for fossils from the Northern Hemisphere (Reimer et al 2020), \code{"shcal20"}, \code{"shcal13shkauri16"} or \code{"shcal13"} for fossils from the Southern Hemisphere (Hogg et al 2020), and \code{"marine20"} or \code{"marine13"} for marine fossils (Heaton et al 2020b). \code{"intcal20"}, \code{"shcal20"} and \code{"marine20"} allow calibrations up to 55000 calendar years Before Present , and the other calibration curves allow calibrations up to 50000 calendar years Before Present.
#'
#' Input data read from working directory using \code{\link[data.table]{fread}}; file extension and separators detected automatically. Data file must contain three or four columns. Columns 1 and 2 must always contain the raw radiocarbon dates and errors, respectively, and column 3 must always contain the dating-lab identifier (for instance, OxA-13579 or UCIAMS-123456). Use a 3-column data file in combination with \code{cal_curve} (for instance, \code{cal_curve = "intcal20"}), if all 14C dates should be calibrated with the same calibration curve. Use column four to include the name of calibration curves when different 14C dates should be calibrated with different calibration curves; argument \code{cal_curve} will be non-functional. Calibration results saved as a text file to working directory using \code{\link[data.table]{fwrite}}.
#'
#' \code{calendar} computes three types of statistical descriptors, namely intercepts, moments and multi-modality:
#'
#' INTERCEPTS: calibrated 14C dates expressed as the calibrated date at which a non-calibrated date intercepts a calibration curve (Telford, Heegaard and Birks 2003).
#' We have implemented the weighted average (mean of calibrated dates with non-zero probability weighted by their PDF probabilities), the median (calibrated date at which 50 percent of the PDF area is accounted for) and the mode (calibrated date at which the PDF probability is highest).
#' Confidence intervals of the weighted average (column \code{calW_sd} in output object \code{cal_intercepts}) can be plus/minus 1 standard deviation (\code{cal_ci = "sd1"} comprising 68.3 percent of the PDF, or 2 standard deviations (\code{cal_ci = "sd1"} comprising 95.4 percent of the pdf. For consistency with the weighted average, confidence intervals of the median (columns \code{calQ_lowerCI} and \code{calQ_upperCI} in output object \code{cal_intercepts}) can have quartile ranges from 0.1585 to 0.8415 (\code{cal_ci = "sd1"} comprising a cumulative 68.3 percent of the PDF area, or from 0.0230 to 0.9770 (\code{cal_ci = "sd2"} comprising a cumulative 95.4 percent of the PDF area.
#'
#' MOMENTS: the shape of the PDF of a 14C date is summarized as the mean, variance, skewness (\code{\link[moments]{skewness}}) and kurtosis (\code{\link[moments]{kurtosis}}).
#'
#' MULTI-MODALITY: the statistically supported modes (highest posterior density regions falling in the shortest intervals of the PDF) for 1 (\code{cal_ci = "sd1"} or 2 (\code{cal_ci = "sd2"} standard deviations using a (\code{summary} of a \code{\link[rcarbon]{calibrate}} object.
#' In output object \code{cal_modality}, column \code{modality_SD1} shows number of statistically supported modes at 1SD, and column \code{modality_SD2} shows number of statistically supported modes at 1SD and 2SD.
#'
#' @param chrono_data Name of text file read from working directory.
#' File must contain three or four columns without row names or column headings. Column 1 = non-calibrated dates, column 2 = non-calibrated errors, column 3 = date identifiers, column 4 (optional) = name of calibration curve per date.
#' Unlimited number of dates (rows).
#' All dates must be Before Present, where Present = 1950.
#' See details.
#'
#' @param cal_curve Character. Sets unique calibration curve for calibrating all radiocarbon dates.
#' Default \code{cal_curve = "intcal20"}.
#' See details.
#'
#' @param upper14C Numeric. Sets the upper temporal boundary for radiocarbon calibrations.
#' Default \code{upper14C = 55000} for calibration curves \code{"intcal20"}, \code{"shcal20"} or \code{"marine20"}.
#' Use \code{upper14C = 50000} for calibration curves \code{"intcal13"}, \code{"shcal13"}, \code{"marine13"}, \code{"intcal13nhpine16"} or \code{"shcal13shkauri16"}.
#' Dates from 1951 to present (post-bomb) not supported.
#' See details.
#'
#' @param cal_ci Character. Determines width of confidence interval of calibrated mean and median.
#' Default \code{cal_ci = "sd1"} for 68.3 percent of the distribution of the calibrated mean (plus/minus 1 standard deviation) or 68.3 percent quartile range of the median.
#' \code{cal_ci = "sd2"} for 95.4 percent (plus/minus 1 standard deviation).
#' See details.
#'
#' @param cal_save Character. Determines if calibrated outputs are saved as text file to working directory.
#' Default \code{cal_save = "intercepts"} saves calibrated intercepts (weighted mean, median, mode).
#' \code{cal_save = "full"} saves PDF moments (mean, variance, skewness, kurtosis), statistically supported modes, and calibrated intercepts.
#' \code{cal_save = "none"} cancels saving of outputs.
#' See details.
#'
#' @return Returns a list containing four elements:
#' \itemize{
#' \item {\code{calibration_arguments}} {Data.frame with 2 columns.
#'  Column 1 = arguments in function \code{calendar}, Column 2 = values selected by user.}\cr
#'  \item {\code{cal_moments}} {Data.table with 6 columns and as many rows as dates in the input data.
#'  Column 1 = date identifiers, Column 2 = PDF mean, Column 3 = PDF variance, Column 4 = PDF skewness, Column 5 = PDF kurtosis, Column = calibration curves.}\cr
#'  \item {\code{cal_modality}} {Data.table with at least 5 columns and as many rows as dates in the input data.
#'  Column 1 = date identifiers, Column 2 = number of PDF modes at +/- 1SD, Column 3 = number of modes at +/- 2SD, last column = calibration curve.
#'  Remaining intermediate columns show upper and lower boundaries of each statistically supported mode in the PDF of each date.}\cr
#'  \item {\code{cal_incercepts}} {Data.table with 10 columns and as many rows as dates in the input data.
#'  Column 1 = date identifiers, Column 2 = non-calibrated date, Column 3 = non-calibrated error,
#'  Column 4 = calibrated weighted mean, Column 5 = calibrated weighted error, Column 6 = calibration curve,
#'  Column 7 = calibrated median (50 percent quartile), Column 8 = calibrated lower quartile, Column 9 = calibrated upper quartile, Column 10 = calibrated mode.}
#' }
#'
#' @examples
#' #None
#'
#' @references
#' Bronk Ramsey, C. (1995) Radiocarbon calibration and analysis of stratigraphy: the OxCal program. Radiocarbon, 37, 425-430.\cr
#'
#' Bronk Ramsey, C. (2017) Methods for summarizing radiocarbon datasets. Radiocarbon, 59, 1809-1833.\cr
#'
#' Crema, E.R. & Bevan, A. (2021) Inference from large sets of radiocarbon dates: software and methods. Radiocarbon, 63, 23-39.\cr
#'
#' Heaton, T.J. et al (2020a) The IntCal20 approach to radiocarbon calibration curve construction: a new methodology using Bayesian splines and errors-in-variables. Radiocarbon, 62, 821-863.
#'
#' Heaton, T.J. et al (2020b) Marine20 - The marine radiocarbon age calibration curve (0-55000 cal BP). Radiocarbon, 62, 779-820.\cr
#'
#' Hogg, A.G. et al (2020) SHCal20 Southern Hemisphere calibration, 0-55000 years cal BP. Radiocarbon, 62, 759-778.\cr
#'
#' Keenan, D.J. (2012) Calibration of a radiocarbon age. Nonlinear Processes in Geophysics, 19, 345-350.\cr
#'
#' Reimer, P.J. et al (2020) The IntCal20 Northern Hemisphere radiocarbon age calibration curve (0-55 cal kBP). Radiocarbon, 62, 725-757.\cr
#'
#' Telford, R.J., Heegaard, E. & Birks, H.J.B. (2004) The intercept is a poor estimate of a calibrated radiocarbon age. Holocene, 14, 296-298.\cr
#'
#' @seealso \code{\link[rcarbon]{calibrate}} (package \bold{rcarbon})
#' \code{\link[data.table]{data.table}} (package \bold{data.table})
#' \code{\link[data.table]{fread}} (package \bold{data.table})
#' \code{\link[data.table]{fwrite}} (package \bold{data.table})
#' \code{\link[moments]{skewness}} (package \bold{moments})
#' \code{\link[moments]{kurtosis}} (package \bold{moments})
#'
#' @importFrom stats approx
#' @importFrom stats var
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#'
#' @export
utils::globalVariables(c("i_date")) #assign non-binding global variables (used as index in loops)
calendar <- function(chrono_data,
                     cal_curve = "intcal20", upper14C = 55000,
                     cal_ci = "sd1", cal_save = "intercepts")
{

  ##############################################
  ## ABORT SESSION IF PARAMETERS ENTERED INCORRECTLY
  if (!cal_ci %in% c("sd1", "sd2")) {
    stop("Parameter cal_ci must be 'sd1' or 'sd2', please check spelling\nComputation aborted\n")
  }
  if (!cal_save %in% c("intercepts", "full", "none")) {
    stop("Parameter cal_save must be 'intercepts', 'full' or 'none', please check spelling\nComputation aborted\n")
  }
  if (!cal_curve %in% c("intcal20", "shcal20", "marine20", "intcal13", "shcal13", "marine13", "intcal13nhpine16", "shcal13shkauri16")) {
    stop("Please check spelling of names of calibration curves\nComputation aborted\n")
  }
  if (upper14C > 50000 & cal_curve %in% c("intcal13", "shcal13", "marine13", "intcal13nhpine16", "shcal13shkauri16")) {
    stop("Calibration series 13 and 16 can only be used up to 50000 years Before Present, please make upper14C no higher than 50000\nComputation aborted\n")
  }

  ##############################################
  ## SAVE selected arguments
  calendar_args <- match.call()
  argsNames <- names(as.list(calendar_args))[-1] #get names of arguments
  argsValues <- as.character(calendar_args)[-1] #get values selected per argument
  calendar_args <- data.frame(argsNames, argsValues) #dataframe arguments and their values
  colnames(calendar_args) <- c("argument", "value") #rename columns
  rm(argsNames, argsValues) #remove objects


  ##############################################
  ## PREPARE DATA

  # LOAD data
  message("\n")
  message("Loading data\n")
  message("Input data assumed to be a text file, without row names or column headings,
          column 1 = dates in years Before Present (integer)
          column 2 = errors in years (integer)
          column 3 = date identifier (character)
          column 4 (optional) = names of calibration curves (character)\n")
  chrono_data <- data.frame(data.table::fread(file = chrono_data))
  dateN <- dim(chrono_data)[1] #number of dates in data

  #rename columns
  if (dim(chrono_data)[2] == 4) { #four columns
    message("Data contains 4 columns (ages, errors, identifiers, calibration curves)\n")
    colnames(chrono_data) <- c("age", "sd", "id", "curve")
  }
  if (dim(chrono_data)[2] == 3) { #three columns
    message("Data contains 3 columns (ages, errors, identifiers)\n")
    chrono_data$curve <- cal_curve #add column with names of calibration curve
    colnames(chrono_data) <- c("age", "sd", "id", "curve")
  }

  ##################################
  ###BENCHMARKING MODULE
  # if (dim(fossil_data)[2] == 2 & isTRUE(doBenchmark == 0)) {
  #   fossil_data$id <- paste("id", c(1:dateN), sep="")
  #   fossil_data$curve <- i_curveM
  #   colnames(fossil_data) <- c("age", "sd", "id", "curve")
  # }
  # if (isTRUE(doBenchmark == 1)) {
  #   age <- sample(c(1:14706), size = 1)
  #   age <- sample(c(age:(age+14706)), size = i_datN, replace = FALSE)
  #   fossil_data <- data.frame(age)
  #   fossil_data$sd <- i_error
  #   fossil_data$id <- paste("date", c(1:dim(fossil_data)[1]), sep="")
  #   fossil_data$curve <- cal_curve
  #   dateN <- dim(fossil_data)[1] #number of dates in data
  #   rm(age)
  # }
  ##################################

  # ABORT calibration if identifiers (column 3) are not unique for each date
  codN <- length(which(table(chrono_data$id) == 1)) #frequency of duplicate ids
  if (codN - dateN != 0) {
    stop("Calibration aborted\nSome identifiers (column 3) have been assigned to more than one date\nTo proceed with calibration, please provide unique identifiers to each date and reload the input data\n")
  }
  if (dim(chrono_data)[1] > 1) {message(paste(dateN, " radiocarbon dates loaded\n", sep = ""))}
  if (dim(chrono_data)[1] == 1) {message(paste(dateN, " radiocarbon date loaded\n", sep = ""))}
  rm(codN, dateN)

  ##############################################
  ## CALIBRATE dates
  message("\n")
  message("Calibrating dates entered in years Before Present\n")

  # CALIBRATE
  #message number of calibration curves
  if (length(table(chrono_data$curve)) > 1) {
    message("Different dates calibrated with different calibration curves\n")
  }
  if (length(table(chrono_data$curve)) == 1) {
    message("All dates calibrated with the same calibration curve = ", row.names(table(chrono_data$curve)), "\n")
  }
  #timeRange = sets temporal window
  #F14C = calibration in F14C domain (radiocarbon concentration) rather than temporal domain (radiocarbon age)
  #eps = range of probability density function considered (eps = 0 captures full PDF)
  calRaw <- rcarbon::calibrate(x = chrono_data[, "age"], errors = chrono_data[, "sd"], ids = chrono_data[, "id"], calCurves = chrono_data[, "curve"], timeRange = c(upper14C,0), F14C = TRUE, eps = 0)
  calib <- summary(calRaw) #summary of calibration output

  message("Calibration in calendar years Before Present completed\n")

  ##############################################
  ## ESTIMATE calibrated intercepts
  message("\n")
  message("Computing destriptive stats of calibration\n")

  #set progress bar
  pb <- txtProgressBar(min = 0, max = length(calRaw), style = 3)

  # Weighted average and standard deviation
  repo_cal <- matrix(NA, nrow = length(calRaw), ncol = 10) #ncol = number of metrics computed

  # # create pdf object
  # if (isTRUE(plot_save)) {pdf("calPlots.pdf", onefile = FALSE)}

  #loop one date at a time
  for (i_date in 1:length(calRaw)) {

    #run progress bar
    setTxtProgressBar(pb = pb, i_date)

    # get dates and probs
    cal_date <- as.numeric(data.frame(calRaw$grids[i_date])[, 1]) #get calibrated values
    cal_probs <- as.numeric(data.frame(calRaw$grids[i_date])[, 2]) #get PDF
    #calculate weighted average and sd
    cal_w_mean <- round(sum(cal_date * cal_probs), 0) #estimate calibrated weighted mean
    if (cal_ci == "sd1") { #calibrated weighted SD x 1 (68.3% CI)
      cal_w_sd <- round(sqrt(sum(((cal_date - cal_w_mean) ^ 2) * cal_probs)), 0)
    }
    if (cal_ci == "sd2") { #calibrated weighted SD x 2 (95.4% CI)
      cal_w_sd <- round(2 * sqrt(sum(((cal_date - cal_w_mean) ^ 2) * cal_probs)), 0)
    }

    # Median
    g <- data.frame(cal_date, cal_probs) #dataframe with calibrated date (1st column) and probability (second column)
    calQ_med <- round(approx(cumsum(g[, "cal_probs"]), g[, "cal_date"], ties = mean, xout = 0.5)$y, 0) #50% quartile = median / interpolates (approx) calibrated date at quartile 0.5 in PDF scaled from 0 to 1

    #quartiles
    if (cal_ci == "sd1") { #estimate calibrated weighted SD (68.3% of the PDF)
      calQ_upperCI <- round(approx(cumsum(g[, "cal_probs"]), g[, "cal_date"], ties = mean, xout = 0.1585)$y, 0) #upper quartile = lower quartile with interpolation (approx)
      calQ_lowerCI <- round(approx(cumsum(g[, "cal_probs"]), g[, "cal_date"], ties = mean, xout = 0.8415)$y, 0) #lower quartile = lower quartile with interpolation (approx)
    }
    if (cal_ci == "sd2") { #estimate calibrated weighted SD
      calQ_upperCI <- round(approx(cumsum(g[, "cal_probs"]), g[, "cal_date"], ties = mean,xout = 0.0230)$y, 0) #upper quartile = lower quartile with interpolation (approx)
      calQ_lowerCI <- round(approx(cumsum(g[, "cal_probs"]), g[, "cal_date"], ties = mean, xout = 0.9770)$y, 0) #lower quartile = upper quartile with interpolation (approx)
    }

    # Mode
    mode_I <- which(g[, "cal_probs"] == max(g[, "cal_probs"]))
    if (length(mode_I == 1)) {cal_mode <- round(g[c(mode_I), "cal_date"], 0)}
    if (length(mode_I > 1)) {cal_mode <- round(mean(g[c(mode_I), "cal_date"]), 0)}

    # Moments
    probs_mom <- c(round(mean(g[, "cal_probs"], na.rm = TRUE), 5), round(var(g[, "cal_probs"], na.rm = TRUE), 7),
                   round(moments::skewness(g[, "cal_probs"], na.rm = TRUE), 1), round(moments::kurtosis(g[, "cal_probs"], na.rm = TRUE),1))

    # Store estimates
    repo_cal[i_date,] <- c(cal_w_mean, cal_w_sd, calQ_med, calQ_lowerCI, calQ_upperCI, cal_mode, probs_mom)

    #   # fill calibration plot in pdf object
    #   if (isTRUE(plot_save)) {
    #     plot(calRaw, main = chrono_data$id, cex.main = 2, col.main = "blue", col.lab = "red",
    #         ylab = "Non-calibrated radiocarbon date (years)", xlab = "Calibrated radiocarbon date (calendar years)")
    #     abline(v = cal_w_mean, col = "red", lwd = 2)
    #     abline(v = cal_w_mean - cal_w_sd, col = "red", lwd = 1, lty = "dashed")
    #     abline(v = cal_w_mean + cal_w_sd, col = "red", lwd = 1, lty = "dashed")
    #     legend("topright", col = "red", legend = c("weighted mean", cal_ci), lwd = c(2, 1), lty = c("solid", "dashed"))
    #   }
  }
  # # save calibration plot in pdf object to working directory
  # if (isTRUE(plot_save)) {
  #   dev.off()
  # }

  #number of modes at +/- 1SD
  indSD1 <- which(grepl("OneSigma", colnames(calib)) == "TRUE")
  if (length(indSD1) > 1) {modality_SD1 <- apply(calib[, indSD1], 1, function(x) length(which(x != "NA to NA")))} #number of non-NAs
  if (length(indSD1) == 1) {modality_SD1 <- 1} #avoid error when there is 1 mode at +/- 1SD

  #number of modes at +/- 2SD
  indSD2 <- which(grepl("TwoSigma", colnames(calib)) == "TRUE")
  if (length(indSD2) > 1) {modality_SD2 <- apply(calib[, indSD2], 1, function(x) length(which(x != "NA to NA")))} #number of non-NAs
  if (length(indSD2) == 1) {modality_SD2 <- 1} #avoid error when there is 1 mode at +/- 2SD

  modality_SD2 <- modality_SD1 + modality_SD2

  message("\n")
  if (cal_ci == "sd1") {
    message("Intercepts include a calibrated mode, a calibrated weighted average with a 1SD error and a calibrated median with 68.3 inter-quartile ranges\n")
  }
  if (cal_ci == "sd2") {
    message("Intercepts include a calibrated mode, a calibrated weighted average with a 2SD error and a calibrated median with 95.4 inter-quartile ranges\n")
  }
  message("Destriptive stats of calibration computed\n")

  #close progress bar
  close(pb)

  ##############################################
  ## STORING OUTPUTS
  message("\n")
  message("Storing calibration outputs\n")

  #multi-modality output
  repo_modal <- data.table::data.table(chrono_data[, "id"], modality_SD1, modality_SD2, calib[, -c(1,2)], chrono_data[, "curve"])
  colnames(repo_modal)[c(1, dim(repo_modal)[2])] <- c("id", "cal_curve")

  #moments and intercepts output
  if (dim(chrono_data)[1] > 1) {
    repo_mom <- data.table::data.table(chrono_data[, "id"], repo_cal[, c(7:10)], chrono_data[, "curve"])
    repo_cal <- data.table::data.table(chrono_data[, "id"], chrono_data[, "age"], chrono_data[, "sd"], chrono_data[, "curve"], repo_cal[, c(1:6)])
  }
  if (dim(chrono_data)[1] == 1) { #re-arrange if only 1 14C date be calibrated
    repo_mom <- data.table::data.table(t(c(chrono_data[, "id"], repo_cal[, c(7:10)], chrono_data[, "curve"])))
    repo_cal <- data.table::data.table(t(c(chrono_data[, "id"], chrono_data[, "age"], chrono_data[, "sd"], chrono_data[, "curve"], repo_cal[, c(1:6)])))
  }

  #rename columns in moments and intercepts output
  colnames(repo_mom) <- c("id", "pdf_mean", "pdf_var", "pdf_skewness", "pdf_kurtosis", "cal_curve")
  colnames(repo_cal) <- c("id", "uncal_date", "uncal_error", "cal_curve", "calW_mean", "calW_sd", "calQ_median", "calQ_lowerCI", "calQ_upperCI", "cal_mode")

  message("Calibration outputs stored\n")

  #save outputs
  if (cal_save == "intercepts") { #save calibrated dates
    data.table::fwrite(x = repo_cal, file = "outputCALIBRATION_intercepts", col.names = TRUE)
    message("Calibrated intercepts saved to working directory (cal_save = TRUE) with file name = output_CALIBRATION_intercepts\n")
  }
  if (cal_save == "full") { #save calibrated dates
    data.table::fwrite(x = data.table::data.table(repo_cal, calib[, -c(1,2)], repo_mom[, c(2:5)]), file = "outputCALIBRATION_full", col.names = TRUE)
    message("Full calibration saved to working directory (cal_save = TRUE) with file name = output_CALIBRATION\n")
  }
  if (cal_save == "none") { #avoid save calibrated dates
    message("No calibration results saved to working directory (cal_save = FALSE)\n")
  }

  rm(indSD1, indSD2)
  rm(cal_date, cal_probs, calQ_med, cal_mode, calQ_lowerCI, calQ_upperCI, calRaw, cal_w_mean, cal_w_sd, g, i_date, mode_I, probs_mom)
  rm(calib, chrono_data, modality_SD1, modality_SD2)

  message("All calibrated times reported in calendar years Before Present (Present = year 1950)\n")

  ##############################################
  ## PRINT RESULTS TO CONSOLE AND CLEAN R ENVIRONMENT
  out <- list(cal_arguments = calendar_args, cal_moments = repo_mom, cal_modality = repo_modal, cal_intercepts = repo_cal)
  return(out)
  rm(list = ls(all.names = TRUE)) #clean R environment
}

## END
##############################################
