calendar <- function(fossil_data,
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
  fossil_data <- data.frame(data.table::fread(file = fossil_data))
  dateN <- dim(fossil_data)[1] #number of dates in data
  
  #rename columns
  if (dim(fossil_data)[2] == 4) { #four columns
    message("Data contains 4 columns (ages, errors, identifiers, calibration curves)\n")
    colnames(fossil_data) <- c("age", "sd", "id", "curve")
  }
  if (dim(fossil_data)[2] == 3) { #three columns
    message("Data contains 3 columns (ages, errors, identifiers)\n")
    fossil_data$curve <- cal_curve #add column with names of calibration curve
    colnames(fossil_data) <- c("age", "sd", "id", "curve")
  } 
  
  # ABORT calibration if identifiers (column 3) are not unique for each date
  codN <- length(which(table(fossil_data$id) == 1)) #frequency of duplicate ids
  if (codN - dateN != 0) {
    stop("Calibration aborted\nSome identifiers (column 3) have been assigned to more than one date\nTo proceed with calibration, please provide unique identifiers to each date and reload the input data\n")
  }
  message(paste(dateN, " radiocarbon dates loaded\n", sep = ""))
  rm(codN, dateN)
  
  ##############################################
  ## CALIBRATE dates
  message("\n")
  message("Calibrating dates entered in years Before Present\n")
  
  # CALIBRATE
  #message number of calibration curves
  if (length(table(fossil_data$curve)) > 1) {
    message("Different dates calibrated with different calibration curves\n")
  }
  if (length(table(fossil_data$curve)) == 1) {
    message("All dates calibrated with the same calibration curve = ", row.names(table(fossil_data$curve)), "\n")
  }
  #timeRange = sets temporal window
  #F14C = calibration in F14C domain (radiocarbon concentration) rather than temporal domain (radiocarbon age)
  #eps = range of probability density function considered (eps = 0 captures full PDF)
  calRaw <- rcarbon::calibrate(x = fossil_data[, "age"], errors = fossil_data[, "sd"], ids = fossil_data[, "id"], calCurves = fossil_data[, "curve"], timeRange = c(upper14C,0), F14C = TRUE, eps = 0)
  calib <- summary(calRaw) #summary of calibration output
  
  message("Calibration in calendar years Before Present completed\n")
  
  ##############################################
  ## ESTIMATE calibrated intercepts
  message("\n")
  message("Computing destriptive stats of calibration\n")
  
  #set progress bar
  pb <- txtProgressBar(min = 1, max = length(calRaw), style = 3)
  
  # Weighted average and standard deviation
  repo_cal <- matrix(NA, nrow = length(calRaw), ncol = 10) #ncol = number of metrics computed
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
    if (length(mode_I == 1)) {cal_mode <- g[c(mode_I), "cal_date"]}
    if (length(mode_I > 1)) {cal_mode <- mean(g[c(mode_I), "cal_date"])}
    
    # Moments
    probs_mom <- c(round(mean(g[, "cal_probs"], na.rm = TRUE), 5), round(var(g[, "cal_probs"], na.rm = TRUE), 7),
                   round(moments::skewness(g[, "cal_probs"], na.rm = TRUE), 1), round(moments::kurtosis(g[, "cal_probs"], na.rm = TRUE),1))
    
    # Store estimates
    repo_cal[i_date,] <- c(cal_w_mean, cal_w_sd, calQ_med, calQ_lowerCI, calQ_upperCI, cal_mode, probs_mom)
  }
  
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
  repo_modal <- data.table::data.table(fossil_data[, "id"], modality_SD1, modality_SD2, calib[, -c(1,2)], fossil_data[, "curve"])
  colnames(repo_modal)[c(1, dim(repo_modal)[2])] <- c("id", "cal_curve")
  
  #moments output
  repo_mom <- data.table::data.table(fossil_data[, "id"], repo_cal[, c(7:10)], fossil_data[, "curve"])
  colnames(repo_mom) <- c("id", 
                          "pdf_mean", "pdf_var", "pdf_skewness", "pdf_kurtosis", 
                          "cal_curve")
  
  #intercepts output
  repo_cal <- data.table::data.table(fossil_data[, "id"], fossil_data[, "age"], fossil_data[, "sd"], fossil_data[, "curve"], repo_cal[, c(1:6)])
  colnames(repo_cal) <- c("id", "uncal_date", "uncal_error", "cal_curve",
                          "calW_mean", "calW_sd", "calQ_median", "calQ_lowerCI", "calQ_upperCI", "cal_mode")
  
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
  rm(calib, fossil_data, modality_SD1, modality_SD2)
  
  message("All calibrated times reported in calendar years Before Present (Present = year 1950)\n")
  
  ##############################################
  ## PRINT RESULTS TO CONSOLE AND CLEAN R ENVIRONMENT
  out <- list(cal_arguments = calendar_args, cal_moments = repo_mom, cal_modality = repo_modal, cal_intercepts = repo_cal)
  return(out)
  rm(list = ls(all = TRUE)) #clean R environment
}

## END
##############################################