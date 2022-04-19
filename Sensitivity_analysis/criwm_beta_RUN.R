criwm <- function(fossil_data, #input data
                  signor_lipps = "ext", biased = TRUE, alpha = 0.05, resample = 10000, #criwm parameters
                  radiocarbon = "all", calibra = TRUE, cal_curve = "intcal20", upper14C = 55000, #calibration parameters
                  cal_save = FALSE, resample_save = FALSE, criwm_save = FALSE) #saving parameters
{
  
  ##############################################
  ## ABORT SESSION IF PARAMETERS ENTERED INCORRECTLY
  if (!signor_lipps %in% c("ext", "inv")) {
    stop("Parameter signor_lipps must be ext or inv, please check spelling\nComputation aborted\n")
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
  criwm_args <- match.call()
  argsNames <- names(as.list(criwm_args))[-1] #get names of arguments
  argsValues <- as.character(criwm_args)[-1] #get values selected per argument
  criwm_args <- data.frame(argsNames, argsValues) #dataframe arguments and their values
  colnames(criwm_args) <- c("argument", "value") #rename columns
  rm(argsNames, argsValues) #remove objects
  
  ##############################################
  ## PREPARE DATA
  
  # LOAD data
  message("\n")
  if (signor_lipps == "ext") {message(paste("########CRIWM-extinction analysis started / Time = ", date(), "\n"))}
  if (signor_lipps == "inv") {message(paste("########CRIWM-invasion analysis started / Time = ", date(), "\n"))}
  message("Loading data\n")
  message("Input data assumed to be a text file, without row names or column headings, column 1 = dates in years Before Present (integers), column 2 = errors in years (integers)\n")
  fossil_data <- data.frame(data.table::fread(file = fossil_data))
  colnames(fossil_data) <- c("age", "sd") #rename columns
  dateN <- dim(fossil_data)[1] #number of dates in data
  dat_uncal <- fossil_data #create object for output file
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
    fossil_data <- data.frame(as.numeric(fossil_data[, 1]), as.numeric(fossil_data[, 2])) #dataframe repository and coerce mean and sd to numeric
    colnames(fossil_data) <- c("age", "sd") #rename columns
    dat_cal <- data.table::data.table(fossil_data)
  }
  #subset non-calibrated radiocarbon dates
  if (isTRUE(calibra) & radiocarbon[1] != 0) { #radiocarbon dates present 
    #subset radiocarbon dates if list includes dates obtained through chronological methods
    if (is.character(radiocarbon[1]) == TRUE) { #all dates are non-calibrated radiocarbon
      fossil_data <- fossil_data
    } 
    if (is.numeric(radiocarbon[1]) == TRUE) { #mixed radiocarbon and non-radiocarbon dates prsent
      dat1 <- fossil_data[-radiocarbon, ] #object contains calibrated dates
      fossil_data <- fossil_data[radiocarbon, ] #object contains non-calibrated dates
    } 
    
    #calibrate
    fossil_data$curve <- cal_curve
    calRaw <- rcarbon::calibrate(fossil_data[, 1], errors = fossil_data[, 2], calCurves = fossil_data[, 3], timeRange = c(upper14C, 0), F14C = TRUE, eps = 0)
    fossil_data <- matrix(NA, nrow = nrow(fossil_data), ncol = 2) #repository of calibrated point estimates of radiocarbon dates
    fossil_pdf <- list() #repository of calibrated point estimates of radiocarbon dates
    for (i_date in 1:length(calRaw)) {
      cal_date <- as.numeric(data.frame(calRaw$grids[i_date])[, 1]) #get calibrated values
      cal_probs <- as.numeric(data.frame(calRaw$grids[i_date])[, 2]) #get PDF
      cal_w_mean <- round(sum(cal_date * cal_probs), 0) #estimate calibrated weighted mean
      cal_w_sd <- round(sqrt(sum(((cal_date - cal_w_mean) ^ 2) * cal_probs)), 0) #estimate calibrated weighted sd
      fossil_data[i_date, ] <- c(cal_w_mean, cal_w_sd) #store cal weighted mean and sd
      fossil_pdf[[i_date]] <- data.frame(calRaw$grids[i_date])
      # message(paste("date", i_date))
    }
    fossil_data <- data.frame(as.numeric(fossil_data[, 1]), as.numeric(fossil_data[, 2])) #dataframe repository and coerce mean and sd to numeric
    colnames(fossil_data) <- c("age", "sd") #rename columns
    
    #rbind calibrated radiocarbon dates with non-radiocarbon dates, re-order by row number and assign data to data.table object
    if (is.numeric(radiocarbon[1]) == TRUE & radiocarbon[1] != 0) {
      rownames(fossil_data) <- radiocarbon #row number for radiocarbon dates
      fossil_data <- data.frame(rbind(fossil_data, dat1)) #rbind radiocarbon and non-radiocarbon dates
      fossil_data <- data.table::data.table(fossil_data[order(as.numeric(rownames(fossil_data))), ]) #get original row order
      rm(dat1) #remove object
    }
    #create object for output file when at least one date is radiocarbon
    dat_cal <- data.table::data.table(fossil_data)
    rm(calRaw, cal_date, cal_probs, cal_w_mean, cal_w_sd) #clean space from R environment
  }
  
  #assign data to data.table object when all input dates are calibrated
  if (isFALSE(calibra)) {
    fossil_data <- data.table::data.table(fossil_data)
    colnames(fossil_data) <- c("age", "sd") #rename columns
    dat_cal <- fossil_data #create object for output file
  }
  
  #assign data to data.table object when no input dates are calibrated or radiocarbon
  if (is.character(radiocarbon[1]) == TRUE | radiocarbon[1] == 0) {
    fossil_data <- data.table::data.table(fossil_data)
  }
  
  # SORT data
  if (signor_lipps == "ext") {fossil_data <- fossil_data[order(fossil_data[, 1], decreasing = FALSE), c(1:2)]} #extinction
  if (signor_lipps == "inv") {fossil_data <- fossil_data[order(fossil_data[, 1], decreasing = TRUE), c(1:2)]} #invasion
  
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
  ## MESSAGE CRIWM PARAMETERS TO CONSOLE
  message("########CRIWM re-sampling started\n")
  if(radiocarbon[1] == 0 & isTRUE(calibra)) {
    message("Calibration requested (calibra == TRUE) but all input dates are calibrated (radiocarbon == 0), therefore re-sampling done from a normal distribution and CRIWM and GRIWM computations and results are identical\n")
  }
  if(radiocarbon[1] == 0 & isFALSE(calibra)) {
    message("Calibration not requested (calibra == FALSE) because all input dates are non-radiocarbon (radiocarbon == 0), therefore re-sampling done from a normal distribution and CRIWM and GRIWM computations and results are identical\n")
  }
  if(isFALSE(calibra) & is.character(radiocarbon[1])) {
    message("Calibration not requested (calibra == FALSE), re-sampling done from a normal distribution and CRIWM and GRIWM computations and results are identical\n")
  }
  if(isFALSE(calibra) & is.numeric(radiocarbon) & radiocarbon[1] != 0) {
    message("Calibration not requested (calibra == FALSE), re-sampling done from a normal distribution and CRIWM and GRIWM computations and results are identical\n")
  }  
  if (signor_lipps == "ext") {
    message(paste("Extinction time (signor_lipps == 'ext') from ", resample, " re-sampled estimates (resample = ", resample, ")", sep = ""), "\n")
  } 
  if (signor_lipps == "inv") {    
    message(paste("Invasion time (signor_lipps == 'inv') from ", resample, " re-sampled estimates (resample = ", resample, ")", sep = ""), "\n")
  }
  if (isTRUE(biased)) {
    message(paste("Biased estimator applied (biased = TRUE) for alpha = ", alpha, " (alpha = ", alpha, ")" , sep = ""), "\n")
  } else {
    message("Unbiased estimator applied (biased = FALSE)\n")
  }
  
  ## COMPUTE CRIWM
  
  # VECTORIZE data, set criwm parameters, define output repositories
  datF <- c(fossil_data[, 1])$age #get dates
  sd.vec <- c(fossil_data[, 2])$sd #get errors
  rm(fossil_data) #remove data object
  k <- length(datF) #get number of dates
  out <- matrix(0, 1, 3) #create repository of final criwm results 
  T.up.vec <- T.mci.vec <- w.T.mci.vec <- rep(0, resample) #create repository of temporal criwm results
  
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
      
      #all dates are non-radiocarbon = re-sampling from normal distribution whether calibration requested or not requested (CRIWM = GRIWM)
      if(radiocarbon[1] == 0) {
        date.samp[i_resample] <- round(rnorm(1, datF[i_resample], sd.vec[i_resample]))
      }
      #all dates are non-calibrated radiocarbon = re-sampling from PDF (calibration requested)
      if(isTRUE(calibra) & is.character(radiocarbon[1])) {
        #date.samp[i_resample] <- as.numeric(dplyr::sample_n(fossil_pdf[[i_resample]], size = 1, weight = fossil_pdf[[i_resample]][,2], replace = TRUE))[1]
        date.samp[i_resample] <- sample(fossil_pdf[[i_resample]][, 1], size = 1, prob = fossil_pdf[[i_resample]][, 2], replace = TRUE)
        
      }
      #some dates are non-calibrated radiocarbon = re-sampling from PDF (calibration requested)
      #some dates are non-radiocarbon or calibrated radiocarbon = re-sampling from normal distribution (calibration avoided)
      if(isTRUE(calibra) & is.numeric(radiocarbon[1]) & radiocarbon[1] != 0) {
        no14C <- setdiff(c(1:k), radiocarbon) #index (row number) for non-radiocarbon dates
        if(i_resample %in% no14C) { #loop over non-radiocarbon dates or calibrated radiocarbon dates
          date.samp[i_resample] <- round(rnorm(1, datF[i_resample], sd.vec[i_resample]))
        }
        if(i_resample %in% radiocarbon) { #loop over non-calibrated radiocarbon dates
          i_resample_I <- which(radiocarbon == i_resample) #index for PDF across non-calibrated radiocarbon dates only
          #date.samp[i_resample] <- as.numeric(dplyr::sample_n(fossil_pdf[[i_resample_I]], size = 1, weight = fossil_pdf[[i_resample_I]][,2], replace = TRUE))[1]
          date.samp[i_resample] <- sample(fossil_pdf[[i_resample_I]][, 1], size = 1, prob = fossil_pdf[[i_resample_I]][, 2], replace = TRUE)
        }
      }
      #all dates are calibrated radiocarbon = re-sampling from normal distribution (calibration not requested)
      if(isFALSE(calibra) & is.character(radiocarbon[1])) {
        date.samp[i_resample] <- round(rnorm(1, datF[i_resample], sd.vec[i_resample]))
      }
      #all dates are calibrated radiocarbon or non-radiocarbon (re-sampling from normal distribution)
      if(isFALSE(calibra) & is.numeric(radiocarbon[1]) & radiocarbon[1] != 0) {
        date.samp[i_resample] <- round(rnorm(1, datF[i_resample], sd.vec[i_resample]))
      }
    }
    
    #set estimate: extinction versus invasion
    if (signor_lipps == "ext") {date.samp <- sort(date.samp)}
    if (signor_lipps == "inv") {date.samp <- sort(date.samp, decreasing = TRUE)}
    
    #estimate weights per date
    last.diff <- 1 / (date.samp - date.samp[1])[-1]
    weight <- last.diff / last.diff[1]
    if (last.diff[1] == Inf) {
      weight <- last.diff / last.diff[2]
      weight <- weight[-1]
    }
    
    #estimate extinction or invasion parameters
    ldate <- length(date.samp)
    T.mci.lst.vec <- rep(0, ldate - 1)
    for (i_criwm in 1:(ldate-1)) { #loop over ldate = k - 1 dates
      date.it <- date.samp[1:(1 + i_criwm)]
      date.age.it <- date.samp[1:(1 + i_criwm)]
      date.mci.it <- rev(max(date.it) + 1 - date.it)
      k <- length(date.it)
      t.n <- date.mci.it[k]
      n <- k
      T.rng <- t.n - date.mci.it[1]
      #biased versus unbiased estimate
      if (isTRUE(biased)) {
        i = t.n - t.n * log(alpha) / n #CRIWM biased estimate (original)
      }
      if (isFALSE(biased)) {
        i = t.n * (n + 1) / n #CRIWM unbiased estimate
      }
      T.mci.lst.vec[i_criwm] <- max(date.it) + 1 - i
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
  
  # CALCULATE median extinction (or invasion) time and 95% CI
  T.wmci.vec.CI <- quantile(na.omit(w.T.mci.vec), probs = c(0.025, 0.5, 0.975));
  out[1] <- round(T.wmci.vec.CI[1], 0) #2.5% quartile (lower CI bound)
  out[2] <- round(T.wmci.vec.CI[2], 0) #50% quartile (median)
  out[3] <- round(T.wmci.vec.CI[3], 0) #97.5% quartile (upper CI bound)
  out <- data.frame(out) #convert to dataframe
  
  #convert BP = Before Present
  outBCE <- out - 1950 #to BCE = Before Current Era
  outBCY <- out  + as.numeric(format(Sys.Date(), "%Y")) - 1950 #to BCY = Before Current Year
  out <- data.frame(rbind(outBCE, out, outBCY)) #data.frame the three time scales
  
  if (signor_lipps == "ext") {rownames(out) <- c("Extinction_Time_BCE", "Extinction_Time_BP", "Extinction_Time_BCY")}
  if (signor_lipps == "inv") {rownames(out) <- c("Invasion_Time_BCE", "Invasion_Time_BP", "Invasion_Time_BCY")}
  colnames(out) <- c("2.5CI", "median", "97.5CI")
  
  # REMOVE CRIWM loop inputs from console
  rm(datF, sd.vec, k, date.samp, last.diff, weight, ldate, date.it, date.age.it, date.mci.it, t.n, n, T.rng, outBCE, outBCY)
  
  # MESSAGE end of computations
  message("\n")
  message("########CRIWM re-sampling ended\n")  
  if (signor_lipps == "ext") {message("CRIWM extinction time reported in calendar years BP, BCE and BCY\n")}
  if (signor_lipps == "inv") {message("CRIWM invasion time reported in calendar years BP, BCE and BCY\n")}
  message("BP = Before Present (Present = year 1950)\n")
  message("BCE = Before Current Era (Current Era starts at year 0)\n")
  message(paste("BCY = Before Current Year (Current Year = ", as.numeric(format(Sys.Date(), "%Y")), ")", sep = ""), "\n")
  if (signor_lipps == "ext") {message("Note = extinction time in BP or BCE will be negative if earlier than 1950 or 0, respectively\n")}
  if (signor_lipps == "inv") {message("Note = invasion time in BP or BCE will be negative if earlier than 1950 or 0, respectively\n")}
  
  ##############################################
  ## SAVE CRIWM MEDIAN EXTINCTION (OR INVASION) TIME AND 95% CI
  if (criwm_save == TRUE & signor_lipps == "ext") {
    data.table::fwrite(x = out, file = "outputCRIWM_medianCI_Extinction", col.names = TRUE)
    message("Median extinction time and 95% confidence interval saved to working directory (criwm_save = TRUE) with file name = outputCRIWM_medianCI_Extinction\n")
  }
  if (criwm_save == TRUE & signor_lipps == "inv") {
    data.table::fwrite(x = out, file = "outputCRIWM_medianCI_Invasion", col.names = TRUE)
    message("Median invasion time and 95% confidence interval saved to working directory (criwm_save = TRUE) with file name = outputCRIWM_medianCI_Invasion\n")
  }
  
  ##############################################
  ## SAVE RE-SAMPLED EXTINCTION (OR INVASION) TIMES    
  if (signor_lipps == "ext" & isTRUE(resample_save)) {
    data.table::fwrite(x = data.frame(w.T.mci.vec), file = "outputCRIWM_resampled_Extinction", col.names = FALSE)
    message("Re-sampled extinction times saved to working directory (resample_save = TRUE) with file name = outputCRIWM_resampled_Extinction\n")
  }
  if (signor_lipps == "inv" & isTRUE(resample_save)) {
    data.table::fwrite(x = data.frame(w.T.mci.vec), file = "outputCRIWM_resampled_Invasion", col.names = FALSE)
    message("Re-sampled invasion times saved to working directory (resample_save = TRUE) with file name = outputCRIWM_resampled_Invasion\n")
  }
  # REMOVE CRIWM loop outputs from console
  rm(T.up.vec, T.mci.vec, w.T.mci.vec, T.wmci.vec.CI, T.mci.lst.vec) #remove CRIWM outputs from console
  
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
    calibration <- data.table::data.table(uncal.age = dat_uncal[, 1], uncal.sd = dat_uncal[, 2], cal = dat_cal[, 1], cal = dat_cal[, 2], type = chrono_type)
    calibration <- calibration[order(calibration$cal.age),] #order from younger to older dates
    if (isTRUE(cal_save) & radiocarbon[1] != 0) { #save non-calibrated and calibrated dates
      data.table::fwrite(x = data.frame(calibration), file = "outputCRIWM_calibration", col.names = TRUE)
      message("Dates in years and calendar years saved to working directory (cal_save = TRUE) with file name = outputCRIWM_calibration\n")
    }
    if (isTRUE(cal_save) & radiocarbon[1] == 0) { #save calibrated dates
      data.table::fwrite(x = data.frame(calibration)[,c(3:5)], file = "outputCRIWM_calibration", col.names = TRUE)
      message("All input dates were in calendar years (radiocarbon = 0) and have been saved to working directory (cal_save = TRUE) with file name = outputCRIWM_calibration\n")
    }
  }
  if (isFALSE(calibra)) { #calibration not requested
    calibration <- data.table::data.table(cal = dat_cal[, 1], cal = dat_cal[, 2], type = chrono_type)
    calibration <- calibration[order(calibration$cal.age),] #order from younger to older dates
    if (isTRUE(cal_save)) { #save calibrated dates
      data.table::fwrite(x = data.frame(calibration), file = "outputCRIWM_calibration", col.names = TRUE)
      message("All input dates are in calendar years (calibra = FALSE) and have been saved to working directory (cal_save = TRUE) with file name = outputCRIWM_calibration\n")
    }
  }
  
  ##############################################
  ## PRINT RESULTS TO CONSOLE AND CLEAN R ENVIRONMENT
  if (signor_lipps == "ext") {message(paste("########CRIWM-extinction analysis ended / Time = ", date(), "\n"))}
  if (signor_lipps == "inv") {message(paste("########CRIWM-invasion analysis ended / Time = ", date(), "\n"))}
  out <- list(criwm_arguments = criwm_args, calibration_BP = calibration, criwm = out)
  return(out)
  rm(list = ls(all = TRUE)) #clean R environment
}

## END
##############################################