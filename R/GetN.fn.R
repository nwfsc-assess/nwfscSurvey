#' Calculate effN input sample sizes
#' based on Stewart, I.J. and O.S. Hamel 2014.
#' Bootstrapping of sample size for legth- or age-composition data used in stock assessment
#' Canadian Journal of Fishery and Aquatic Science, 71: 581-588.
#'
#' @param dir directory
#' @param dat data frame to be passed in
#' @param type specify whether doing "length" or "age". Used to read associatied excel sheets
#' @param species species specific value to determine the number of unique samples per tow (flatfish, shelfrock, sloperock, thorny, others, all)
#' @param printfolder name of the folder to create and save files. Location will be paste0(dir, printfolder)
#' @param output default = NULL will return only a vector of samples sizes, summary will return a table of observations
#' by year and sex
#' @param verbose opt to print out message statements
#'
#' @author Chantel Wetzel
#' @export

GetN.fn <- function(dir = NULL, dat, type, species = NULL, printfolder = "forSS", output = NULL, verbose = TRUE) {
  n.unq <- NA
  if (species == "flatfish") {
    n.unq <- 3.09
  }
  if (species == "shelfrock") {
    n.unq <- 2.43
  }
  if (species == "sloperock") {
    n.unq <- 2.43
  }
  if (species == "thorny") {
    n.unq <- 6.91
  }
  if (species == "others") {
    n.unq <- 2.38
  }
  if (species == "all") {
    n.unq <- 2.73
  }
  if (is.na(n.unq)) {
    if (verbose) {
      message("\nThe species input does not match one of the following options; flatfish, shelfrock, sloperock, thorny, others, or all\n")
    }
  }

  if (verbose) {
    message("\nThe effN sample size is calculated using the ", species, " multiplier of ", n.unq, ". This number is multiplied by the number of tows in each year.\n")
  }


  if (type == "length") {
    temp <- dat[!is.na(dat$Length_cm), ]
    nSamp <- table(temp$Year, !duplicated(as.character(temp$Trawl_id)))[, "TRUE"]
  }

  if (type == "age") {
    temp <- dat[!is.na(dat$Age), ]
    nSamp <- table(temp$Year, !duplicated(as.character(temp$Trawl_id)))[, "TRUE"]
  }

  if (length(nSamp) == 1) {
    yr <- unique(temp$Year)
  }
  if (length(nSamp) > 1) {
    yr <- as.numeric(names(nSamp))
  }


  n <- floor(n.unq * nSamp)
  fish <- table(temp$Year, temp$Sex)

  if ("F" %in% colnames(fish)) {
    female <- as.numeric(fish[, "F"])
  } else {
    female <- rep(0, dim(fish)[1])
  }

  if ("M" %in% colnames(fish)) {
    male <- as.numeric(fish[, "M"])
  } else {
    male <- rep(0, dim(fish)[1])
  }

  if ("U" %in% colnames(fish)) {
    unsex <- as.numeric(fish[, "U"])
  } else {
    unsex <- rep(0, dim(fish)[1])
  }

  # Add check to cap input N to not be greater than total fish sampled
  ind <- n > (male + female + unsex)
  n[ind] <- male[ind] + female[ind] + unsex[ind]

  if (sum(ind) > 0) {
    if (verbose) {
      cat("\nInput sample size exceded the number of fish for", yr[ind], "and has been capped equal to number of fish.\n")
    }
  }

  samples <- data.frame(
    Year = yr,
    Tows = nSamp,
    All_Fish = female + male + unsex,
    Sexed_Fish = female + male,
    Unsexed_Fish = unsex,
    Sample_Size = n
  )

  # save output as a csv
  if (!is.null(dir)) {
    plotdir <- file.path(dir, printfolder)
    plotdir.isdir <- file.info(plotdir)$isdir
    if (is.na(plotdir.isdir) | !plotdir.isdir) {
      dir.create(plotdir)
    }
    write.csv(samples, file = file.path(plotdir, paste0(type, "_SampleSize.csv", sep = "")), row.names = FALSE)
  }
  if (is.null(output)) {
    return(n)
  } else {
    return(samples)
  }
}
