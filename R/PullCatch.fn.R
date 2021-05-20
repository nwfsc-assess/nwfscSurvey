#' Pull catch data from the NWFSC data warehouse
#'
#' Pull catch data from the
#' [NWFSC data warehouse](https://www.webapps.nwfsc.noaa.gov/data)
#' for a single species or all observed species, where the latter is specified
#' by leaving both `Name = NULL` and `SciName = NULL`.
#'
#' @details
#' The data available in the warehouse are cleaned pior to being downloaded
#' with the intent that they provide the best available information for use
#' in an index-standardization procedure. The removed samples may be of use
#' to others with a less-restrictive goal than producing an index of abundance.
#' For example, life-stage samples are excluded because they are not collected
#' using the same protocols as standard samples.
#' To download all data, we currently recommend going to the
#' [NWFSC data warehouse](https://www.webapps.nwfsc.noaa.gov/data)
#' and using the csv link to extract data for a single species at a time.
#' In the future, we hope to add functionality to this package such that
#' downloading all data can be done easily within this function.
#' See [Issue #43](https://github.com/nwfsc-assess/nwfscSurvey/issues/43)
#' for more information.
#'
#' @param Name A character entry with the desired common name of the
#' species you want to pull data for from the data warehouse.
#' Use a vector of names if you want information for more than one species or
#' if the desired species is included in the database using more than one name,
#' e.g., vermilion rockfish (see the example below).
#' Use the `SciName` argument if you know the latin name.
#' @param SciName A character entry with the desired scientific name of the
#' species you want to pull data for from the data warehouse.
#' Use a vector of names if you want information for more than one species or
#' if the desired species is included in the database using more than one name,
#' e.g., vermilion rockfish (see the example below).
#' Use the `Name` argument if you know the common name.
#' @param YearRange An integer vector of length two with the
#' range of years to pull data for.
#' @param SurveyName A character entry from one of the following options that
#' specifies which survey to pull the data for:
#'   * Triennial,
#'   * AFSC.Slope,
#'   * NWFSC.Combo,
#'   * NWFSC.Slope,
#'   * NWFSC.Shelf,
#'   * NWFSC.Hypoxia,
#'   * NWFSC.Santa.Barb.Basin,
#'   * NWFSC.Shelf.Rockfish (NWFSC.Hook.Line but both are not working),
#'   * NWFSC.Video.
#' Currently, you must pull data one survey at a time, though we are working on
#' allowing for a vector of survey names and
#' `NWFSC.Shelf.Rockfish` and `NWFSC.Hook.Line` are not supported.
#' The default of `NULL` is a placeholder that must be replaced with an entry.
#' @param SaveFile A logical value specifying whether or not the the data should
#' be saved to a file in `Dir`. Must change from the default of `FALSE` to save a file.
#' @param Dir If `SaveFile = TRUE`, then one must specify the directory where you want
#' the resulting file to be saved. The directory where the file should be saved.
#' The name of the file within `Dir` will start with Catch_ and end with .rda.
#' @param verbose A logical value specifying whether or not to pring out
#' message statements to the screen while using this function.
#' The default is to print all messages, i.e., `verbose = TRUE`.
#'
#' @author Chantel Wetzel based on code by John Wallace
#' @export
#' @md
#'
#' @import jsonlite
#' @import chron
#' @importFrom stringr str_replace_all
#' @importFrom dplyr left_join rename
#'
#' @examples
#' \dontrun{
#' # SurveyName is only arg that has to be specified
#' dat <- PullCatch.fn(SurveyName = "NWFSC.Combo")
#'
#' # Example with specified common name
#' catch_dat <- PullCatch.fn(Name = "vermilion rockfish",
#' SurveyName = "NWFSC.Combo")
#'
#' # Example with specified scientific name
#' catch_dat <- PullCatch.fn(SciName = "Eopsetta jordani",
#' SurveyName = "NWFSC.Combo")
#'
#' # Example with multiple names
#' catch_dat <- PullBio.fn(Name = c("vermilion rockfish",
#' "vermilion and sunset rockfish"), SurveyName = "NWFSC.Combo")
#'
# catch_dat <- PullCatch.fn(SciName = c("Sebastes miniatus",
# "Sebastes sp. (crocotulus)","Sebastes sp. (miniatus / crocotulus)"),
# SurveyName = "NWFSC.Combo")
#' }
#'
PullCatch.fn <- function(Name = NULL, SciName = NULL, YearRange = c(1980, 5000), SurveyName = NULL, SaveFile = FALSE, Dir = NULL, verbose = TRUE) {
  if (SurveyName %in% c("NWFSC.Shelf.Rockfish", "NWFSC.Hook.Line")) {
    stop("The catch pull currently does not work for hook & line data.",
      "\nPull directly from the warehouse https://www.webapp.nwfsc.noaa.gov/data")
  }

  if (SaveFile) {
    if (is.null(Dir)) {
      stop("The Dir input needs to be specified in order to save output file.")
    }
    if (!file.exists(Dir)) {
      stop(
        "The Dir argument leads to a location",
        ",\ni.e., ", Dir, ", that doesn't exist."
      )
    }
  }

  if (is.null(Name)) {
    var.name <- "scientific_name"
    Species <- SciName
    new.name <- "Scientific_name"
    outName <- Name
  }
  if (is.null(SciName)) {
    var.name <- "common_name"
    Species <- Name
    new.name <- "Common_name"
    outName <- SciName
  }
  if (is.null(SciName) & is.null(Name)) {
    var.name <- "common_name"
    Species <- "pull all"
    new.name <- "Common_name"
  } # stop("Need to specifiy Name or SciName to pull data!")}

  # Survey options available in the data warehouse
  surveys <- createMatrix()

  # Check the input survey name against available options
  if (!SurveyName %in% surveys[, 1]) {
    stop(
      "The SurveyName argument does not match one of the available options:\n",
      paste(surveys[, 1], collapse = "\n")
    )
  }

  # Find the long project name to extract data from the warehouse
  for (i in 1:dim(surveys)[1]) {
    if (SurveyName == surveys[i, 1]) {
      project <- surveys[i, 2]
      projectShort <- surveys[i, 1]
    }
  }

  if (length(YearRange) == 1) {
    YearRange <- c(YearRange, YearRange)
  }


  # Pull data for the specific species for the following variables
  Vars <- c(
    var.name, "year", "subsample_count", "subsample_wt_kg", "project", "cpue_kg_per_ha_der",
    "total_catch_numbers", "total_catch_wt_kg", "vessel", "tow", "operation_dim$legacy_performance_code",
    "statistical_partition_dim$statistical_partition_type", "station_code"
  )

  Vars.short <- c(
    var.name, "year", "subsample_count", "subsample_wt_kg", "project", "cpue_kg_per_ha_der",
    "total_catch_numbers", "total_catch_wt_kg", "vessel", "tow", "station_code"
  )

  # symbols here are generally: %22 = ", %2C = ",", %20 = " "
  species_str <- paste0("%22",stringr::str_replace_all(Species[1]," ","%20"),"%22")
  if(length(Species) > 1) {
    for(i in 2:length(Species)) {
      species_str <- paste0(species_str, "%2C", paste0("%22",stringr::str_replace_all(Species[i]," ","%20"),"%22"))
    }
  }

  UrlText <- paste0(
    "https://www.webapps.nwfsc.noaa.gov/data/api/v1/source/trawl.catch_fact/selection.json?filters=project=", paste(strsplit(project, " ")[[1]], collapse = "%20"), ",",
    "station_invalid=0,",
    "performance=Satisfactory,", "depth_ftm>=30,depth_ftm<=700,",
    "field_identified_taxonomy_dim$", var.name, "|=[", species_str,"]",
    ",date_dim$year>=", YearRange[1], ",date_dim$year<=", YearRange[2],
    "&variables=", paste0(Vars, collapse = ",")
  )

  if (Species[1] == "pull all") {
    UrlText <- paste0(
      "https://www.webapps.nwfsc.noaa.gov/data/api/v1/source/trawl.catch_fact/selection.json?filters=project=", paste(strsplit(project, " ")[[1]], collapse = "%20"), ",",
      "station_invalid=0,",
      "performance=Satisfactory,", "depth_ftm>=30,depth_ftm<=700,",
      "date_dim$year>=", YearRange[1], ",date_dim$year<=", YearRange[2],
      "&variables=", paste0(Vars, collapse = ",")
    )
  }

  if (verbose) {
    message("Pulling catch data. This can take up to ~ 30 seconds (or more).")
  }
  # Pull data from the warehouse
  DataPull <- try(jsonlite::fromJSON(UrlText))
  if (!is.data.frame(DataPull)) {
    stop(cat("\nNo data returned by the warehouse for the filters given.\n Make sure the year range is correct for the project selected and the input name is correct,\n otherwise there may be no data for this species from this project.\n"))
  }

  # Remove water hauls
  fix <- is.na(DataPull[, "operation_dim$legacy_performance_code"])
  if (sum(fix) > 0) {
    DataPull[fix, "operation_dim$legacy_performance_code"] <- -999
  }
  # Whether values are NA or "NA" varies based on the presence of "Life Stage" samples
  if (sum(is.na(DataPull[, "statistical_partition_dim$statistical_partition_type"])) != dim(DataPull)[1]) {
    keep <- DataPull[, "statistical_partition_dim$statistical_partition_type"] == "NA"
    DataPull <- DataPull[keep, ]
  }

  keep <- DataPull[, "operation_dim$legacy_performance_code"] != 8
  DataPull <- DataPull[keep, ]
  DataPull <- DataPull[, Vars.short]

  Data <- dplyr::rename(DataPull,
    Year = year, Subsample_count = subsample_count,
    Subsample_wt_kg = subsample_wt_kg, Project = project,
    CPUE_kg_per_ha = cpue_kg_per_ha_der, Subsample_count = subsample_count,
    Subsample_wt_kg = subsample_wt_kg, Vessel = vessel, Tow = tow
  )

  names(Data)[which(names(Data) == "scientific_name")] <- "Scientific_name"
  names(Data)[which(names(Data) == "common_name")] <- "Common_name"

  # Pull all tow data (includes tows where the species was not observed)
  Vars <- c("project", "year", "vessel", "pass", "tow", "datetime_utc_iso", "depth_m", "longitude_dd", "latitude_dd", "area_swept_ha_der", "trawl_id", "operation_dim$legacy_performance_code", "station_code")
  Vars.short <- c("project", "year", "vessel", "pass", "tow", "datetime_utc_iso", "depth_m", "longitude_dd", "latitude_dd", "area_swept_ha_der", "trawl_id", "station_code")

  UrlText <- paste0(
    "https://www.webapps.nwfsc.noaa.gov/data/api/v1/source/trawl.operation_haul_fact/selection.json?filters=project=", paste(strsplit(project, " ")[[1]], collapse = "%20"), ",",
    "station_invalid=0,",
    "performance=Satisfactory,",
    "depth_ftm>=30,depth_ftm<=700,",
    "date_dim$year>=", YearRange[1], ",date_dim$year<=", YearRange[2],
    "&variables=", paste0(Vars, collapse = ",")
  )

  All.Tows <- jsonlite::fromJSON(UrlText)

  # Remove water hauls
  fix <- is.na(All.Tows[, "operation_dim$legacy_performance_code"])
  if (sum(fix) > 0) {
    All.Tows[fix, "operation_dim$legacy_performance_code"] <- -999
  }
  keep <- All.Tows[, "operation_dim$legacy_performance_code"] != 8
  All.Tows <- All.Tows[keep, ]
  All.Tows <- All.Tows[, Vars.short]

  All.Tows <- dplyr::rename(All.Tows,
    Project = project, Trawl_id = trawl_id, Year = year,
    Pass = pass, Vessel = vessel, Tow = tow, Date = datetime_utc_iso,
    Depth_m = depth_m, Longitude_dd = longitude_dd, Latitude_dd = latitude_dd,
    Area_Swept_ha = area_swept_ha_der
  )

  All.Tows <- All.Tows[
    !duplicated(paste(All.Tows$Year, All.Tows$Pass, All.Tows$Vessel, All.Tows$Tow)),
    c("Project", "Trawl_id", "Year", "Pass", "Vessel", "Tow", "Date", "Depth_m", "Longitude_dd", "Latitude_dd", "Area_Swept_ha", "station_code")
  ]

  # Link each data set together based on trawl_id
  if ("Common_name" %in% names(Data)) {
    grid <- expand.grid(
      "Trawl_id" = unique(All.Tows$Trawl_id), "Common_name" = unique(Data$Common_name),
      stringsAsFactors = FALSE
    )
  } else {
    grid <- expand.grid(
      "Trawl_id" = unique(All.Tows$Trawl_id), "Scientific_name" = unique(Data$Scientific_name),
      stringsAsFactors = FALSE
    )
  }
  Out <- dplyr::left_join(grid, All.Tows)
  Out <- dplyr::left_join(Out, Data)
  # Out = dplyr::left_join(All.Tows, Data)

  # Fill in zeros where needed
  Out$total_catch_wt_kg[is.na(Out$total_catch_wt_kg)] <- 0

  Out$CPUE_kg_per_ha[is.na(Out$CPUE_kg_per_ha)] <- 0

  Out$Subsample_count[is.na(Out$Subsample_count)] <- 0

  Out$Subsample_wt_kg[is.na(Out$Subsample_wt_kg)] <- 0

  Out$total_catch_numbers[is.na(Out$total_catch_numbers)] <- 0

  # Need to check what this is doing
  noArea <- which(is.na(Out$Area_Swept_ha))
  if (length(noArea) > 0) {
    if (verbose) {
      print(cat("\nThere are", length(noArea), "records with no area swept calculation. These record will be filled with the mean swept area across all tows.\n"))
      print(Out[noArea, c("Trawl_id", "Year", "Area_Swept_ha", "CPUE_kg_per_ha", "total_catch_numbers")])
    }
    Out[noArea, "Area_Swept_ha"] <- mean(Out$Area_Swept_ha, trim = 0.05, na.rm = TRUE)
  }

  # Scientific Name is missing after the matching when Total_sp_wt_kg is zero
  # if (!is.null(Name)) {
  #   Out$Common_name <- Species
  # }
  # if (!is.null(SciName)) {
  #   Out$Scientific_name <- Species
  # }

  Out$Date <- chron::chron(format(as.POSIXlt(Out$Date, format = "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d"), format = "y-m-d", out.format = "YYYY-m-d")

  Out$Project <- projectShort

  Out$Trawl_id <- as.character(Out$Trawl_id)

  # Convert the CPUE into km2
  Out$cpue_kg_km2 <- Out$CPUE_kg_per_ha * 100
  remove <- "CPUE_kg_per_ha"
  Out <- Out[, !(names(Out) %in% remove)]

  if (SaveFile) {
    time <- Sys.time()
    time <- substring(time, 1, 10)
    # save(Out, file = paste0(Dir, "/Catch_", outName, "_", SurveyName, "_",  time, ".rda"))
    save(Out, file = file.path(Dir, paste("Catch_", outName, "_", SurveyName, "_", time, ".rda", sep = "")))
    if (verbose) {
      message(paste("Catch data file saved to following location:", Dir))
    }
  }

  return(Out)
}
