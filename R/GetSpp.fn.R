#' Get Species Information
#'
#' Get the scientific name, common name, and strata group for
#' a vector of species.
#'
#' @param species A vector of species names. The words can be separated
#' using spaces or underscores. Full names are not required but they will
#' increase the consistency of the results should partial matches return
#' multiple matches.
#' @param unident A logical entry with the default value of \code{FALSE},
#' to match historical output that did not include unidentified groups.
#'
#' @export
#' @return A data frame with the
#' scientific name in the latin column and in the scientific_name column
#' that has spaces replaced with underscores;
#' common name in the common column and in the common_name column
#' that has spaces replaced with underscores;
#' species values used as input in the input column; and
#' strata used to assess its status in the strata column.
#' @seealso See
#' \code{\link{PullSpp.fn}} for information on common and scientific names;
#' and \code{\link{GetStrata.fn}} for the different stratifications.
#' @author Kelli Faye Johnson
#' @examples
#' GetSpp.fn(c("sablefish", "petrale"))
#' # Expect a warning because vermilion doesn't have an strata assigned to it
#' testthat::expect_warning(
#'   GetSpp.fn(c("vermilion"))
#' )
#' # Dusky returns multiple entries
#' testthat::expect_warning(
#'   GetSpp.fn(c("dusky"))
#' )
#'
GetSpp.fn <- function(species, unident = FALSE) {

  # background information
  sppnames <- PullSpp.fn()
  if (!unident) {
    sppnames <- sppnames[!grepl("unident", sppnames[["common"]]), ]
  }

  spplist <- t(data.frame(
    c("north_south", "shortbelly_rockfish"),
    c("wa_or", "pacific_cod"),
    c("ca", "brown_rockfish"),
    c("ca", "copper_rockfish"),
    c("coast", "dover_sole"),
    c("sablefish", "sablefish"),
    c("coast", "petrale_sole"),
    c("ca", "bank_rockfish"),
    c("coast", "pacific_sanddab"),
    c("coast", "arrowtooth_flounder"),
    c("coast", "aurora_rockfish"),
    c("ca", "blackgill_rockfish"),
    c("north_south", "bocaccio"),
    c("coast", "canary_rockfish"),
    c("north_south", "chilipepper"),
    c("coast", "darkblotched_rockfish"),
    c("coast", "english_sole"),
    c("wa_or", "flathead_sole"),
    c("coast", "greenspotted_rockfish"),
    c("coast", "greenstriped_rockfish"),
    c("coast", "lingcod"),
    c("coast", "longspine_thornyhead"),
    c("coast", "pacific_spiny_dogfish"),
    c("coast", "redbanded_rockfish"),
    c("coast", "rex_sole"),
    c("coast", "rougheye_rockfish"),
    c("coast", "sharpchin_rockfish"),
    c("coast", "shortspine_thornyhead"),
    c("north_south", "splitnose_rockfish"),
    c("north_south", "yellowtail_rockfish"),
    c("coast", "big_skate"),
    c("coast", "longnose_skate"),
    c("north_south", "pacific_ocean_perch"),
    c("coast", "widow_rockfish"),
    c("ca", "copper_rockfish"),
    c("ca", "squarespot_rockfish"),
    c("coast", "quillback_rockfish")
  ))
  row.names(spplist) <- NULL

  # Match species name
  index <- lapply(species, function(y) {
    out <- which(apply(apply(sppnames, 1, grepl, pattern = y, ignore.case = TRUE), 2, any))
    names(out) <- NULL
    return(out)
  })
  if (any(lapply(index, length) == 0)) {
    bad <- which(lapply(index, length) == 0)
    warning("The following species were not found in the look up table\n",
      "generated from the data warehouse (webapps.nwfsc.noaa.gov)\n",
      "and only information for species in the table will be returned:\n",
      paste(species[bad], collapse = ", "),
      call. = FALSE
    )
    species <- species[-1 * bad]
    index <- index[-1 * bad]
  }
  if (length(species) == 0) {
    stop("No matches were found for your input species.")
  }
  out <- sppnames[unlist(index), ]
  out[, "input"] <- rep(species, times = vapply(index, length, FUN.VALUE = 1L))

  # Match strata
  index <- match(tolower(out[, "common_name"]), tolower(spplist[, 2]))
  if (any(is.na(index))) {
    bad <- which(is.na(index))
    warning("The following species were not found in the look up table\n",
      "stored in GetSpp for default strata, please self-assign strata:\n",
      paste(unique(out[bad, "input"]), collapse = ", "),
      call. = FALSE
    )
  }

  out[, "strata"] <- ifelse(
    test = is.na(index),
    yes = "coast",
    no = spplist[index, 1]
    )

  return(out)
}
