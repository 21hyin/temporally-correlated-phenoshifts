# Run with R version 4.2.0
# This script has labeled code sections for easier navigation
# All figures and tables are commented as well, just search for one desired
#  (ex: "Figure 6")

library(dplyr) # version 1.1.4
# library(plyr)
library(forcats)# version 1.0.0
library(reshape2) # version 1.4.4
library(robustbase) # version 0.99-2
library(mgcv) # version 1.9-1
library(ggplot2) # version 3.5.1

# Specify directories (change as needed)
DATA_IN <- "datasets"
DATA_OUT <- "out"
DATA_OUT_RDATA <- "out/saves_RData"

# Quick start using .RData file to reproduce results ---------------------------
# after loading, it's possible to run only portions of the script
# ctrl+f or look through labeled code sections to reproduce specific results
load(file.path(DATA_OUT_RDATA, "real_community_phenology.RData"))

# Dataset helper fcns & constants ----------------------------------------------
# specifies number of years over which to calculate phenophase timings in the
#  first and last years of datasets
FIRSTLAST_NYEARS <- 10
# minimum number of years over which data exists for a species' phenophase
#  within the first and last FIRSTLAST_NYEARS of datasets
MIN_YEARS_TRACKED <- 7
# minimum number of species per community
MIN_NUM_SP <- 6
# minimum consecutive days phenophase centers must span in a community
MIN_CENTER_RANGE <- 45

# mean_n() returns the mean of values in vector x if x has >= n non-null values
#  and NA otherwise
mean_n <- function(x, n) {
    if (sum(!is.na(x)) >= n) {
        return (mean(x, na.rm = TRUE))
    }
    return (NA)
}

# firstlast_nyears_decadeshift() returns a data frame of phenological shifts
#  using species documented for at least first_n of the first n years and
#  last_n of the last n years, in the years firstyear though lastyear in the
#  data frame grouped_data
# notes on function parameters
#  year_name: column name for the year in grouped_data
#  doy_name: column name for the day of year of the phenophase of interest in
#   grouped_data
#  firstyear: year to consider as the first year of the data frame grouped_data
#  lastyear: year to consider as the last year of the data frame grouped_data
# returned data frame columns include
#  column(s) used to group grouped_data
#  firstnyears_meandoy: average day of year of the phenophase for the first
#   first_n to n years in grouped_data (depending on data availability), used to
#   compute rate of phenological shift
#  lastnyears_meandoy: average day of year of the phenophase for the last
#   last_n to n years in grouped_data (depending on data availability), used to
#   compute rate of phenological shift
#  first_year: value of firstyear
#  last_year: value of lastyear
#  [doy_name]_days_shift_per_decade: standardized decadal rate of phenological
#   shift
firstlast_nyears_decadeshift <- function(
    grouped_data,
    year_name,
    doy_name,
    firstyear,
    lastyear,
    n,
    first_n,
    last_n
) {
    first <- na.omit(summarise(
        filter(
            grouped_data,
            .data[[year_name]] >= firstyear & (
                .data[[year_name]] <= firstyear + (n - 1)
            )
        ),
        firstnyears_meandoy = mean_n(.data[[doy_name]], n = first_n)
    ))
    last <- na.omit(summarise(
        filter(
            grouped_data,
            (
                .data[[year_name]] >= lastyear - (n - 1)
            ) & .data[[year_name]] <= lastyear
        ),
        lastnyears_meandoy = mean_n(.data[[doy_name]], n = last_n)
    ))
    firstlast <- merge(first, last)
    if (nrow(firstlast) > 0) {
        firstlast$first_year <- firstyear
        firstlast$last_year <- lastyear
        firstlast[, paste0(doy_name, "_days_shift_per_decade")] <- (
            firstlast$lastn - firstlast$firstn
        ) / ((lastyear - firstyear + 1) / 10)
    }
    return (firstlast)
}

# UK Butterfly flight periods --------------------------------------------------
BFLYWINGS <- read.csv(
    file.path(DATA_IN, "ukbmsphenology2021.csv"),
    header = TRUE
)
str(BFLYWINGS)
anyNA(BFLYWINGS)

names(BFLYWINGS) <- tolower(names(BFLYWINGS))
BFLYWINGS$flightperiod_range <- BFLYWINGS$flightperiod_range + 1
BFLYWINGS$siteno <- as.factor(BFLYWINGS$siteno)
BFLYWINGS <- BFLYWINGS[order(BFLYWINGS$year), ]

# verify relationship between siteno and sitename
#  to see how to define communities
length(unique(BFLYWINGS$siteno)) > length(unique(BFLYWINGS$sitename))
# use siteno since more siteno than sitename

# select data with species and communities spanning at least 30 years,
#  considering all flight period data over a given year per species (brood = 0)
#  for our annual cycle, works best anyway since this specification makes up
#  most of the dataset and covers the most species
bflywings <- filter(
    group_by(subset(BFLYWINGS, subset = (brood == 0)), species_name, siteno),
    max(year, na.rm = TRUE) - min(year, na.rm = TRUE) + 1 >= 30
)
bflywings$siteno <- droplevels(bflywings$siteno)
bflywings_sp_siteno <- distinct(
    bflywings[, c("species_name", "siteno")]
)
# lots of species per community
hist(count(ungroup(bflywings), siteno)$n)

# quickly visualize the number of years species get tracked within the first and
#  last 10 years needed for analysis
first_year_bflywings <- min(bflywings$year, na.rm = TRUE)
hist(count(subset(
    bflywings,
    year >= first_year_bflywings & (year <= first_year_bflywings + 9)
), species_name, siteno)$n, breaks = 1:10)
hist(count(subset(
    bflywings,
    year >= first_year_bflywings + 10 & (year <= first_year_bflywings + 19)
), species_name, siteno)$n, breaks = 1:10)
last_year_bflywings <- max(bflywings$year, na.rm = TRUE)
hist(count(subset(
    bflywings,
    (year >= last_year_bflywings - 9) & (year <= last_year_bflywings)
), species_name, siteno)$n, breaks = 1:10)

# get mean flight date and duration across years
# check if need all these columns
bflywings_by_site <- summarise(
    group_by(bflywings, species_name, siteno),
    mean_mean_flight_date = mean(mean_flight_date, na.rm = TRUE),
    mean_duration = mean(flightperiod_range, na.rm = TRUE),
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    num_years_observed = length(unique(year))
)

# first and last years differ across communities
bflywings_firstlastyears <- summarise(
    group_by(bflywings_by_site, siteno),
    first_year = min(first_year, na.rm = TRUE),
    last_year = max(last_year, na.rm = TRUE)
)

# compute standardized phenological shifts for species with at least 7 years of
#  data within the first and separately for the last 10 years
bflywings_firstlasts <- by(
    bflywings,
    bflywings$siteno,
    function(x) {
        bflywings_add_siteno <- firstlast_nyears_decadeshift(
            grouped_data = group_by(x, species_name),
            year_name = "year",
            doy_name = "mean_flight_date",
            firstyear = min(x$year, na.rm = TRUE),
            lastyear = max(x$year, na.rm = TRUE),
            n = FIRSTLAST_NYEARS,
            first_n = MIN_YEARS_TRACKED,
            last_n = MIN_YEARS_TRACKED
        )
        if (nrow(bflywings_add_siteno) > 0) {
            bflywings_add_siteno$siteno <- unique(x$siteno)
        }
        return(bflywings_add_siteno)
    }
)
bflywings_firstlasts <- bind_rows(Filter(Negate(is.null), bflywings_firstlasts))
bflywings_by_site <- merge(subset(
    bflywings_by_site, select = -c(first_year, last_year, num_years_observed)
), bflywings_firstlasts)

# only keep data meeting requirements for analysis of either correlation with
#  phenophase centers or with phenophase durations
bflywings_by_site <- mutate(group_by(bflywings_by_site, siteno), num_sp = n())
bflywings_by_site <- bflywings_by_site[bflywings_by_site$num_sp >= MIN_NUM_SP, ]
bflywings_by_site <- mutate(
    group_by(bflywings_by_site, siteno),
    mean_mean_flight_date_range = max(
        mean_mean_flight_date,
        na.rm = TRUE
    ) - min(mean_mean_flight_date, na.rm = TRUE),
    mean_duration_range = max(mean_duration, na.rm = TRUE) - min(
        mean_duration,
        na.rm = TRUE
    )
)
bflywings_by_site <- filter(
    group_by(bflywings_by_site, siteno),
    (mean_mean_flight_date_range >= MIN_CENTER_RANGE & (
        median(mean_mean_flight_date, na.rm = TRUE) > min(
            mean_mean_flight_date,
            na.rm = TRUE
        ) + unique(
            mean_mean_flight_date_range
        ) / 4 & median(mean_mean_flight_date, na.rm = TRUE) < max(
            mean_mean_flight_date,
            na.rm = TRUE
        ) - unique(mean_mean_flight_date_range) / 4
    )) | (median(mean_duration, na.rm = TRUE) > min(
        mean_duration,
        na.rm = TRUE
    ) + unique(
            mean_duration_range
        ) / 4 & median(mean_duration, na.rm = TRUE) < max(
            mean_duration,
            na.rm = TRUE
        ) - unique(mean_duration_range) / 4
    )
)

# Subalpine flowering plants ---------------------------------------------------
# CSV renamed to "rmbl_phensyn.csv" from "phensyn-data-2022.csv" to better
#  distinguish from other datasets
RMBL <- read.csv(file.path(DATA_IN, "rmbl_phensyn.csv"), header = TRUE)
str(RMBL)
anyNA(RMBL)

# common species locations
rmbl_sp_loc <- list(
    Speyeria_mormonia = "Gunnison Basin",
    diptera = "Gothic Research Meadow",
    hymenoptera = "Gothic Research Meadow",
    Callospermophilus_lateralis = "Gunnison Basin",
    Ambystoma_mavortium_nebulosum = "Mexican Cut Nature Preserve"
)

# standardize names of communities based on metadata and personal communication
#  with Brian Inouye for analysis
RMBL$site_name <- gsub("[0-9.]", "", RMBL$plot)
RMBL$habitat <- ifelse(RMBL$site_name == "VR", "EM", RMBL$site_name)
RMBL$habitat <- ifelse(
    RMBL$plot == "RM8" | RMBL$plot == "RM9",
    # wet, underneath aspen forest canopy
    "AS",
    RMBL$habitat
)

# Select data with species and communities spanning at least 30 years
rmbl <- filter(
    group_by(RMBL, species, event, habitat),
    max(year, na.rm = TRUE) - min(year, na.rm = TRUE) + 1 >= 30
)
anyNA(rmbl)
rmbl <- rmbl[
    !is.na(rmbl$habitat) & !is.na(rmbl$species) & !is.na(rmbl$event) & !is.na(
      rmbl$plot) & !is.na(rmbl$year),
]
# only the flowering phenophase tracked at community level
length(unique(rmbl$event))

# can assess correlation with phenophase timing and duration by using all
#  filtered species in communities, no need to create separate data frames to
#  assess each of the two correlations
any(c(is.na(rmbl$first.doy), is.na(rmbl$last.doy), is.na(rmbl$median.doy)))

# quickly visualize the number of years species get tracked within the first and
#  last 10 years needed for analysis, cross-reference for next step
rmbl_sp_year_habitat <- distinct(rmbl[, c("species", "habitat", "year")])
first_year_rmbl <- min(rmbl$year, na.rm = TRUE)
hist(count(subset(
    rmbl_sp_year_habitat,
    year >= first_year_rmbl & (year <= first_year_rmbl + 9)
), species, habitat)$n, breaks = 1:10)
last_year_rmbl <- max(rmbl$year, na.rm = TRUE)
hist(count(subset(
    rmbl_sp_year_habitat,
    (year >= last_year_rmbl - 9) & (year <= last_year_rmbl)
), species, habitat)$n, breaks = 1:10)

# Pick single day of phenophase across plots within a habitat because we treat
#  plots under the same habitat as one community
rmbl <- summarise(
    group_by(rmbl, species, event, habitat, year),
    median.doy = median(median.doy, na.rm = TRUE),
    first.doy = median(first.doy, na.rm = TRUE),
    last.doy = median(last.doy, na.rm = TRUE)
)

# get mean of median flowering day and flowering duration across years
rmbl_by_habitat <- summarise(
    group_by(rmbl, species, event, habitat),
    mean_mediandoy = mean(median.doy, na.rm = TRUE),
    mean_duration = mean(last.doy - first.doy + 1, na.rm = TRUE),
    first_year = min(year, na.rm = TRUE),
    last_year = max(year, na.rm = TRUE),
    num_years_observed = length(unique(year))
)

# Compute standardized phenological shifts for species with at least 7 years of
#  data within the first and separately for the last 10 years

# all the same years, easy to deal with
rmbl_firstlastyears <- summarise(
    group_by(rmbl_by_habitat, event, habitat),
    first_year = min(first_year, na.rm = TRUE),
    last_year = max(last_year, na.rm = TRUE)
)

rmbl_by_habitat <- merge(subset(
    rmbl_by_habitat, select = -c(first_year, last_year, num_years_observed)
), firstlast_nyears_decadeshift(
    grouped_data = group_by(rmbl, species, habitat, event),
    year_name = "year",
    doy_name = "median.doy",
    firstyear = unique(rmbl_firstlastyears$first_year),
    lastyear = unique(rmbl_firstlastyears$last_year),
    n = FIRSTLAST_NYEARS,
    first_n = MIN_YEARS_TRACKED,
    last_n = MIN_YEARS_TRACKED
))

# only keep data meeting requirements for analysis of either correlation with
#  phenophase centers or with phenophase durations
rmbl_by_habitat <- mutate(group_by(rmbl_by_habitat, habitat), num_sp = n())
rmbl_by_habitat <- rmbl_by_habitat[rmbl_by_habitat$num_sp >= MIN_NUM_SP, ]
rmbl_by_habitat <- mutate(
    group_by(rmbl_by_habitat, habitat),
    mean_mediandoy_range = max(mean_mediandoy, na.rm = TRUE) - min(
        mean_mediandoy,
        na.rm = TRUE
    ),
    mean_duration_range = max(mean_duration, na.rm = TRUE) - min(
        mean_duration,
        na.rm = TRUE
    )
)
rmbl_by_habitat <- filter(
    group_by(rmbl_by_habitat, habitat),
    (mean_mediandoy_range >= MIN_CENTER_RANGE & (
        median(mean_mediandoy, na.rm = TRUE) > min(
            mean_mediandoy,
            na.rm = TRUE
        ) + unique(mean_mediandoy_range) / 4 & median(
            mean_mediandoy,
            na.rm = TRUE
        ) < max(mean_mediandoy, na.rm = TRUE) - unique(mean_mediandoy_range) / 4
    )) | (median(mean_duration, na.rm = TRUE) > min(
        mean_duration,
        na.rm = TRUE
    ) + unique(mean_duration_range) / 4 & median(
        mean_duration,
        na.rm = TRUE
    ) < max(mean_duration, na.rm = TRUE) - unique(mean_duration_range) / 4)
)

# Boreal forest plant phenophases ----------------------------------------------
# Rdata file renamed to barguzin.Rdata from "PhenoDataCleaned_2.0.Rdata" to
#  better distinguish from other datasets
load(file.path(DATA_IN, "barguzin.Rdata"))
BARGUZIN <- pheno
str(BARGUZIN)
anyNA(BARGUZIN)
unique(BARGUZIN$Phenophase_English)

# remove duplicate rows (likely errors)
BARGUZIN <- BARGUZIN[which(!duplicated(BARGUZIN)), ]

# remove stages with data on "onset" without a corresponding "end", since we
#  cannot compute durations from them
BARGUZIN <- subset(
    BARGUZIN,
    Phenophase_English != "onset of flower budding" &
        Phenophase_English != "seed set visible" &
        Phenophase_English != "onset of fruiting" &
        Phenophase_English != "onset of fruit dispersion" &
        Phenophase_English != "onset of leaf unfolding" &
        Phenophase_English != "onset of autumn colouring" &
        Phenophase_English != "full autumn colouring of leaves"
)

barguzin <- filter(
    group_by(BARGUZIN, Species, Phenophase_English),
    max(Year, na.rm = TRUE) - min(Year, na.rm = TRUE) + 1 >= 30
)
barguzin <- barguzin[
    !is.na(barguzin$Site) & !is.na(barguzin$Species) & !is.na(
        barguzin$Year
    ) & !is.na(barguzin$Phenophase_English),
]
stopifnot(!anyNA(barguzin))

# Consider only one instance of event for each year-site-sp but note for some
#  species this is not always true (multi-modal distributions)

phenophase_cnt <- count(barguzin, Year, Site, Species, Phenophase_English)
phenophase_dups <- phenophase_cnt[phenophase_cnt$n > 1, ]
# ex: for Antennaria dioica (L.) Gaertn., mass blooming may occur twice a year
unique(phenophase_cnt$Phenophase_English)

# get doys for sp-year across sites while selecting earliest dates for
#  consistency
barguzin <- subset(mutate(
    group_by(barguzin, Species, Year, Phenophase_English, Site),
    Day = min(Day, na.rm = TRUE)
), select = -Date)

# Pick single day of phenophase across sites because we treat sites as one
#  community
barguzin <- summarise(
    group_by(barguzin, Species, Year, Phenophase_English),
    Day = median(Day, na.rm = TRUE)
)

# Standardize names of phenophases for analysis

barguzin$Phenophase_English <- fct_recode(
    barguzin$Phenophase_English,
    "end of leaf fall" = "leaf fall end"
)

barguzin$event <- ""
barguzin$event[c(
    grep("onset", barguzin$Phenophase_English),
    grep("withering away", barguzin$Phenophase_English)
)] <- vapply(
    X = barguzin$Phenophase_English[c(
        grep("onset", barguzin$Phenophase_English),
        grep("withering away", barguzin$Phenophase_English)
    )],
    FUN = function(x) {
        if (x == "withering away") {
            return ("vegetative growth in herbaceous plants")
        }
        return (sub("onset of (\\w*)", "\\1", x))
    },
    FUN.VALUE = "character",
    USE.NAMES = FALSE
)
barguzin$event[c(grep("end", barguzin$Phenophase_English))] <- vapply(
    X = barguzin$Phenophase_English[
        c(grep("end", barguzin$Phenophase_English))
    ],
    FUN = function(x) { return (sub("end of (\\w*)", "\\1", x)) },
    FUN.VALUE = "character",
    USE.NAMES = FALSE
)
barguzin$event[c(grep("mass", barguzin$Phenophase_English))] <- vapply(
    X = barguzin$Phenophase_English[
        c(grep("mass", barguzin$Phenophase_English))
    ],
    FUN = function(x) { return (sub("mass (\\w*)", "\\1", x)) },
    FUN.VALUE = "character",
    USE.NAMES = FALSE
)
barguzin$event <- as.factor(barguzin$event)

# assign first & last for duration, plus center
barguzin$event_type <- ""
barguzin$event_type[c(grep("onset", barguzin$Phenophase_English))] <- "firstday"
barguzin$event_type[c(
    grep("end", barguzin$Phenophase_English),
    grep("full", barguzin$Phenophase_English),
    grep("withering away", barguzin$Phenophase_English)
)] <- "lastday"
barguzin$event_type[c(grep("mass", barguzin$Phenophase_English))] <- "center"
barguzin$event_type <- as.factor(barguzin$event_type)

# Prepare data frames for data that can only work for analyses on correlation
#  with phenophase centers versus with phenophase durations

barguzin_cast <- dcast(
    # subset(barguzin, select = -c(Phenophase_English, Phenophase_Abb)),
    subset(barguzin, select = -Phenophase_English),
    ... ~ event_type,
    fun.aggregate = NULL,
    value.var = "Day"
)
stopifnot(!any(rowSums(is.na(
    barguzin_cast[, c("firstday", "lastday", "center")]
)) == 3))

barguzin_cast$duration <- barguzin_cast$lastday - barguzin_cast$firstday + 1
# assumption for phenophases without documented central dates, these phenophases
#  end up getting filtered out anyway
barguzin_cast[is.na(barguzin_cast$center), ]$center <- (
    barguzin_cast[is.na(barguzin_cast$center), ]$firstday + barguzin_cast[is.na(
        barguzin_cast$center
    ), ]$lastday
) / 2

# Create data frame for phenophase centers and associated phenological shifts
#  while treating species of the same phenophase as a community

# when first and/or last day is NA even when attempting to compute a midpoint
barguzin_cast_centers <- barguzin_cast[!is.na(barguzin_cast$center), ]

# quickly visualize the number of years species get tracked within the first and
#  last 10 years needed for analysis
first_year_barguzin_center <- min(barguzin_cast_centers$Year, na.rm = TRUE)
hist(count(subset(
    barguzin_cast_centers,
    Year >= first_year_barguzin_center & (
        Year <= first_year_barguzin_center + 9
    )
), Species, event)$n, breaks = 1:10)
hist(count(subset(
    barguzin_cast_centers,
    Year >= first_year_barguzin_center + 10 & (
        Year <= first_year_barguzin_center + 19
    )
), Species, event)$n, breaks = 1:10)
hist(count(subset(
    barguzin_cast_centers,
    Year >= first_year_barguzin_center + 20 & (
        Year <= first_year_barguzin_center + 29
    )
), Species, event)$n, breaks = 1:10)
last_year_barguzin_center <- max(barguzin_cast_centers$Year, na.rm = TRUE)
hist(count(subset(
    barguzin_cast_centers,
    (Year >= last_year_barguzin_center - 9) & (
        Year <= last_year_barguzin_center
    )
), Species, event)$n, breaks = 1:10)

# start with data frame of mean phenophase centers across years
# first_year_barguzin_center + 10 to try to capture more communities
barguzin_center_summary <- summarise(
    group_by(barguzin_cast_centers[
        barguzin_cast_centers$Year >= first_year_barguzin_center + 10,
    ], Species, event),
    mean_center = mean(center, na.rm = TRUE),
    first_year = min(Year, na.rm = TRUE),
    last_year = max(Year, na.rm = TRUE),
    num_years_observed = length(unique(Year))
)

# last years differ across communities
barguzin_center_firstlastyears <- summarise(
    group_by(barguzin_center_summary, event),
    first_year = min(first_year, na.rm = TRUE),
    last_year = max(last_year, na.rm = TRUE)
)

# compute standardized phenological shifts for species with at least 7 years of
#  data within the first and separately for the last 10 years
barguzin_center_firstlasts <- by(
    barguzin_cast_centers[
        barguzin_cast_centers$Year >= first_year_barguzin_center + 10,
    ],
    barguzin_cast_centers[
        barguzin_cast_centers$Year >= first_year_barguzin_center + 10,
    ]$event,
    function(x) {
        barguzin_center_add_event <- firstlast_nyears_decadeshift(
            grouped_data = group_by(x, Species),
            year_name = "Year",
            doy_name = "center",
            firstyear = min(x$Year, na.rm = TRUE),
            lastyear = max(x$Year, na.rm = TRUE),
            n = FIRSTLAST_NYEARS,
            first_n = MIN_YEARS_TRACKED,
            last_n = MIN_YEARS_TRACKED
        )
        if (nrow(barguzin_center_add_event) > 0) {
            barguzin_center_add_event$event <- unique(x$event)
        }
        return(barguzin_center_add_event)
    }
)
barguzin_center_firstlasts <- bind_rows(Filter(
    Negate(is.null),
    barguzin_center_firstlasts)
)
barguzin_center <- merge(subset(
    barguzin_center_summary, select = -c(
        first_year,
        last_year,
        num_years_observed
    )
), barguzin_center_firstlasts)
unique(barguzin_center$event)

# only keep data meeting requirements for analysis
barguzin_center <- mutate(group_by(barguzin_center, event), num_sp = n())
barguzin_center <- barguzin_center[barguzin_center$num_sp >= MIN_NUM_SP, ]
barguzin_center <- mutate(
    group_by(barguzin_center, event),
    mean_center_range = max(mean_center, na.rm = TRUE) - min(
        mean_center,
        na.rm = TRUE
    )
)
barguzin_center <- barguzin_center[
    barguzin_center$mean_center_range >= MIN_CENTER_RANGE,
]
stopifnot(
    median(barguzin_center$mean_center) > min(
        barguzin_center$mean_center,
        na.rm = TRUE
    ) + unique(barguzin_center$mean_center_range) / 4 & median(
        barguzin_center$mean_center,
        na.rm = TRUE
    ) < max(barguzin_center$mean_center, na.rm = TRUE) - unique(
        barguzin_center$mean_center_range
    ) / 4
)
unique(barguzin_center$event)
barguzin_center$event <- droplevels(barguzin_center$event)

# Create data frame for phenophase durations and associated phenological shifts
#  while treating species of the same phenophase as a community

barguzin_cast_durations <- barguzin_cast[!is.na(barguzin_cast$duration), ]

# quickly visualize the number of years species get tracked within the first and
#  last 10 years needed for analysis
first_year_barguzin_duration <- min(barguzin_cast_durations$Year, na.rm = TRUE)
hist(count(subset(
    barguzin_cast_durations,
    Year >= first_year_barguzin_duration & (
        Year <= first_year_barguzin_duration + 9
    )
), Species, event)$n, breaks = 1:10)
hist(count(subset(
    barguzin_cast_durations,
    Year >= first_year_barguzin_duration + 10 & (
        Year <= first_year_barguzin_duration + 19
    )
), Species, event)$n, breaks = 1:10)
hist(count(subset(
    barguzin_cast_durations,
    Year >= first_year_barguzin_duration + 20 & (
        Year <= first_year_barguzin_duration + 29
    )
), Species, event)$n, breaks = 1:10)
last_year_barguzin_duration <- max(barguzin_cast_durations$Year, na.rm = TRUE)
hist(count(subset(
    barguzin_cast_durations,
    Year >= last_year_barguzin_duration - 9 & (
        Year <= last_year_barguzin_duration
    )
), Species, event)$n, breaks = 1:10)

# start with data frame of mean phenophase durations across years
# first_year_barguzin_duration + 10 to try to capture more communities
barguzin_duration_summary <- summarise(
    group_by(barguzin_cast_durations[
        barguzin_cast_durations$Year >= first_year_barguzin_duration + 10,
    ], Species, event),
    mean_duration = mean(duration, na.rm = TRUE),
    first_year = min(Year, na.rm = TRUE),
    last_year = max(Year, na.rm = TRUE),
    num_years_observed = length(unique(Year))
)

# last years differ across communities
barguzin_duration_firstlastyears <- summarise(
    group_by(barguzin_duration_summary, event),
    first_year = min(first_year, na.rm = TRUE),
    last_year = max(last_year, na.rm = TRUE)
)

barguzin_duration_firstlasts <- by(
    barguzin_cast_durations[
        barguzin_cast_durations$Year >= first_year_barguzin_duration + 10,
    ],
    barguzin_cast_durations[
        barguzin_cast_durations$Year >= first_year_barguzin_duration + 10,
    ]$event,
    function(x) {
        barguzin_duration_add_event <- firstlast_nyears_decadeshift(
            grouped_data = group_by(x, Species),
            year_name = "Year",
            doy_name = "duration",
            firstyear = min(x$Year, na.rm = TRUE),
            lastyear = max(x$Year, na.rm = TRUE),
            n = FIRSTLAST_NYEARS,
            first_n = MIN_YEARS_TRACKED,
            last_n = MIN_YEARS_TRACKED
        )
        if (nrow(barguzin_duration_add_event) > 0) {
            barguzin_duration_add_event$event <- unique(x$event)
        }
        return(barguzin_duration_add_event)
    }
)
barguzin_duration_firstlasts <- bind_rows(Filter(
    Negate(is.null),
    barguzin_duration_firstlasts
))
barguzin_duration <- merge(subset(
    barguzin_duration_summary, select = -c(
        first_year,
        last_year,
        num_years_observed
    )
), barguzin_duration_firstlasts)
unique(barguzin_duration$event)

# only keep data meeting requirements for analysis
barguzin_duration <- mutate(group_by(barguzin_duration, event), num_sp = n())
barguzin_duration <- barguzin_duration[barguzin_duration$num_sp >= MIN_NUM_SP, ]
barguzin_duration <- mutate(
    group_by(barguzin_duration, event),
    mean_duration_range = max(mean_duration, na.rm = TRUE) - min(
        mean_duration,
        na.rm = TRUE
    )
)
barguzin_duration <- filter(
    group_by(barguzin_duration, event),
    median(mean_duration, na.rm = TRUE) > min(
        mean_duration,
        na.rm = TRUE
    ) + unique(mean_duration_range) / 4 & median(
        mean_duration,
        na.rm = TRUE
    ) < max(mean_duration, na.rm = TRUE) - unique(mean_duration_range) / 4
)
unique(barguzin_duration$event)
barguzin_duration$event <- droplevels(barguzin_duration$event)

# Consolidate datasets ---------------------------------------------------------
bflywings_by_site$community_id <- sapply(
    bflywings_by_site$siteno, \(x) paste0("ukbms", x)
)
bflywings_by_site <- bflywings_by_site[,
    !(colnames(bflywings_by_site) %in% c("siteno"))
]
bflywings_by_site <- rename(
    bflywings_by_site,
    sp = species_name,
    center_phenophase = mean_mean_flight_date,
    duration_phenophase = mean_duration,
    days_shift_per_decade = mean_flight_date_days_shift_per_decade,
    center_phenophase_range = mean_mean_flight_date_range,
    duration_phenophase_range = mean_duration_range
)
bflywings_by_site$dataset_name <- "ukbms"
bflywings_by_site$event <- "flight"
bflywings_by_site$year_start <- min(bflywings_by_site$first_year)
bflywings_by_site$year_end <- max(bflywings_by_site$last_year)
bflywings_by_site$dataset_years <- max(bflywings_by_site$last_year) - min(
    bflywings_by_site$first_year
) + 1

rmbl_by_habitat$community_id <- sapply(
    rmbl_by_habitat$habitat,
    \(x) paste0("rmbl", x)
)
rmbl_by_habitat <- rmbl_by_habitat[,
    !(colnames(rmbl_by_habitat) %in% c("mean_peakdoy", "habitat"))
]
rmbl_by_habitat <- rename(
    rmbl_by_habitat,
    sp = species,
    center_phenophase = mean_mediandoy,
    duration_phenophase = mean_duration,
    days_shift_per_decade = median.doy_days_shift_per_decade,
    center_phenophase_range = mean_mediandoy_range,
    duration_phenophase_range = mean_duration_range
)
rmbl_by_habitat$dataset_name <- "rmbl"
rmbl_by_habitat$year_start <- min(rmbl_by_habitat$first_year)
rmbl_by_habitat$year_end <- max(rmbl_by_habitat$last_year)
rmbl_by_habitat$dataset_years <- max(rmbl_by_habitat$last_year) - min(
    rmbl_by_habitat$first_year
) + 1

barguzin_center$community_id <- sapply(
    barguzin_center$event,
    \(x) paste0("barguzin-", x)
)
barguzin_center <- rename(
    barguzin_center,
    sp = Species,
    center_phenophase = mean_center,
    days_shift_per_decade = center_days_shift_per_decade,
    center_phenophase_range = mean_center_range
)
barguzin_center$dataset_name <- "barguzin_centers"
barguzin_center$duration_phenophase <- NA
barguzin_center$duration_phenophase_range <- NA
barguzin_center$year_start <- min(barguzin_center$first_year)
barguzin_center$year_end <- max(barguzin_center$last_year)
barguzin_center$dataset_years <- max(barguzin_center$last_year) - min(
    barguzin_center$first_year
) + 1

barguzin_duration$community_id <- sapply(
    barguzin_duration$event, \(x) paste0("barguzin-", x)
)
barguzin_duration <- rename(
    barguzin_duration,
    sp = Species,
    duration_phenophase = mean_duration,
    days_shift_per_decade = duration_days_shift_per_decade,
    duration_phenophase_range = mean_duration_range
)
barguzin_duration$dataset_name <- "barguzin_durations"
barguzin_duration$center_phenophase <- NA
barguzin_duration$center_phenophase_range <- NA
barguzin_duration$year_start <- min(barguzin_duration$first_year)
barguzin_duration$year_end <- max(barguzin_duration$last_year)
barguzin_duration$dataset_years <- max(barguzin_duration$last_year) - min(
    barguzin_duration$first_year
) + 1

communities <- rbind(
    bflywings_by_site,
    rmbl_by_habitat,
    barguzin_center,
    barguzin_duration
)
communities$community_id <- as.factor(communities$community_id)

# total number of UKMBS communities noted in methods
length(unique(communities[communities$dataset_name == "ukbms", ]$community_id))

# store data for analyses on correlation with phenophase centers
center_D <- rbind(
    filter(
        group_by(
            communities[!grepl("barguzin", communities$dataset_name), ],
            community_id
        ),
        center_phenophase_range >= MIN_CENTER_RANGE & (
            median(center_phenophase) > min(center_phenophase) + unique(
                center_phenophase_range
            ) / 4 & median(center_phenophase) < max(center_phenophase) - unique(
                center_phenophase_range
            ) / 4
        )
    ),
    barguzin_center
)
center_D <- center_D[!is.na(center_D$center_phenophase), ]
center_D$community_id <- as.factor(center_D$community_id)
center_D$dataset_name <- fct_recode(
    factor(center_D$dataset_name, levels = c(
        "barguzin_centers",
        "rmbl",
        "ukbms"
    )),
    "barguzin" = "barguzin_centers"
)

# store data for analyses on correlation with phenophase durations
duration_D <- rbind(
    filter(
        group_by(
            communities[!grepl("barguzin", communities$dataset_name), ],
            community_id
        ),
        median(duration_phenophase) > min(duration_phenophase) + unique(
            duration_phenophase_range
        ) / 4 & median(duration_phenophase) < max(duration_phenophase) - unique(
            duration_phenophase_range
        ) / 4
    ),
    barguzin_duration
)
duration_D <- duration_D[!is.na(duration_D$duration_phenophase), ]
duration_D$community_id <- as.factor(duration_D$community_id)
duration_D$dataset_name <- fct_recode(
    factor(duration_D$dataset_name, levels = c(
        "barguzin_durations",
        "rmbl",
        "ukbms"
    )),
    "barguzin" = "barguzin_durations"
)

# Show data for Appendix S1: Table S2 & Appendix S1: Fig. S7 -------------------
# for part of Appendix S1: Table S2
summarise(
    group_by(center_D, dataset_name),
    first_year = min(first_year),
    last_year = max(last_year),
    num_communities = length(unique(community_id))
)
summarise(
    group_by(duration_D, dataset_name),
    first_year = min(first_year),
    last_year = max(last_year),
    num_communities = length(unique(community_id))
)

# Appendix S1: Figure S7a
count(center_D[center_D$dataset_name == "rmbl", ], community_id)
count(duration_D[duration_D$dataset_name == "rmbl", ], community_id)

# Appendix S1: Figure S7b
png(
    filename = file.path(DATA_OUT, "figure_s7b.png"),
    width = 8.5,
    height = 8.5,
    units = "cm",
    pointsize = 9,
    res = 600,
    bg = "transparent"
)
hist(
    count(communities[communities$dataset_name == "ukbms", ], community_id)$n,
    breaks = 4:28,
    main = "",
    xlim = c(5, 30),
    xlab = "# of species per community",
    ylab = "# of communities"
)
dev.off()

# Appendix S1: Figure S7c
count(center_D[center_D$dataset_name == "barguzin", ], event)
count(duration_D[duration_D$dataset_name == "barguzin", ], event)


# Analysis across datasets -----------------------------------------------------
# significance levels to test non-zero slope of linear correlation and nonlinear
#  correlation
ALPHA_LEVEL <- 0.05

# determine and characterize linear or nonlinear correlation between
#  phenological shift and phenophase center/timing
communities_center_shapeshift <- lapply(
    split(center_D, center_D$community_id, drop = TRUE),
    function(df) {
        # robust linear regression
        m_linear <- lmrob(
            days_shift_per_decade ~ center_phenophase,
            data = df,
            # recommended setting
            setting = "KS2014"
        )
        linear_resids <- residuals(m_linear)
        # generalized additive modeling
        m_nonlinear <- gam(
            days_shift_per_decade ~ s(center_phenophase, bs = "tp", k = 6),
            data = df,
            method = "REML"
        )
        # effective degrees of freedom from gam to check for linearity
        edf <- summary(m_nonlinear)$s.table[1, c("edf")]
        # general gam diagnostics
        k_check <- k.check(m_nonlinear)
        # prevent gam's sensitivity to outliers that could turn a linear
        #  relationship into a nonlinear one
        if (round(edf) != 1) {
            q1 <- quantile(linear_resids, probs = 0.25)
            q3 <- quantile(linear_resids, probs = 0.75)
            iqr <- q3 - q1
            df_test <- df[which(
                linear_resids <= q3 + (
                    1.5 * iqr
                ) & linear_resids >= q1 - (1.5 * iqr)
            ), ]
            m_nonlinear_test <- gam(
                days_shift_per_decade ~ s(center_phenophase, bs = "tp", k = 6),
                data = df_test,
                method = "REML"
            )
            edf_test <- summary(m_nonlinear_test)$s.table[1, c("edf")]
            if (round(edf_test) == 1) {
                edf <- edf_test
                print(paste(
                    df$community_id[1],
                    "has",
                    dim(df)[1] - dim(df_test)[1],
                    "outlier(s),",
                    "changed from nonlinear to linear"
                ))
            }
        }
        # Kendall's tau especially helpful for nonlinear correlation
        tau <- cor.test(
            df$days_shift_per_decade,
            df$center_phenophase,
            method = "kendall"
        )
        return (c(
            round(coef(m_linear)[2], 3),
            round(coef(m_linear)[1], 3),
            confint(m_linear)[2, ],
            summary(m_linear)$coefficients[2, 4],
            min(m_linear$fitted.values),
            max(m_linear$fitted.values),
            k_check[1, c("p-value")],
            # roughly evaluate heterosckedascity
            shapiro.test(residuals.gam(m_nonlinear, type = "response"))$p.value,
            edf,
            min(m_nonlinear$fitted.values),
            max(m_nonlinear$fitted.values),
            tau$estimate,
            tau$p.value
        ))
    }
)
communities_center_shapeshift_stacked <- stack(communities_center_shapeshift)
communities_center_shapeshift_stacked$param <- rep(
    c(
        "slope",
        "intercept",
        "slope_95ci_lower",
        "slope_95ci_upper",
        "slope_p",
        "slope_min_fitted_y",
        "slope_max_fitted_y",
        "k_check_p",
        "gamresids_sw_p",
        "edf",
        "gam_min_fitted_y",
        "gam_max_fitted_y",
        "kendall",
        "kendall_p"
    ),
    times = dim(communities_center_shapeshift_stacked)[1] / 14
)
communities_center_shapeshift <- unstack(
    communities_center_shapeshift_stacked,
    form = values ~ param
)
communities_center_shapeshift$community_id <- unique(
    communities_center_shapeshift_stacked$ind
)
communities_center_shapeshift$nonzero_slope <- ifelse(
    communities_center_shapeshift$slope_p < ALPHA_LEVEL,
    "yes",
    "no"
)
communities_center_shapeshift$nonzero_slope <- factor(
    communities_center_shapeshift$nonzero_slope,
    levels = c("yes", "no")
)
communities_center_shapeshift$correlated <- ifelse(
    communities_center_shapeshift$kendall_p < ALPHA_LEVEL,
    "yes",
    "no"
)
communities_center_shapeshift$correlated <- factor(
    communities_center_shapeshift$correlated,
    levels = c("yes", "no")
)
communities_center_shapeshift <- mutate(
    communities_center_shapeshift,
    shift_case_num = case_when(
        round(edf) != 1 & gam_min_fitted_y >= 0 & gam_max_fitted_y > 0 ~ "nonlinear: shift later",
        round(edf) != 1 & gam_min_fitted_y < 0 & gam_max_fitted_y <= 0 ~ "nonlinear: shift earlier",
        round(edf) != 1 ~ "nonlinear: some shift later, others shift earlier",
        round(edf) == 1 & slope_min_fitted_y >= 0 & slope_max_fitted_y > 0 ~ "linear: shift later",
        round(edf) == 1 & slope_min_fitted_y < 0 & slope_max_fitted_y <= 0 ~ "linear: shift earlier",
        round(edf) == 1 ~ "linear: some shift later, others shift earlier"
    )
)
communities_center_shapeshift$shift_case_num <- factor(
    communities_center_shapeshift$shift_case_num
)
communities_center_shapeshift$dataset_name <- communities_center_shapeshift$community_id
communities_center_shapeshift <- mutate(
    communities_center_shapeshift,
    dataset_name = case_when(
        grepl("ukbms", dataset_name) ~ "ukbms",
        grepl("rmbl", dataset_name) ~ "rmbl",
        grepl("barguzin", dataset_name) ~ "barguzin",
        .default = dataset_name
    ),
    dataset_abbrev = case_when(
        grepl("ukbms", dataset_name) ~ "UKBMS",
        grepl("rmbl", dataset_name) ~ "RMBL",
        grepl("barguzin", dataset_name) ~ "BNP",
        .default = dataset_name
    )
)

# determine and characterize linear or nonlinear correlation between
#  phenological shift and phenophase duration
communities_duration_shapeshift <- lapply(
    split(duration_D, duration_D$community_id, drop = TRUE),
    function(df) {
        # robust linear regression
        m_linear <- lmrob(
            days_shift_per_decade ~ duration_phenophase,
            data = df,
            # recommended setting
            setting = "KS2014"
        )
        linear_resids <- residuals(m_linear)
        # generalized additive modeling
        m_nonlinear <- gam(
            days_shift_per_decade ~ s(duration_phenophase, bs = "tp", k = 6),
            data = df,
            method = "REML"
        )
        # effective degrees of freedom from gam to check for linearity
        edf <- summary(m_nonlinear)$s.table[1, c("edf")]
        # general gam diagnostics
        k_check <- k.check(m_nonlinear)
        # prevent gam's sensitivity to outliers that could turn a linear
        #  relationship into a nonlinear one
        if (round(edf) != 1) {
            q1 <- quantile(linear_resids, probs = 0.25)
            q3 <- quantile(linear_resids, probs = 0.75)
            iqr <- q3 - q1
            df_test <- df[which(
                linear_resids <= q3 + (
                    1.5 * iqr
                ) & linear_resids >= q1 - (1.5 * iqr)
            ), ]
            m_nonlinear_test <- gam(
                days_shift_per_decade ~ s(duration_phenophase, bs = "tp", k = 6),
                data = df_test,
                method = "REML"
            )
            edf_test <- summary(m_nonlinear_test)$s.table[1, c("edf")]
            if (round(edf_test) == 1) {
                edf <- edf_test
                print(paste(
                    df$community_id[1],
                    "has",
                    dim(df)[1] - dim(df_test)[1],
                    "outlier(s),",
                    "changed from nonlinear to linear"
                ))
            }
        }
        # Kendall's tau especially helpful for nonlinear correlation
        tau <- cor.test(
            df$days_shift_per_decade,
            df$duration_phenophase,
            method = "kendall"
        )
        return (c(
            round(coef(m_linear)[2], 3),
            round(coef(m_linear)[1], 3),
            confint(m_linear)[2, ],
            summary(m_linear)$coefficients[2, 4],
            min(m_linear$fitted.values),
            max(m_linear$fitted.values),
            k_check[1, c("p-value")],
            # roughly evaluate heterosckedascity
            shapiro.test(residuals.gam(m_nonlinear, type = "response"))$p.value,
            summary(m_nonlinear)$s.table[1, c("edf")],
            min(m_nonlinear$fitted.values),
            max(m_nonlinear$fitted.values),
            tau$estimate,
            tau$p.value,
            unique(df$duration_phenophase_range)
        ))
    }
)
communities_duration_shapeshift_stacked <- stack(communities_duration_shapeshift)
communities_duration_shapeshift_stacked$param <- rep(
    c(
        "slope",
        "intercept",
        "slope_95ci_lower",
        "slope_95ci_upper",
        "slope_p",
        "slope_min_fitted_y",
        "slope_max_fitted_y",
        "k_check_p",
        "gamresids_sw_p",
        "edf",
        "gam_min_fitted_y",
        "gam_max_fitted_y",
        "kendall",
        "kendall_p",
        "duration_phenophase_range"
    ),
    times = dim(communities_duration_shapeshift_stacked)[1] / 15
)
communities_duration_shapeshift <- unstack(
    communities_duration_shapeshift_stacked,
    form = values ~ param
)
communities_duration_shapeshift$community_id <- unique(
    communities_duration_shapeshift_stacked$ind
)
communities_duration_shapeshift$nonzero_slope <- ifelse(
    communities_duration_shapeshift$slope_p < ALPHA_LEVEL,
    "yes",
    "no"
)
communities_duration_shapeshift$nonzero_slope <- factor(
    communities_duration_shapeshift$nonzero_slope,
    levels = c("yes", "no")
)
communities_duration_shapeshift$correlated <- ifelse(
    communities_duration_shapeshift$kendall_p < ALPHA_LEVEL,
    "yes",
    "no"
)
communities_duration_shapeshift$correlated <- factor(
    communities_duration_shapeshift$correlated,
    levels = c("yes", "no")
)
communities_duration_shapeshift <- mutate(
    communities_duration_shapeshift,
    shift_case_num = case_when(
        round(edf) != 1 & gam_min_fitted_y >= 0 & gam_max_fitted_y > 0 ~ "nonlinear: shift later",
        round(edf) != 1 & gam_min_fitted_y < 0 & gam_max_fitted_y <= 0 ~ "nonlinear: shift earlier",
        round(edf) != 1 ~ "nonlinear: some shift later, others shift earlier",
        round(edf) == 1 & slope_min_fitted_y >= 0 & slope_max_fitted_y > 0 ~ "linear: shift later",
        round(edf) == 1 & slope_min_fitted_y < 0 & slope_max_fitted_y <= 0 ~ "linear: shift earlier",
        round(edf) == 1 ~ "linear: some shift later, others shift earlier"
    )
)
communities_duration_shapeshift$shift_case_num <- factor(
    communities_duration_shapeshift$shift_case_num
)
communities_duration_shapeshift$dataset_name <- communities_duration_shapeshift$community_id
communities_duration_shapeshift <- mutate(
    communities_duration_shapeshift,
    dataset_name = case_when(
        grepl("ukbms", dataset_name) ~ "ukbms",
        grepl("rmbl", dataset_name) ~ "rmbl",
        grepl("barguzin", dataset_name) ~ "barguzin",
        .default = dataset_name
    ),
    dataset_abbrev = case_when(
        grepl("ukbms", dataset_name) ~ "UKBMS",
        grepl("rmbl", dataset_name) ~ "RMBL",
        grepl("barguzin", dataset_name) ~ "BNP",
        .default = dataset_name
    )
)
communities_duration_shapeshift <- merge(
    communities_duration_shapeshift,
    distinct(duration_D[, c("num_sp", "community_id")]),
    by = "community_id",
    all.x = TRUE
)

# Plot analyses from across datasets (Fig. 6, Appendix S1: Fig. S8) ------------
mytheme <- theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black")
)

# Figure 6a
# manually add in the lighter shade of color for statistically insignificant
#  results
png(
    filename = file.path(DATA_OUT, "figure_6a.png"),
    width = 18,
    height = 16.615,
    units = "cm",
    bg = "transparent",
    res = 600
)
ggplot(subset(
    communities_center_shapeshift,
    round(edf) == 1
), aes(
    slope,
    fill = shift_case_num,
    alpha = nonzero_slope
)) +
    geom_histogram(binwidth = 0.05) +
    # linewidth = 0.5 for thinner dashed line
    geom_vline(xintercept = 0, linetype = 2, linewidth = 1) +
    scale_fill_manual(values = c(
        "linear: shift earlier" = unname(palette.colors()[2]),
        "linear: shift later" = unname(palette.colors()[3]),
        "linear: some shift later, others shift earlier" = unname(
            palette.colors()[4]
        )
    ), labels = c(
        "shift earlier",
        "shift later",
        "some shift later, others shift earlier"
    )) +
    scale_alpha_manual(values = c("yes" = 1, "no" = 0.33), guide = "none") +
    coord_cartesian(ylim = c(0, 30)) +
    scale_x_continuous(breaks = seq(-0.1, 0.2, 0.05)) +
    annotate(
        geom = "text",
        x = -0.125,
        y = 29,
        label = paste0(dim(subset(
            communities_center_shapeshift,
            round(edf) == 1
        ))[1]," communities\nwith linear fits"),
        hjust = 0,
        size = 6
    ) +
    labs(
        x = "Slope of linear fit",
        y = "# of communities",
        fill = "direction(s) of shift"
    ) +
    # toggle this to get the transparent legend key items
    # guides(fill = guide_legend(override.aes = list(alpha = 0.33))) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        legend.position = "top",
        legend.direction = "vertical"
    )
dev.off()

# Figure 6c
png(
    filename = file.path(DATA_OUT, "figure_6c.png"),
    width = 18,
    height = 13.5,
    units = "cm",
    bg = "transparent",
    res = 600
)
ggplot(subset(communities_center_shapeshift, round(edf) != 1), aes(
    kendall,
    fill = shift_case_num,
    alpha = correlated
)) +
    geom_histogram(binwidth = 0.1) +
    # linewidth = 0.5625 for thinner dashed line
    geom_vline(xintercept = 0, linetype = 2, linewidth = 1) +
    annotate(
        geom = "text",
        x = -0.25,
        y = 4.5,
        label = paste0(
            dim(subset(communities_center_shapeshift, round(edf) != 1))[1],
            " communities\nwith nonlinear fits"
        ),
        hjust = 0,
        size = 6
    ) +
    scale_fill_manual(values = c(
        "nonlinear: shift earlier" = unname(palette.colors()[2]),
        "nonlinear: shift later" = unname(palette.colors()[3]),
        "nonlinear: some shift later, others shift earlier" = unname(
            palette.colors()[4]
        )
    ), labels = c(
        "shift earlier",
        "shift later",
        "some shift later,\nothers shift earlier"
    ), drop = FALSE) +
    scale_alpha_manual(values = c("yes" = 1, "no" = 0.33), guide = "none") +
    scale_x_continuous(breaks = seq(-0.2, 0.6, 0.1)) +
    labs(
        x = "Kendall's tau (-1 to +1)",
        y = "# of communities",
        fill = "direction of shift"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        legend.position = "top",
        legend.direction = "vertical"
    )
dev.off()

# Figure 6b
png(
    filename = file.path(DATA_OUT, "figure_6b.png"),
    width = 18,
    height = 16.615,
    units = "cm",
    bg = "transparent",
    res = 600
)
ggplot(subset(
    communities_duration_shapeshift,
    round(edf) == 1
), aes(
    slope,
    fill = shift_case_num,
    alpha = nonzero_slope
)) +
    geom_histogram(binwidth = 0.05, show.legend = TRUE) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 1) +
    scale_fill_manual(values = c(
        "linear: shift earlier" = unname(palette.colors()[2]),
        "linear: shift later" = unname(palette.colors()[3]),
        "linear: some shift later, others shift earlier" = unname(palette.colors()[4])
    ), labels = c(
        "shift earlier",
        "shift later",
        "some shift later, others shift earlier"
    )) +
    scale_alpha_manual(values = c("yes" = 1, "no" = 0.33), guide = "none") +
    scale_x_continuous(breaks = seq(-0.1, 0.2, 0.05)) +
    coord_cartesian(ylim = c(0, 60)) +
    annotate(
        geom = "text",
        x = -0.125,
        y = 57.5,
        label = paste0(dim(subset(
            communities_duration_shapeshift,
            round(edf) == 1
        ))[1]," communities\nwith linear fits"),
        hjust = 0,
        size = 6
    ) +
    labs(
        x = "Slope of linear fit",
        y = "# of communities",
        fill = "direction(s) of shift"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        legend.position = "top",
        legend.direction = "vertical"
    )
dev.off()

# Figure 6d
png(
    filename = file.path(DATA_OUT, "figure_6d.png"),
    width = 18,
    height = 13.5,
    units = "cm",
    bg = "transparent",
    res = 600
)
ggplot(subset(communities_duration_shapeshift, round(edf) != 1), aes(
    kendall,
    fill = shift_case_num,
    alpha = correlated
)) +
    geom_histogram(binwidth = 0.1) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 1) +
    annotate(
        geom = "text",
        x = -0.45,
        y = 7.5,
        label = paste0(
            dim(subset(communities_duration_shapeshift, round(edf) != 1))[1],
            " communities\nwith nonlinear fits"
        ),
        hjust = 0,
        size = 6
    ) +
    scale_fill_manual(values = c(
        "nonlinear: shift earlier" = unname(palette.colors()[2]),
        "nonlinear: shift later" = unname(palette.colors()[3]),
        "nonlinear: some shift later, others shift earlier" = unname(palette.colors()[4])
    ), labels = c(
        "shift earlier",
        "shift later",
        "some shift later,\nothers shift earlier"
    )) +
    scale_alpha_manual(values = c("yes" = 1, "no" = 0.33), guide = "none") +
    scale_x_continuous(breaks = seq(-0.4, 0.3, 0.1)) +
    labs(
        x = "Kendall's tau (-1 to +1)",
        y = "# of communities",
        fill = "shape of shift",
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2)),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        legend.position = "top",
        legend.direction = "vertical"
    )
dev.off()

# Appendix S1: Figure S8
png(
    filename = file.path(DATA_OUT, "figure_s8.png"),
    width = 8.5,
    height = 5,
    units = "cm",
    bg = "transparent",
    pointsize = 10,
    res = 600
)
ggplot(
    communities_duration_shapeshift,
    aes(
        x = duration_phenophase_range,
        y = kendall,
        colour = correlated,
        shape = dataset_abbrev,
        size = num_sp
    )
) +
    geom_point(alpha = 0.66) +
    labs(
        x = "Range of phenophase duration",
        y = "Kendall's tau",
        colour = "days per decade shifted\ncorrelated with duration\n(significance at 0.05 level)",
        shape = "dataset",
        size = "scale for species richness"
    ) +
    scale_shape_manual(values = c(17, 16, 15), guide = guide_legend(
        order = 1
    )) +
    scale_size(
        breaks = c(6, 28, 56),
        labels = c("6 species", "28 species", "56 species"),
        range = c(0.5, 2.5),
        guide = guide_legend(order = 3)
    ) +
    scale_colour_manual(
        values = rev(unname(palette.colors(n = 2))),
        guide = guide_legend(order = 2)
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(
            size = rel(0.6),
            margin = margin(0, 0, 0, 0, unit = "cm")
        ),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.margin = margin(-0.15, 0, -0.15, 0, unit = "cm"),
        legend.text = element_text(size = rel(0.6)),
        legend.key.spacing.y = unit(-2.5, 'mm'),
        axis.title = element_text(size = rel(0.7)),
        axis.text = element_text(size = rel(0.7))
    )
dev.off()

# Save to .RData files (as needed) ---------------------------------------------
save(communities, center_D, duration_D, file = file.path(
    DATA_OUT_RDATA,
    "real_community_phenology.RData"
))
