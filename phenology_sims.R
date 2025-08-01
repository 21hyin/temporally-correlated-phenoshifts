# Run with R version 4.2.0
# This script has labeled code sections for easier navigation
# All figures and tables are commented as well, just search for one desired
#  (ex: "Figure 3")

library(ggplot2) # version 3.5.1
library(pracma) # version 2.4.4
library(dplyr) # version 1.1.4
# library(plyr)
library(purrr) # version 1.0.2
library(reshape2) # version 1.4.4

# Specify directories (change as needed)
DATA_IN <- "datasets"
DATA_OUT <- "out"
DATA_OUT_RDATA <- "out/saves_RData"

# Quick start using .RData files to reproduce figures & tables -----------------
# after loading, it's possible to run only portions of the script
# ctrl+f or look through labeled code sections to reproduce specific results
# it's also possible for data frames in the RData files to differ from
#  those produced by running the script from scratch, due to randomness in
#  sampling procedures; however, averages computed for our results should stay
#  robust to this randomness
load(file.path(DATA_OUT_RDATA, "anuran_phenology.RData"))
load(file.path(DATA_OUT_RDATA, "phenology_sims_samples.RData"))
load(file.path(DATA_OUT_RDATA, "phenology_sims_samples_variance.RData"))
load(file.path(DATA_OUT_RDATA, "anuran_phenology_sims.RData"))
load(file.path(DATA_OUT_RDATA, "synthetic_phenology_sims.RData"))


# General code applicable across simulations & outputs -------------------------
# ggplot2 settings
mytheme <- theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(colour = "black")
)

# species_names() returns the genus abbreviation and full species name tied to
#  a vector of species abbreviations
species_names <- function(sp_col) {
    factor(
        sp_col,
        levels = c(
            "AC",
            "BV",
            "BW",
            "GC",
            "HC",
            "HV",
            "PC",
            "PT",
            "RCA",
            "RCL",
            "RP",
            "RS"
        ),
        labels = c(
            "A. crepitans",
            "B. valliceps",
            "B. woodhouseii",
            "G. carolinensis",
            "H. cineria",
            "H. versicolor",
            "P. crucifer",
            "P. triseriata",
            "R. catesbeiana",
            "R. clamitans",
            "R. palustris",
            "R. sphenocephala"
        )
    )
}

# fshift_x() returns a function that shifts function f by h units
# h > 0: shifts function f to the right by h units
# h < 0: shifts function f to the left by h units
# Helps create variation in phenophases in empirical community simulations and
#  simulate phenological shifts
fshift_x <- function(f, h) {
    return (function(x) {f(x - h)})
}

# interaction_potential() returns a data frame of interaction potentials, given
#  a list of probability density functions pdfs and the time interval indicated
#  by days_per_integration
#  returned data frame columns include
#   sp: species
#   lower lim, upper lim: range of days of year for each time interval
#    interaction potential gets computed for
#   interval: numbers the time intervals in chronological order
#    proportion of integrated areas
#   prob: integrated area or probability of phenophase,
#    denoted a_i in manuscript
#   prop_prob_squared, mean_area, evenness_simpson: intermediate computations
#    for interaction potential
#   ip: interaction potential
interaction_potential <- function(
    pdfs,
    days_per_integration = 5
) {
    limits_of_integration <- seq.int(
        from = 1,
        to = 366,
        by = days_per_integration
    )
    fiveday <- data.frame(
        lower_lim = rep.int(
            limits_of_integration[-length(limits_of_integration)],
            times = length(pdfs)
        ),
        upper_lim = rep.int(
            limits_of_integration[-1],
            times = length(pdfs)
        ),
        interval = rep.int(
            1:(length(limits_of_integration) - 1),
            times = length(pdfs)
        )
    )
    probs <- sapply(pdfs, function(kde) {
        map2(
            limits_of_integration[-length(limits_of_integration)],
            limits_of_integration[-1],
            function(a, b, pdf = kde) {
                if (all(pdf(a:b) == 0)) {
                    return (0)
                }
                auc <- integral(
                    pdf,
                    xmin = a,
                    xmax = b,
                    reltol = .Machine$double.eps^0.25
                    # arguements for troubleshooting
                    # subdivisions = 200L,
                    # rel.tol = .Machine$double.eps^0.26
                    # stop.on.error = FALSE
                )
                # cap at 1 due to approximation
                return (min(c(auc, 1)))
            }
        )
    })
    fiveday$prob <- unlist(probs)
    f_prop_prob_squared <- function(p) {
        if (sum(p) > 0) {
            return ((p / sum(p)) ^ 2)
        } else {
            return (0)
        }
    }
    fiveday <- mutate(
        group_by(fiveday, interval),
        prop_prob_squared = f_prop_prob_squared(prob)
    )
    f_evenness_simpson <- function(p, l) {
        if (sum(p) > 0) {
            return (1 / (l * sum(p)))
        } else {
            return (0)
        }
    }
    fiveday_interaction_tibble <- summarize(
        group_by(fiveday, interval),
        mean_area = mean(prob),
        evenness_simpson = f_evenness_simpson(prop_prob_squared, length(pdfs))
    )
    return(mutate(
        group_by(fiveday_interaction_tibble, interval),
        ip = mean_area * evenness_simpson
    ))
}

sample_size <- 100
set.seed(0)

# Empirical community simulation setup (Appendix S1: Fig. S2, part of Fig. 2 & Appendix: Fig. S2, Appendix S1: Fig. S9) ----
NIGHTLYCALLS <- read.csv(file.path(DATA_IN, "nightlycalls.csv"), header = TRUE)
str(NIGHTLYCALLS)
NIGHTLYCALLS$pond <- as.factor(NIGHTLYCALLS$pond)
NIGHTLYCALLS$date <- as.Date(strptime(NIGHTLYCALLS$date, "%Y-%m-%d"))
# technically the year starts in 2000 and not in 2001
# however, sampling always straddles two days since monitoring occurs over 6
#  hours with midnight at the midpoint; the day of sampling always refers to the
#  day the sampling started, so the first sampling day of 2000-12-31 crosses
#  over into new years day, and the last sampling day on 2022-12-30 makes up a
#  complete set of days across years with no more or less days to spare
min(NIGHTLYCALLS$date)
max(NIGHTLYCALLS$date)

# Appendix S1: Figure S2
nightlycalls <- NIGHTLYCALLS
nightlycalls$species_name <- species_names(nightlycalls$sp)
png(
    filename = file.path(DATA_OUT, "figure_s2.png"),
    width = 8.5,
    height = 10,
    units = "cm",
    pointsize = 10,
    res = 300
)
ggplot(subset(nightlycalls, subset = (nightlysum > 0)), aes(doy)) +
    mytheme +
    stat_density(
        alpha = 0.75,
        adjust = 1/5,
        aes(ymax = after_stat(density),  ymin = -after_stat(density)),
        fill = "darkslategray",
        colour = "darkslategray",
        geom = "ribbon",
        position = "identity"
    ) +
    facet_grid(species_name ~.  , switch = "y") +
    scale_x_continuous(breaks = seq(0, 350, 50)) +
    theme(
        axis.title = element_text(size = 10),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linewidth = 1),
        strip.text = element_text(size = 8),
        strip.text.y = element_text(angle = 180, face = "italic"),
        strip.text.y.left = element_text(angle = 0)
    ) +
    ylab(NULL) +
    xlab("Calendar day of year")
dev.off()

# Calls between days of year 1-21 and 289-366 to be ignored
#  based on calling activity of Rana sphenocephala (see Appendix S1: Section S1
#  "Creating temporal structure of empirical communities")
window_first <- 289
window_last <- 21

# First day of re-centered calendar year
#  relative to R. sphenocephala calling activity
frog_doy1 <- floor((window_first + 366 + window_last) / 2)

# Split RS into "spring" and "fall" using calendar day of year 165 (see
#  Appendix S1: Section S1 "Creating temporal structure of empirical
#  communities")
nightlycalls_RS_spring <- subset(NIGHTLYCALLS, subset = (sp == "RS"))
nightlycalls_RS_spring$sp <- "RS_spring"
nightlycalls_RS_spring[
    nightlycalls_RS_spring$doy >= 165 & nightlycalls_RS_spring$doy <= frog_doy1,
]$nightlysum <- 0
nightlycalls_RS_fall <- subset(NIGHTLYCALLS, subset = (sp == "RS"))
nightlycalls_RS_fall$sp <- "RS_fall"
nightlycalls_RS_fall[
    nightlycalls_RS_fall$doy < 165 | nightlycalls_RS_fall$doy >= frog_doy1,
]$nightlysum <- 0
nightlycalls_RS_split <- rbind(
    NIGHTLYCALLS,
    nightlycalls_RS_spring,
    nightlycalls_RS_fall
)

# Re-center calendar year relative to R. sphenocephala calling activity
nightlycalls_RS_split$frog_doy <- as.numeric(
    strftime(
        nightlycalls_RS_split$date + (366 - frog_doy1),
        format = "%j"
    )
)
nightlycalls_RS_split$frog_year <- as.numeric(
    strftime(
        nightlycalls_RS_split$date + (366 - frog_doy1),
        format = "%Y"
    )
)

# More pruning before getting calling probability density functions

# cumsumxNA() returns the cumulative sum of vector x
#  while ignoring NAs in the vector
cumsumxNA <- function(x) {
    x[which(is.na(x))] <- 0
    return(cumsum(x))
}
nightlycalls_RS_split <- mutate(
    group_by(nightlycalls_RS_split, pond, sp, frog_year),
    frog_cumsum = cumsumxNA(nightlysum)
)
# plyr version
# nightlycalls_RS_split <- ddply(
#     nightlycalls_RS_split,
#     .(pond, sp, frog_year),
#     transform,
#     frog_cumsum  = cumsumxNA(nightlysum)
# )

# find years with complete data relative to annual breeding season
#  which we re-centered our data on
first_frog_year <- min(nightlycalls_RS_split$frog_year)
first_full_frog_year <- first_frog_year
if(min(nightlycalls_RS_split[
    nightlycalls_RS_split$frog_year == first_frog_year,
]$frog_doy) != 1) {
    first_full_frog_year <- first_frog_year + 1
}
last_frog_year <- max(nightlycalls_RS_split$frog_year)
last_full_frog_year <- last_frog_year
if(max(nightlycalls_RS_split[
    nightlycalls_RS_split$frog_year == last_frog_year,
]$frog_doy) < 365) {
    last_full_frog_year <- last_frog_year - 1
}

d <- subset(
    nightlycalls_RS_split,
    select = c(
        "sp",
        "nightlysum",
        "pond",
        "date",
        "doy",
        "frog_doy",
        "frog_year",
        "frog_cumsum"
    )
)
d <- subset(d, subset = (nightlysum > 0))

d_frog <- d[
    d$frog_year <= last_full_frog_year & d$frog_year >= first_full_frog_year,
]

# use calls to ignored between calendar days of year 1-21 and 289-366
#  as starting point to reduce more noise
#  in first and last 1% of cumulative calls
#  within each species-pond-year combination
nightlycalls_RS_split$nightlysum_window <- nightlycalls_RS_split$nightlysum
nightlycalls_RS_split$nightlysum_window[which(
    nightlycalls_RS_split$doy >= window_first | nightlycalls_RS_split$doy <= window_last
)] <- 0
nightlycalls_RS_split <- mutate(
    group_by(nightlycalls_RS_split, pond, sp, frog_year),
    frog_cumsum_window  = cumsumxNA(nightlysum_window)
)
# plyr version
# nightlycalls <- ddply(
#     nightlycalls,
#     .(pond, sp, frog_year),
#     transform,
#     frog_cumsum_window  = cumsumxNA(nightlysum_window)
# )

# find days when calls make up the first and last 1% of cumulative calls
#  within each species-pond-year combination and get median day
# warnings are ok, accounted for right after
callperiods_frog <- dplyr::summarize(
    group_by(dplyr::arrange(
        group_by(nightlycalls_RS_split, pond, sp, frog_year),
        frog_doy,
        .by_group = TRUE
    ), pond, sp, frog_year),
    firstdate = frog_doy[min(which(
        frog_cumsum_window > (max(frog_cumsum_window) * 0.01)
    ))],
    meddate = frog_doy[min(which(
        frog_cumsum_window >= (max(frog_cumsum_window) / 2)
    ))],
    lastdate = frog_doy[max(which(
        frog_cumsum_window < (max(frog_cumsum_window) * 0.99)
    ))]
)
# correct issues created by NAs
callperiods_frog$meddate[callperiods_frog$firstdate %in% NA] <- NA
callperiods_frog$duration <- as.numeric(
    callperiods_frog$lastdate - callperiods_frog$firstdate
)
callperiods_frog <- na.omit(callperiods_frog)

p_frog <- subset(callperiods_frog, subset = (duration > 0))

# exclude days when calls make up the first and last 1% of cumulative calls
#  within each species-pond-year combination (see Appendix S1: Section S1
#  "Creating temporal structure of empirical communities")
d_frog <- mutate(
    # d_frog has nightlysum > 0, more efficient than using the full nightlycalls
    #   dataframe
    # left join (
    #   left dataframe = abundances,
    #   right dataframe = callperiods_frog
    # )
    #  works but leaves sp-pond-years with nightlysum = 0 as NA if using
    #   nightlycalls
    # merge()'s default inner join avoids extra step of removing NA post join
    group_by(merge(d_frog, p_frog), pond, sp, frog_year),
    nightlysum_trim = ifelse(
        frog_doy < firstdate | frog_doy > lastdate,
        0,
        nightlysum
    )
)
d_frog <- subset(d_frog, subset = nightlysum_trim > 0)

# Create data frame to store data on calling probability density functions
#  using the Gaussian kernel

# dplyr version
base_pdfs <- summarise(
    group_by(d_frog, sp),
    firstdate = min(frog_doy),
    meddate = quantile(frog_doy, probs = c(0.5)),
    lastdate = max(frog_doy)
)
# plyr version
# base_pdfs <- ddply(
#     d_frog,
#     .(sp),
#     summarise,
#     firstdate = min(frog_doy),
#     meddate = quantile(frog_doy, probs = c(0.5)),
#     lastdate = max(frog_doy)
# )
base_pdfs$duration <- base_pdfs$lastdate - base_pdfs$firstdate
base_pdfs$kde <- tapply(d_frog$frog_doy, d_frog$sp, function(frog_doy) {
    return(approxfun(density(frog_doy, bw = 5, cut = 0), yleft = 0, yright = 0))
})

# Part of Figure 2 & Figure S1
plot_base_pdfs <- apply(
    subset(base_pdfs, subset = (sp != "RS_spring" & sp != "RS_fall")),
    1,
    function(row) { return (data.frame(
        sp = row["sp"],
        frog_doy = 1:366,
        kde = row["kde"][[1]](1:366)
    ))}
)
plot_base_pdfs <- bind_rows(plot_base_pdfs)
plot_base_pdfs$species_name <- species_names(plot_base_pdfs$sp)
png(
    filename = file.path(DATA_OUT, "figure_2_base_calling_pdfs.png"),
    width = 18,
    height = 18,
    units = "cm",
    res = 600,
    bg = "transparent"
)
ggplot(
    plot_base_pdfs,
    aes(x = frog_doy, y = kde)
) +
    geom_area(fill = "gray", position = 'identity') +
    geom_line(colour = "black") +
    facet_grid(species_name ~ ., switch = "y") +
    scale_x_continuous(breaks = seq(0, 366, 50)) +
    scale_y_continuous(position = "right") +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line.x  = element_line(),
        axis.title.x = element_text(size = rel(1.5)),
        axis.text.x  = element_text(size = rel(1.5)),
        panel.spacing.y = unit(10, "pt"),
        strip.text.y.left = element_text(angle = 0),
        strip.text.y = element_text(size = rel(1.5), face = "italic"),
        axis.ticks.y = element_blank()
    ) +
    ylab("Density") +
    xlab("Day of year")
dev.off()

# Select the data from which to subsample median days and durations for
#  realistic phenophases that vary in their timing and duration

# high variation
# get extremes only 0.02 through (0.02 + 0.5 / 2), (0.98 - 0.5 / 2) through 0.98
p_frog_extremes50 <- filter(
    group_by(p_frog, sp),
    (
        (
            (meddate >= quantile(meddate, probs = c(0.02))) & (
                meddate < quantile(meddate, probs = c(0.02 + 0.5 / 2))
            )
        ) | (
            (meddate > quantile(meddate, probs = c(0.98 - 0.5 / 2))) & (
                meddate <= quantile(meddate, probs = c(0.98))
            )
        )
    ) & (
        (
            (duration >= quantile(duration, probs = c(0.02))) & (
                duration < quantile(duration, probs = c(0.02 + 0.5 / 2))
            )
        ) | (
            (duration > quantile(duration, probs = c(0.98 - 0.5 / 2))) & (
                duration <= quantile(duration, probs = c(0.98))
            )
        )
    )
)
# low variation (did not end up using for main results)
p_frog_mid50 <- filter(
    group_by(p_frog, sp),
    ((meddate >= quantile(
        meddate,
        probs = c(0.25)
    )) & (meddate <= quantile(meddate, probs = c(0.75)))) & ((
        duration >= quantile(duration, probs = c(0.25))) & (
            duration <= quantile(duration, probs = c(0.75))
        )
    )
)
# determines the number of days phenophases shift and helps set bounds for
#  synthetic communities in next section
p_frog_mid96 <- filter(
    group_by(p_frog, sp),
    ((meddate >= quantile(meddate, probs = c(0.02))) & (
        meddate <= quantile(meddate, probs = c(0.98))
    )) & ((duration >= quantile(duration, probs = c(0.02))) & (
        duration <= quantile(duration, probs = c(0.98))
    ))
)

# Sample median days and durations

# no simple way of checking repeats, but the chance of getting repeats is
#  extremely low (< 2 ^ -31)
#  exact probability of getting repeats from p_frog_extremes50 would be:
#  1 / (
#       prod(sapply(group_rows(group_by(
#           p_frog_extremes50[p_frog_extremes50$sp != "RS", ],
#           sp
#       )), length))
#  )
p_frog_extremes50_samples <- slice_sample(
    group_by(
        subset(p_frog_extremes50, subset = (sp != "RS")),
        sp
    ),
    n = sample_size,
    replace = TRUE
)
p_frog_extremes50_samples$var <- "high"
p_frog_extremes50_samples$scheme <- "extremes50"

p_frog_mid50_samples <- slice_sample(
    group_by(
        subset(p_frog_mid50, subset = (sp != "RS")),
        sp
    ),
    n = sample_size,
    replace = TRUE
)
p_frog_mid50_samples$var <- "low"
p_frog_mid50_samples$scheme <- "mid50"

calling_phenophases <- rbind(
    p_frog_mid50_samples,
    p_frog_extremes50_samples
)

# enumerate trials where each trial is a community
calling_phenophases <- mutate(
    group_by(calling_phenophases, sp, var, scheme),
    trial = seq_along(sp)
)

# use var as the standardized column name for variation in phenophases
calling_phenophases <- calling_phenophases[,
    !(names(calling_phenophases) %in% c("scheme"))
]

# Appendix S1: Figure S9
calling_phenophases_corr <- summarise(
    group_by(calling_phenophases[calling_phenophases$var == "high", ], trial),
    kendall_corr = unname(
        cor.test(meddate, duration, method = "kendall")$estimate
    )
)
png(
    filename = file.path(DATA_OUT, "figure_s9.png"),
    width = 8.5,
    height = 8.5,
    units = "cm",
    pointsize = 10,
    res = 600
)
hist(
    calling_phenophases_corr$kendall_corr,
    main = "",
    xlab = "Kendall rank correlation coefficient",
    ylab = "# of anuran communities",
    xlim = c(-0.4, 0.8),
    ylim = c(0, 25)
)
dev.off()

# mean reported in Appendix S1: Section S2
mean(calling_phenophases_corr$kendall_corr)

# Get the refined shapes of the calling phenophase activities in communities
#  with high variation (sampled from p_frog_extremes50)

# pdf_resize() returns the re-scaled function of function f_pdf
#  where f_pdf(x) = 0 when
#  x < center - range_new / 2 and x > center + range_new / 2
#  for 0 <= x <= 367
# used to re-scale calling probability distributions and
# equivalent to Equation # in Appendix S1
pdf_resize <- function(f_pdf, center, range_new) {
    max_xrange <- 0:367
    nonzeros_max_xrange <- which(f_pdf(max_xrange) > 0)
    first_day_original <- max_xrange[min(nonzeros_max_xrange)]
    duration_original <- max_xrange[
        max(nonzeros_max_xrange) - 1
    ] - first_day_original + 1
    first_day_resized <- -1 * ((
        center - first_day_original
    ) / duration_original * range_new - center)
    return (function(x) {duration_original / range_new * f_pdf(
        first_day_original + duration_original / range_new * (
            x - first_day_resized
        )
    )})
}
# merge_pdfs() returns a function that outputs the maximum values of the list of
#  functions pdfs across the domain 1 through 366
# great for merging calling probability distributions of spring and fall
#  R. sphenocephala
merge_pdfs <- function(pdfs) {
    pdf_outputs <- as.data.frame(cbind(sapply(pdfs, function(f) f(1:366))))
    f <- function(x) {
        if (x < 1 || x > 366) {
            return (0)
        } else {
            probs <- apply(pdf_outputs, 1, max)
            # normalize to integrate to 1
            probs <- probs / sum(probs)
            return (probs[x])
        }
    }
    return (Vectorize(f, USE.NAMES = FALSE))
}

# create realistic calling phenophases based on observed instances of calling at
#  ponds and years per species
calling_phenophases$kde <- apply(
    calling_phenophases[, c("sp", "meddate", "duration")],
    1,
    function(x) { return (
        pdf_resize(
            fshift_x(
                base_pdfs$kde[[x["sp"]]],
                as.numeric(x["meddate"]) - base_pdfs[
                    base_pdfs[, "sp"] == x["sp"],
                    "meddate"
                ]
            ),
            as.numeric(x["meddate"]),
            as.numeric(x["duration"])
        )
    )}
)
# get the R. sphenocephala calling distribution for the full year
calling_phenophases_RS <- subset(calling_phenophases, subset = (
    sp == "RS_spring" | sp == "RS_fall"))
calling_phenophases_RS <- summarise(
    group_by(calling_phenophases_RS, trial, var),
    kde = c(merge_pdfs(kde))
)
calling_phenophases_RS$sp <- "RS"

calling_phenophases_no_shift <- rbind(
    subset(
        calling_phenophases,
        select = c("trial", "var", "kde", "sp"),
        subset = (sp != "RS" & sp != "RS_spring" & sp != "RS_fall")
    ),
    calling_phenophases_RS
)
calling_phenophases_no_shift <- calling_phenophases_no_shift[
    calling_phenophases_no_shift$var == "high",
]

# Synthetic community simulation setup (Appendix S1: Fig. S5) ------------------
# Create data frame to store phenophases

num_sp <- length(unique(NIGHTLYCALLS$sp))

synthetic_phenophases <- expand.grid(
    trial = 1:sample_size,
    sp = LETTERS[1:num_sp],
    center_var = c("low", "high"),
    spread_var = c("low", "high")
)

# rnsbeta() returns random samples from the non-standard Beta distribution
#  within the range defined by the min and max parameters
# wrapper function around rbeta() in R stats package
# n: sample size; if length(n) > 1, the length is taken to be the number
#  required
# alpha, beta: non-negative shape parameters of the Beta distribution
rnsbeta <- function(n, alpha, beta, min = 0, max = 1) {
    stopifnot(min < max, alpha >= 0, beta >= 0)
    return (rbeta(n, shape1 = alpha, shape2 = beta) * (max - min) + min)
}

# specify symmetric Beta distributions to sample from to create high and low
#  temporal variation in phenophases within communities
beta_sym_params <- list(
    # high: uniform distribution
    #  (middle 50% of samples covers 50% of probability distribution)
    high = rep.int(1, times = 2),
    # low: middle 50% of samples covers >99% of probability distribution
    low = rep(11.8, times = 2)
)

# bound possible median days and durations to those used in
#  empirical community simulations and in terms of mean and standard deviation
#  of truncated Normal distribution used to represent phenophase distributions
center_min <- min(subset(p_frog_mid96, subset = (sp != "RS"))$meddate)
center_max <- max(subset(p_frog_mid96, subset = (sp != "RS"))$meddate)
spread_min <- min(subset(p_frog_mid96, subset = (sp != "RS"))$duration) / 6
spread_max <- max(subset(p_frog_mid96, subset = (sp != "RS"))$duration) / 6

# sample phenophases
synthetic_phenophases <- mutate(
    group_by(synthetic_phenophases, center_var, spread_var),
    center = rnsbeta(
        n = num_sp * sample_size,
        alpha = beta_sym_params[[unique(as.character(center_var))]][1],
        beta = beta_sym_params[[unique(as.character(center_var))]][2],
        min = center_min,
        max = center_max
    ),
    spread = rnsbeta(
        n = num_sp * sample_size,
        alpha = beta_sym_params[[unique(as.character(spread_var))]][1],
        beta = beta_sym_params[[unique(as.character(spread_var))]][2],
        min = spread_min,
        max = spread_max
    )
)
synthetic_phenophases$kde <- apply(
    synthetic_phenophases[, c("center", "spread")],
    1,
    function(params) {
        # some truncated and normalized version of dnorm
        a <- params["center"] - 3 * params["spread"]
        b <- params["center"] + 3 * params["spread"]
        # make sure truncation integrates to 1
        c <- 1 / (
            pnorm(
                q = b,
                mean = params["center"],
                sd = params["spread"]
            ) - pnorm(q = a, mean = params["center"], sd = params["spread"])
        )
        return (approxfun(
            x = a:b,
            y = c * dnorm(a:b, mean = params["center"], sd = params["spread"]),
            yleft = 0,
            yright = 0
        ))
    }
)

# Appendix S1: Figure S5
synthetic_phenophases_xy <- apply(
    synthetic_phenophases,
    1,
    function(phenophase) { return (data.frame(
        trial = phenophase["trial"],
        sp = phenophase["sp"],
        center_var = phenophase["center_var"],
        spread_var = phenophase["spread_var"],
        frog_doy = 1:366,
        kde = phenophase["kde"][[1]](1:366)
    ))}
)
synthetic_phenophases_xy <- bind_rows(synthetic_phenophases_xy)
synthetic_phenophases_xy <- synthetic_phenophases_xy[
    synthetic_phenophases_xy$trial <= 3,
]
synthetic_phenophases_xy$trial <- as.factor(synthetic_phenophases_xy$trial)
centers_spreads <- list(
    c("high", "high"),
    c("high", "low"),
    c("low", "high"),
    c("low", "low")
)
synthetic_phenophase_samples_plots <- map(
    centers_spreads,
    function(cs) {
        ggplot(
            synthetic_phenophases_xy[
                synthetic_phenophases_xy$center_var == cs[1] & synthetic_phenophases_xy$spread_var == cs[2],
            ],
            aes(x = frog_doy, y = kde, colour = trial)
        ) +
            mytheme +
            geom_line() +
            facet_grid(sp ~ ., scales = "free_y") +
            scale_y_continuous(n.breaks = 3) +
            theme(
                plot.background = element_rect(
                    fill = "transparent",
                    color = NA
                ),
                axis.ticks = element_blank(),
                axis.title = element_text(size = rel(0.75)),
                axis.text.x  = element_text(size = rel(0.75)),
                axis.text.y = element_text(size = rel(0.5)),
                legend.title = element_text(
                    size = rel(0.75),
                    margin = margin(0, 0, 0, 0, unit = "cm")
                ),
                legend.background = element_rect(
                    fill = "transparent",
                    color = NA
                ),
                legend.text = element_text(size = rel(0.75)),
                legend.key.spacing.y = unit(-2, 'mm'),
                strip.text.y = element_text(size = rel(0.75))
            ) +
            scale_colour_manual(values = unname(palette.colors(n = 3))) +
            guides(colour = guide_legend(override.aes = list(
                linewidth = 1.05
            ))) +
            labs(colour = "community id") +
            ylab("Probability density") +
            xlab("Day of year")
    }
)
synthetic_phenophase_samples_paths <- map(
    centers_spreads,
    \(cs) {
        file.path(
            DATA_OUT,
            paste0(
                "figure_s5_",
                cs[1],
                "center_",
                cs[2],
                "spread.png"
            )
        )
    }
)
walk2(
    synthetic_phenophase_samples_paths,
    synthetic_phenophase_samples_plots,
    \(path, plot) {
        png(
            filename = path,
            width = 8.5,
            height = 10.2,
            units = "cm",
            res = 600,
            bg = "transparent"
        )
        print(plot)
        dev.off()
    }
)

# Variance in timing & duration of phenophases (Appendix S1: Table S1) ---------
calling_phenophases_var_stats <- summarise(
    group_by(
        summarise(
            group_by(
                subset(calling_phenophases, subset = (sp != "RS")),
                trial,
                var
            ),
            meddate_var = var(meddate),
            duration_var = var(duration)
        ),
        var
    ),
    meddate_meanvar = mean(meddate_var),
    duration_meanvar = mean(duration_var)
)

synthetic_phenophases_var_stats <- summarise(
    group_by(summarise(
        group_by(synthetic_phenophases, trial, center_var, spread_var),
        center_v = var(center),
        duration_v = var(spread * 6)
    ), center_var, spread_var),
    center_meanvar = mean(center_v),
    duration_meanvar = mean(duration_v)
)

# Initial interaction potentials -----------------------------------------------
# note: this section can take several minutes

# Empirical communities

calling_ip_no_shift <- map(
    group_split(group_by(calling_phenophases_no_shift, trial)),
    function(group) {
        ip_trial <- interaction_potential(group$kde)
        ip_trial$trial <- group$trial[1]
        return (ip_trial)
    }
)
calling_ip_no_shift <- bind_rows(calling_ip_no_shift)
calling_ip_no_shift$var <- "high"

calling_annotate_box <- summarize(
    ungroup(calling_ip_no_shift),
    q75 = quantile(ip, probs = c(0.75)),
    q75_min_interval = min(interval[which(ip >= q75)]),
    q75_max_interval = max(interval[which(ip >= q75)])
)

# Synthetic communities

synthetic_ip_no_shift <- map(
    group_split(group_by(
        synthetic_phenophases,
        trial,
        center_var,
        spread_var
    )),
    function(group) {
        ip_trial <- interaction_potential(group$kde)
        ip_trial$trial <- group$trial[1]
        ip_trial$center_var <- group$center_var[1]
        ip_trial$spread_var <- group$spread_var[1]
        return (ip_trial)
    }
)
synthetic_ip_no_shift <- bind_rows(synthetic_ip_no_shift)

synthetic_annotate_box <- summarize(
    group_by(synthetic_ip_no_shift, center_var, spread_var),
    q75 = quantile(ip, probs = c(0.75)),
    q75_min_interval = min(interval[which(ip >= q75)]),
    q75_max_interval = max(interval[which(ip >= q75)])
)

# Plot initial interaction potentials (part of Fig. 4) -------------------------
# Figure 4a left side
g <- ggplot(
    calling_ip_no_shift,
    aes(
        x = interval * 5,
        y = ip
    )
) +
    geom_line(
        mapping = aes(group = trial),
        colour = "black",
        alpha = 0.15
    ) +
    geom_rect(
        data = calling_annotate_box,
        mapping = aes(
            xmin = q75_min_interval * 5,
            xmax = q75_max_interval * 5,
            ymin = -Inf,
            ymax = Inf,
            x = NULL,
            y = NULL
        ),
        alpha = 0.25,
        fill = "gray"
    ) +
    scale_x_continuous(
        breaks = c(1, 122, 244, 366),
        labels = as.character(c(1, 122, 244, 366)),
        limits = c(1, 366)
    ) +
    labs(
        x = "Day of year (community-adjusted)",
        y = "Initial interaction potential"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(2.25)),
        axis.text = element_text(size = rel(2.25))
    )
ggsave(
    plot = g,
    filename = file.path(DATA_OUT, "figure_4a_left.png"),
    width = 18,
    height = 18,
    units = "cm",
    dpi = 600,
    bg = "transparent"
)
# Figure 4b left side
g <- ggplot(
    subset(synthetic_ip_no_shift, center_var == "high" & spread_var == "high"),
    aes(
        x = interval * 5,
        y = ip
    )
) +
    geom_line(
        mapping = aes(group = trial),
        colour = "black",
        alpha = 0.15
    ) +
    geom_rect(
        data = subset(
            synthetic_annotate_box, center_var == "high" & spread_var == "high"
        ),
        mapping = aes(
            xmin = q75_min_interval * 5,
            xmax = q75_max_interval * 5,
            ymin = -Inf,
            ymax = Inf,
            x = NULL,
            y = NULL
        ),
        alpha = 0.25,
        fill = "gray"
    ) +
    scale_x_continuous(
        breaks = c(1, 122, 244, 366),
        labels = as.character(c(1, 122, 244, 366)),
        limits = c(1, 366)
    ) +
    labs(
        x = "Day of year",
        y = "Initial interaction potential"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(2.25)),
        axis.text = element_text(size = rel(2.25))
    )
ggsave(
    plot = g,
    filename = file.path(DATA_OUT, "figure_4b_left.png"),
    width = 18,
    height = 18,
    units = "cm",
    dpi = 600,
    bg = "transparent"
)

# Shift phenophases ------------------------------------------------------------
# note: this section can take several minutes

# shift_phenophase_by() returns the number of days a phenophase shifts given
#  x: numeric value of the phenophase's median day or duration
#  predictor: "median" or "duration" to specify what input x refers to,
#   or "none" to indicate an uncorrelated phenological shift
#  case: integer indicating which of the 9 types of phenological shift
#  max_shift: numeric indicating the maximum number of days a phenophase can
#   across all types of phenological shifts examined
#  median_global_max: maximum expected median day of any given phenophase
#  median_global_min: minimum expected median day of any given phenophase
#  duration_global_max: maximum expected duration of any given phenophase
#  duration_global_min: minimum expected duration of any given phenophase
shift_phenophase_by <- function(
    x,
    predictor,
    case,
    max_shift = 24,
    median_global_max = max(subset(
        p_frog_mid96,
        subset = (sp != "RS")
    )$meddate),
    median_global_min = min(subset(
        p_frog_mid96,
        subset = (sp != "RS")
    )$meddate),
    duration_global_max = max(subset(
        p_frog_mid96,
        subset = (sp != "RS")
    )$duration),
    duration_global_min = min(subset(
        p_frog_mid96,
        subset = (sp != "RS")
    )$duration)
) {
    stopifnot(is.numeric(x))
    stopifnot(is.integer(case))
    stopifnot(
        predictor == "duration" || predictor == "median" || predictor == "none"
    )
    # x is median day or duration
    m <- max_shift / (median_global_max - median_global_min)
    h <- median_global_min
    if (predictor == "duration") {
        m <- max_shift / (duration_global_max - duration_global_min)
        h <- duration_global_min
    }
    b <- max_shift
    return (
        switch(
            case,
            m * (x - h),
            m * (x - h) - b,
            m * (x - h) - b / 2,
            -m * (x - h) + b,
            -m * (x - h),
            -m * (x - h) + b / 2,
            runif(n = 1, min = 0, max = b),
            runif(n = 1, min = -b, max = 0),
            runif(n = 1, min = -b / 2, max = b / 2)
        )
    )
}

# Data frames of shifted phenophases
# note: max_shift = 6 not used in main results, for sensitivity analysis

# Empirical communities
# chose only communities with high variation (appendix section)

calling_phenophases_shift_corr <- expand.grid(
    trial = 1:sample_size,
    sp = unique(calling_phenophases_no_shift$sp),
    shift_case_num = 1:6,
    predictor = c("median", "duration"),
    max_shift = c(24, 6)
)
calling_phenophases_shift_uncorr <- expand.grid(
    trial = 1:sample_size,
    sp = unique(calling_phenophases_no_shift$sp),
    shift_case_num = 7:9,
    predictor = c("none"),
    max_shift = c(24, 6)
)

calling_phenophases_shift_corr_no_RS <- subset(
    calling_phenophases_shift_corr,
    subset = (sp != "RS")
)
calling_phenophases_shift_corr_no_RS$kde <- apply(
    calling_phenophases_shift_corr_no_RS,
    1,
    function(x, x_trial_sp = calling_phenophases[
        calling_phenophases$var == "high",
    ]) {
        single_sample <- subset(
            x_trial_sp,
            subset = (trial == as.integer(x["trial"]) & sp == x["sp"][1])
        )
        predictor <- "duration"
        if (x["predictor"][1] == "median") {
            predictor <- "meddate"
        }
        return (fshift_x(
            single_sample$kde[[1]],
            # find how much to shift by
            shift_phenophase_by(
                as.integer(single_sample[, predictor][[1]]),
                x["predictor"][1],
                as.integer(x["shift_case_num"]),
                max_shift = as.integer(x["max_shift"])
            )
        ))
    }
)
calling_phenophases_shift_corr_RS <- subset(
    calling_phenophases_shift_corr,
    subset = (sp == "RS")
)
calling_phenophases_shift_corr_RS$kde <- apply(
    calling_phenophases_shift_corr_RS,
    1,
    function(x, x_trial = subset(
        calling_phenophases,
        subset = ((sp == "RS_spring" | sp == "RS_fall") & var == "high")
    )) {
        single_trial <- subset(
            x_trial,
            subset = (trial == as.integer(x["trial"]))
        )
        RS_spring_trial <- single_trial[single_trial$sp == "RS_spring", ]
        RS_fall_trial <- single_trial[single_trial$sp == "RS_fall", ]
        predictor <- "duration"
        if (x["predictor"][1] == "median") {
            predictor <- "meddate"
        }
        return (merge_pdfs(c(
            fshift_x(
                RS_spring_trial$kde[[1]],
                # find how much to shift by
                shift_phenophase_by(
                    as.integer(RS_spring_trial[, predictor][[1]]),
                    x["predictor"][1],
                    as.integer(x["shift_case_num"]),
                    max_shift = as.integer(x["max_shift"])
                )
            ),
            fshift_x(
                RS_fall_trial$kde[[1]],
                # find how much to shift by
                shift_phenophase_by(
                    as.integer(RS_fall_trial[, predictor][[1]]),
                    x["predictor"][1],
                    as.integer(x["shift_case_num"]),
                    max_shift = as.integer(x["max_shift"])
                )
            )
        )))
    }
)
calling_phenophases_shift_corr <- rbind(
    calling_phenophases_shift_corr_no_RS,
    calling_phenophases_shift_corr_RS
)

# repeat for noncorrelated shifts
calling_phenophases_shift_uncorr_no_RS <- subset(
    calling_phenophases_shift_uncorr,
    subset = (sp != "RS")
)
calling_phenophases_shift_uncorr_no_RS$kde <- apply(
    calling_phenophases_shift_uncorr_no_RS,
    1,
    function(x, x_trial_sp = calling_phenophases[
        calling_phenophases$var == "high",
    ]) {
        single_sample <- subset(
            x_trial_sp,
            subset = (trial == as.integer(x["trial"]) & sp == x["sp"][1])
        )
        predictor <- "duration"
        if (x["predictor"][1] == "median") {
            predictor <- "meddate"
        }
        return (fshift_x(
            single_sample$kde[[1]],
            # find how much to shift by
            shift_phenophase_by(
                -1, # param that doesn't actually do anything here
                x["predictor"][1],
                as.integer(x["shift_case_num"]),
                max_shift = as.integer(x["max_shift"])
            )
        ))
    }
)
calling_phenophases_shift_uncorr_RS <- subset(
    calling_phenophases_shift_uncorr,
    subset = (sp == "RS")
)
calling_phenophases_shift_uncorr_RS$kde <- apply(
    calling_phenophases_shift_uncorr_RS,
    1,
    function(x, x_trial = subset(
        calling_phenophases,
        subset = ((sp == "RS_spring" | sp == "RS_fall") & var == "high")
    )) {
        single_trial <- subset(
            x_trial,
            subset = (trial == as.integer(x["trial"]))
        )
        RS_spring_trial <- single_trial[single_trial$sp == "RS_spring", ]
        RS_fall_trial <- single_trial[single_trial$sp == "RS_fall", ]
        predictor <- "duration"
        if (x["predictor"][1] == "median") {
            predictor <- "meddate"
        }
        return (merge_pdfs(c(
            fshift_x(
                RS_spring_trial$kde[[1]],
                # find how much to shift by
                shift_phenophase_by(
                    -1, # param that doesn't actually do anything here
                    x["predictor"][1],
                    as.integer(x["shift_case_num"]),
                    max_shift = as.integer(x["max_shift"])
                )
            ),
            fshift_x(
                RS_fall_trial$kde[[1]],
                # find how much to shift by
                shift_phenophase_by(
                    -1, # param that doesn't actually do anything here
                    x["predictor"][1],
                    as.integer(x["shift_case_num"]),
                    max_shift = as.integer(x["max_shift"])
                )
            )
        )))
    }
)
calling_phenophases_shift_uncorr <- rbind(
    calling_phenophases_shift_uncorr_no_RS,
    calling_phenophases_shift_uncorr_RS
)

calling_phenophases_shift <- rbind(
    calling_phenophases_shift_corr,
    calling_phenophases_shift_nocorr
)
calling_phenophases_shift$var <- "high"

# Synthetic communities

synthetic_phenophases_shift_corr <- expand.grid(
    trial = 1:sample_size,
    sp = LETTERS[1:num_sp],
    center_var = c("high", "low"),
    spread_var = c("high", "low"),
    shift_case_num = 1:6,
    predictor = c("median", "duration"),
    max_shift = c(24, 6)
)
synthetic_phenophases_shift_corr$kde <- apply(
    synthetic_phenophases_shift_corr,
    1,
    function(x, x_trial_sp = synthetic_phenophases) {
        single_sample <- subset(
            x_trial_sp,
            subset = (
                trial == as.integer(
                    x["trial"]
                ) & sp == x[
                    "sp"
                ][1] & center_var == x[
                    "center_var"
                ][1] & spread_var == x[
                    "spread_var"
                ][1]
            )
        )
        p <- NULL
        if (x["predictor"][1] == "median") {
            p <- single_sample$center
        }
        if (x["predictor"][1] == "duration") {
            p <- single_sample$spread * 6
        }
        return (fshift_x(
            single_sample$kde[[1]],
            # find how much to shift by
            shift_phenophase_by(
                p,
                x["predictor"][1],
                as.integer(x["shift_case_num"]),
                max_shift = as.integer(x["max_shift"])
            )
        ))
    }
)

# repeat for uncorrelated shifts
synthetic_phenophases_shift_uncorr <- expand.grid(
    trial = 1:sample_size,
    sp = LETTERS[1:num_sp],
    center_var = c("high", "low"),
    spread_var = c("high", "low"),
    shift_case_num = 7:9,
    predictor = c("none"),
    max_shift = c(24, 6)
)
synthetic_phenophases_shift_uncorr$kde <- apply(
    synthetic_phenophases_shift_uncorr,
    1,
    function(x, x_trial_sp = synthetic_phenophases) {
        single_sample <- subset(
            x_trial_sp,
            subset = (
                trial == as.integer(
                    x["trial"]
                ) & sp == x["sp"][1] & center_var == x[
                    "center_var"
                ][1] & spread_var == x[
                    "spread_var"
                ][1]
            )
        )
        p <- NULL
        if (x["predictor"][1] == "median") {
            p <- single_sample$center
        }
        if (x["predictor"][1] == "duration") {
            p <- single_sample$spread * 6
        }
        return (fshift_x(
            single_sample$kde[[1]],
            # find how much to shift by
            shift_phenophase_by(
                -1, # param that doesn't actually need to do anything here
                x["predictor"][1],
                as.integer(x["shift_case_num"]),
                max_shift = as.integer(x["max_shift"])
            )
        ))
    }
)

synthetic_phenophases_shift <- rbind(
    synthetic_phenophases_shift_corr,
    synthetic_phenophases_shift_uncorr
)

# Interaction potentials after simulated phenological shifts -------------------
# note: this section can take several hours

# empirical communities
calling_ip_shift <- map(
    group_split(group_by(
        calling_phenophases_shift,
        trial,
        shift_case_num,
        predictor,
        max_shift
    )),
    function(group) {
        ip_trial <- interaction_potential(group$kde)
        ip_trial$trial <- group$trial[1]
        ip_trial$shift_case_num <- group$shift_case_num[1]
        ip_trial$predictor <- group$predictor[1]
        ip_trial$max_shift <- group$max_shift[1]
        return (ip_trial)
    }
)
calling_ip_shift <- bind_rows(calling_ip_shift)
calling_ip_shift$var <- "high"

# synthetic communities
synthetic_ip_shift <- map(
    group_split(group_by(
        synthetic_phenophases_shift,
        trial,
        center_var,
        spread_var,
        shift_case_num,
        predictor,
        max_shift
    )),
    function(group) {
        ip_trial <- interaction_potential(group$kde)
        ip_trial$trial <- group$trial[1]
        ip_trial$center_var <- group$center_var[1]
        ip_trial$spread_var <- group$spread_var[1]
        ip_trial$shift_case_num <- group$shift_case_num[1]
        ip_trial$predictor <- group$predictor[1]
        ip_trial$max_shift <- group$max_shift[1]
        return (ip_trial)
    }
)
synthetic_ip_shift <- bind_rows(synthetic_ip_shift)

# Compute change in interaction potentials -------------------------------------

# Empirical communities

calling_ip <- merge(
    calling_ip_no_shift[, c("trial", "interval", "ip")],
    calling_ip_shift[, c(
        "trial",
        "interval",
        "shift_case_num",
        "predictor",
        "max_shift",
        "ip"
    )],
    by = c("trial", "interval")
)
calling_ip <- dplyr::rename(
    calling_ip,
    ip_no_shift = ip.x,
    ip_shift = ip.y
)
calling_ip$shift_case_num <- as.factor(calling_ip$shift_case_num)
calling_ip <- calling_ip[order(calling_ip$interval), ]
calling_ip$var <- "high"

calling_ip$ip_change <- calling_ip$ip_shift - calling_ip$ip_no_shift
calling_ip$direction_of_change <- ifelse(
    calling_ip$ip_change >= 0,
    "increase",
    "decrease"
)
calling_ip <- mutate(
    group_by(calling_ip, trial, shift_case_num, predictor, max_shift),
    ip_change_frac_init = ip_change / max(ip_no_shift)
)

# annual change
calling_ip_change_annual_direction <- summarise(
    group_by(
        calling_ip,
        trial,
        shift_case_num,
        predictor,
        max_shift,
        direction_of_change
    ),
    total = sum(ip_change),
    frac_init_total = sum(ip_change_frac_init)
)
calling_ip_change_annual_direction$var <- "high"
calling_ip_change_annual_direction_xy <- merge(
    x = rename(
        subset(
            calling_ip_change_annual_direction,
            direction_of_change == "decrease",
            select = -direction_of_change
        ),
        total_x = total,
        frac_init_total_x = frac_init_total
    ),
    y = rename(
        subset(
            calling_ip_change_annual_direction,
            direction_of_change == "increase",
            select = -direction_of_change
        ),
        total_y = total,
        frac_init_total_y = frac_init_total
    ),
    all = TRUE
)
calling_ip_change_annual_direction_xy_meanci95 <- summarise_at(
    group_by(
        calling_ip_change_annual_direction_xy,
        shift_case_num,
        predictor,
        max_shift
    ),
    c("total_x", "total_y", "frac_init_total_x", "frac_init_total_y"),
    list(
        mean = ~ mean(.),
        upper = ~ mean_cl_boot(.)$ymax,
        lower = ~ mean_cl_boot(.)$ymin
    )
)
calling_ip_change_annual_direction_xy_meanci95$var <- "high"

# Synthetic communities

synthetic_ip <- merge(
    synthetic_ip_no_shift[, c(
        "trial",
        "interval",
        "center_var",
        "spread_var",
        "ip"
    )],
    synthetic_ip_shift[, c(
        "trial",
        "interval",
        "center_var",
        "spread_var",
        "shift_case_num",
        "predictor",
        "max_shift",
        "ip"
    )],
    by = c("trial", "interval", "center_var", "spread_var")
)
synthetic_ip <- rename(
    synthetic_ip,
    ip_no_shift = ip.x,
    ip_shift = ip.y
)
synthetic_ip$shift_case_num <- as.factor(synthetic_ip$shift_case_num)
synthetic_ip <- synthetic_ip[order(synthetic_ip$interval), ]

synthetic_ip$ip_change <- synthetic_ip$ip_shift - synthetic_ip$ip_no_shift
synthetic_ip$direction_of_change <- ifelse(
    synthetic_ip$ip_change >= 0,
    "increase",
    "decrease"
)
synthetic_ip <- mutate(
    group_by(
        synthetic_ip,
        trial,
        shift_case_num,
        predictor,
        center_var,
        spread_var,
        max_shift
    ),
    ip_change_frac_init = ip_change / max(ip_no_shift)
)

# annual change
synthetic_ip_change_annual_direction <- summarise(
    group_by(
        synthetic_ip,
        trial,
        center_var,
        spread_var,
        shift_case_num,
        predictor,
        max_shift,
        direction_of_change
    ),
    total = sum(ip_change),
    frac_init_total = sum(ip_change_frac_init)
)
synthetic_ip_change_annual_direction_xy <- merge(
    x = rename(
        subset(
            synthetic_ip_change_annual_direction,
            direction_of_change == "decrease",
            select = -direction_of_change
        ),
        total_x = total,
        frac_init_total_x = frac_init_total
    ),
    y = rename(
        subset(
            synthetic_ip_change_annual_direction,
            direction_of_change == "increase",
            select = -direction_of_change
        ),
        total_y = total,
        frac_init_total_y = frac_init_total
    ),
    all = TRUE
)
synthetic_ip_change_annual_direction_xy_meanci95 <- summarise_at(
    group_by(
        synthetic_ip_change_annual_direction_xy,
        center_var,
        spread_var,
        shift_case_num,
        predictor,
        max_shift
    ),
    c("total_x", "total_y", "frac_init_total_x", "frac_init_total_y"),
    list(
        mean = ~ mean(.),
        upper = ~ mean_cl_boot(.)$ymax,
        lower = ~ mean_cl_boot(.)$ymin
    )
)

# store results of annual effect of variation in phenophase timing & duration
#  relative to high variation scenario
synthetic_ip_change_annual_direction_xy_meanci95$diff_from <- paste0(
    synthetic_ip_change_annual_direction_xy_meanci95$center_var,
    synthetic_ip_change_annual_direction_xy_meanci95$spread_var
)
synthetic_ip_loss_var_annual_direction_xy_meanci95 <- lapply(
    c("lowhigh", "highlow", "lowlow"), FUN = function(variation) {
        summarise(
            group_by(
                synthetic_ip_change_annual_direction_xy_meanci95,
                shift_case_num,
                predictor,
                max_shift
            ),
            "{variation}_minus_highhigh_total_x" := total_x_mean[
                match(variation, diff_from)
            ] - total_x_mean[match("highhigh", diff_from)],
            "frac_init_{variation}_minus_highhigh_total_x" := frac_init_total_x_mean[
                match(variation, diff_from)
            ] - frac_init_total_x_mean[match("highhigh", diff_from)],
            "{variation}_minus_highhigh_total_y" := total_y_mean[
                match(variation, diff_from)
            ] - total_y_mean[match("highhigh", diff_from)],
            "frac_init_{variation}_minus_highhigh_total_y" := frac_init_total_y_mean[
                match(variation, diff_from)
            ] - frac_init_total_y_mean[match("highhigh", diff_from)]
        )
    }
)
synthetic_ip_loss_var_annual_direction_xy_meanci95 <- merge(
    x = synthetic_ip_loss_var_annual_direction_xy_meanci95[[1]],
    y = merge(
        x = synthetic_ip_loss_var_annual_direction_xy_meanci95[[2]],
        y = synthetic_ip_loss_var_annual_direction_xy_meanci95[[3]],
        all = TRUE
    ),
    all = TRUE
)
synthetic_ip_loss_var_annual_direction_x_meanci95_melt <- merge(
    x = mutate(melt(
        select(
            synthetic_ip_loss_var_annual_direction_xy_meanci95,
            c(
                shift_case_num,
                predictor,
                max_shift,
                contains("total_x") & !(contains("frac_init"))
            )
        ),
        id.vars = c("shift_case_num", "predictor", "max_shift"),
        variable.name = "diff_from",
        value.name = "diff_from_highhigh_total_x"
    ), diff_from = sub("(\\w*)_minus.*", "\\1", diff_from)),
    y = mutate(melt(
        select(
            synthetic_ip_loss_var_annual_direction_xy_meanci95,
            c(
                shift_case_num,
                predictor,
                max_shift,
                contains("total_x") & contains("frac_init")
            )
        ),
        id.vars = c("shift_case_num", "predictor", "max_shift"),
        variable.name = "diff_from",
        value.name = "frac_init_diff_from_highhigh_total_x"
    ), diff_from = sub("frac_init_(\\w*)_minus.*", "\\1", diff_from)),
    all = TRUE
)
synthetic_ip_loss_var_annual_direction_y_meanci95_melt <- merge(
    x = mutate(melt(
        select(
            synthetic_ip_loss_var_annual_direction_xy_meanci95,
            c(
                shift_case_num,
                predictor,
                max_shift,
                contains("total_y") & !(contains("frac_init"))
            )
        ),
        id.vars = c("shift_case_num", "predictor", "max_shift"),
        variable.name = "diff_from",
        value.name = "diff_from_highhigh_total_y"
    ), diff_from = sub("(\\w*)_minus.*", "\\1", diff_from)),
    y = mutate(melt(
        select(
            synthetic_ip_loss_var_annual_direction_xy_meanci95,
            c(
                shift_case_num,
                predictor,
                max_shift,
                contains("total_y") & contains("frac_init")
            )
        ),
        id.vars = c("shift_case_num", "predictor", "max_shift"),
        variable.name = "diff_from",
        value.name = "frac_init_diff_from_highhigh_total_y"
    ), diff_from = sub("frac_init_(\\w*)_minus.*", "\\1", diff_from)),
    all = TRUE
)
synthetic_ip_loss_var_annual_direction_xy_meanci95_melt <- merge(
    x = synthetic_ip_loss_var_annual_direction_x_meanci95_melt,
    y = synthetic_ip_loss_var_annual_direction_y_meanci95_melt,
    all = TRUE
)

# Plot change in interaction potentials (Fig. 3-5) -----------------------------
# Figure 3a
calling_ip_change_annual_direction_xy_meanci95$shift_case_num <- as.integer(
    calling_ip_change_annual_direction_xy_meanci95$shift_case_num
)
png(
    filename = file.path(DATA_OUT, "figure_3a.png"),
    width = 18,
    height = 7.5,
    units = "cm",
    res = 600,
    bg = "transparent"
)
ggplot(
    subset(
        calling_ip_change_annual_direction_xy_meanci95,
        max_shift == 24 & shift_case_num < 7
    ),
    aes(
        x = abs(frac_init_total_x_mean),
        y = abs(frac_init_total_y_mean),
        shape = factor(shift_case_num)
    )
) +
    facet_wrap(
        facets = vars(predictor),
        labeller = labeller(predictor = c(
            `median` = "shift as function of median day",
            `duration` = "shift as function of duration"
        ))
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.33) +
    geom_point(size = 1.75, show.legend = TRUE) +
    geom_errorbar(aes(
        ymin = abs(frac_init_total_y_lower),
        ymax = abs(frac_init_total_y_upper)
    ), linewidth = 0.33, width = 0) +
    geom_errorbarh(aes(
        xmin = abs(frac_init_total_x_lower),
        xmax = abs(frac_init_total_x_upper)
    ), linewidth = 0.33, height = 0) +
    geom_point(data = subset(
        calling_ip_change_annual_direction_xy_meanci95,
        max_shift == 24 & shift_case_num >= 7,
        select = -predictor
    ), size = 2, show.legend = TRUE) +
    geom_errorbar(data = subset(
        calling_ip_change_annual_direction_xy_meanci95,
        max_shift == 24 & shift_case_num >= 7,
        select = -predictor
    ), aes(
        ymin = abs(frac_init_total_y_lower),
        ymax = abs(frac_init_total_y_upper)
    ), linewidth = 0.33, width = 0) +
    geom_errorbarh(data = subset(
        calling_ip_change_annual_direction_xy_meanci95,
        max_shift == 24 & shift_case_num >= 7,
        select = -predictor
    ), aes(
        xmin = abs(frac_init_total_x_lower),
        xmax = abs(frac_init_total_x_upper)
    ), linewidth = 0.33, height = 0) +
    annotate(
        "point",
        x = 0,
        y = 0,
        size = 2,
        colour = "black",
        fill = "gray",
        shape = 21
    ) +
    annotate(
        "text",
        x = 0,
        y = 0,
        label = "no change",
        colour = "black",
        vjust = -1,
        hjust = 0,
        size = 3.25
    ) +
    coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 4.5)) +
    scale_shape_manual(
        drop = FALSE,
        labels = c(
            "1: slope > 0; shift later",
            "2: slope > 0; shift earlier",
            "3: slope > 0; some shift later,\nsome shift earlier",
            "4: slope < 0; shift later",
            "5: slope < 0; shift earlier",
            "6: slope < 0; some shift later,\nsome shift earlier",
            "7: uncorrelated; shift later",
            "8: uncorrelated; shift earlier",
            "9: uncorrelated; some shift later,\nsome shift earlier"
        ),
        limits = as.character(1:9),
        values = c(2, 1, 7, 5, 6, 13, 15, 16, 18)
    ) +
    labs(
        shape = "shape of shift function",
        x = "Cumulative annual decrease in interaction potential",
        y = "Cumulative annual increase in interaction potential"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(0.75)),
        axis.text = element_text(size = rel(0.75), colour = "black"),
        axis.line = element_line(linewidth = 0.33),
        strip.text.x = element_text(
            size = rel(0.9),
            margin = margin(0.1, 0, 0.1, 0, unit = "cm")
        ),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(
            size = rel(0.75),
            margin = margin(0, 0, 0, 0, unit = "cm")
        ),
        legend.text = element_text(size = rel(0.75)),
        legend.key.spacing.y = unit(-1, 'mm')
    ) +
    guides(shape = guide_legend(override.aes = list(size = 2.5)))
dev.off()

# Figure 3b
synthetic_ip_change_annual_direction_xy_meanci95$shift_case_num <- as.integer(
    synthetic_ip_change_annual_direction_xy_meanci95$shift_case_num
)
png(
    filename = file.path(DATA_OUT, "figure_3b.png"),
    width = 18,
    height = 7.5,
    units = "cm",
    res = 600,
    bg = "transparent"
)
ggplot(
    subset(
        synthetic_ip_change_annual_direction_xy_meanci95,
        center_var == "high" & spread_var == "high" & max_shift == 24 & shift_case_num < 7,
        select = -c(max_shift, center_var, spread_var)
    ),
    aes(
        x = abs(frac_init_total_x_mean),
        y = abs(frac_init_total_y_mean),
        shape = factor(shift_case_num)
    )
) +
    facet_wrap(
        facets = vars(predictor),
        labeller = labeller(predictor = c(
            `median` = "shift as function of median day",
            `duration` = "shift as function of duration"
        ))
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.33) +
    geom_point(size = 1.75, show.legend = TRUE) +
    geom_errorbar(aes(
        ymin = abs(frac_init_total_y_lower),
        ymax = abs(frac_init_total_y_upper)
    ), linewidth = 0.33, width = 0) +
    geom_errorbarh(aes(
        xmin = abs(frac_init_total_x_lower),
        xmax = abs(frac_init_total_x_upper)
    ), linewidth = 0.33, height = 0) +
    geom_point(data = subset(
        synthetic_ip_change_annual_direction_xy_meanci95,
        center_var == "high" & spread_var == "high" & max_shift == 24 & shift_case_num >= 7,
        select = -c(max_shift, center_var, spread_var, predictor)
    ), size = 2, show.legend = TRUE) +
    geom_errorbar(data = subset(
        synthetic_ip_change_annual_direction_xy_meanci95,
        center_var == "high" & spread_var == "high" & max_shift == 24 & shift_case_num >= 7,
        select = -c(max_shift, center_var, spread_var, predictor)
    ), aes(
        ymin = abs(frac_init_total_y_lower),
        ymax = abs(frac_init_total_y_upper),
    ), linewidth = 0.33, width = 0) +
    geom_errorbarh(data = subset(
        synthetic_ip_change_annual_direction_xy_meanci95,
        center_var == "high" & spread_var == "high" & max_shift == 24 & shift_case_num >= 7,
        select = -c(max_shift, center_var, spread_var, predictor)
    ), aes(
        xmin = abs(frac_init_total_x_lower),
        xmax = abs(frac_init_total_x_upper)
    ), linewidth = 0.33, height = 0) +
    annotate(
        "point",
        x = 0,
        y = 0,
        size = 2,
        colour = "black",
        fill = "gray",
        shape = 21
    ) +
    annotate(
        "text",
        x = 0,
        y = 0,
        label = "no change",
        colour = "black",
        vjust = -1,
        hjust = 0,
        size = 3.25
    ) +
    coord_cartesian(xlim = c(0, 4.5), ylim = c(0, 4.5)) +
    scale_shape_manual(
        labels = c(
            "1: slope > 0; shift later",
            "2: slope > 0; shift earlier",
            "3: slope > 0; some shift later,\nsome shift earlier",
            "4: slope < 0; shift later",
            "5: slope < 0; shift earlier",
            "6: slope < 0; some shift later,\nsome shift earlier",
            "7: uncorrelated; shift later",
            "8: uncorrelated; shift earlier",
            "9: uncorrelated; some shift later,\nsome shift earlier"
        ),
        limits = as.character(1:9),
        values = c(2, 1, 7, 5, 6, 13, 15, 16, 18)
    ) +
    labs(
        shape = "shape of shift function",
        x = "Cumulative annual decrease in interaction potential",
        y = "Cumulative annual increase in interaction potential"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.title = element_text(size = rel(0.75)),
        axis.text = element_text(size = rel(0.75), colour = "black"),
        axis.line = element_line(linewidth = 0.33),
        strip.text.x = element_text(
            size = rel(0.9),
            margin = margin(0.1, 0, 0.1, 0, unit = "cm")
        ),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(
            size = rel(0.75),
            margin = margin(0, 0, 0, 0, unit = "cm")
        ),
        legend.text = element_text(size = rel(0.75)),
        legend.key.spacing.y = unit(-1, 'mm')
    ) +
    guides(shape = guide_legend(override.aes = list(size = 2.5)))
dev.off()

# Figure 4a right
png(
    filename = file.path(DATA_OUT, "figure_4a_right.png"),
    width = 18,
    height = 10.5,
    units = "cm",
    res = 600,
    bg = "transparent",
)
ggplot(
    calling_ip[calling_ip$max_shift == 24, ],
    aes(
        x = interval * 5,
        y = ip_change_frac_init
    )
) +
    geom_hline(
        yintercept = 0,
        colour = "black",
        linetype = 2,
        linewidth = 0.5
    ) +
    geom_smooth(aes(colour = predictor), linewidth = 0.5) +
    facet_wrap(
        facets = vars(shift_case_num),
        labeller = labeller(shift_case_num = c(
            "1" = "1: slope > 0; shift later",
            "2" = "2: slope > 0; shift earlier",
            "3" = "3: slope > 0; some shift later,\nsome shift earlier",
            "4" = "4: slope < 0; shift later",
            "5" = "5: slope < 0; shift earlier",
            "6" = "6: slope < 0; some shift later,\nsome shift earlier",
            "7" = "7: uncorrelated; shift later",
            "8" = "8: uncorrelated; shift earlier",
            "9" = "9: uncorrelated; some shift later,\nsome shift earlier"
        ))
    ) +
    # divide x-axis into thirds for tick labels
    scale_x_continuous(
        breaks = c(1, 122, 244, 366),
        labels = as.character(c(1, 122, 244, 366)),
        limits = c(1, 366)
    ) +
    # colorblind-friendly
    scale_color_manual(
        labels = c("median day", "duration", "none"),
        values = unname(palette.colors(n = 5))[2:4]
    ) +
    guides(colour = guide_legend(override.aes = list(linewidth = 0.75))) +
    labs(
        colour = "predictor of \nphenological shift",
        x = "Day of year (community adjusted)",
        y = "Change in interaction potential\nout of max initial (decrease < 0, increase > 0)",
    ) +
    geom_rect(
        data = calling_annotate_box,
        mapping = aes(
            xmin = q75_min_interval * 5,
            xmax = q75_max_interval * 5,
            ymin = -Inf,
            ymax = Inf,
            x = NULL,
            y = NULL
        ),
        alpha = 0.25,
        fill = "gray"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line = element_line(linewidth = 0.5),
        strip.text.x = element_text(
            size = rel(0.85),
            margin = margin(t = 0.1, r = 0, b = 0.1, l = 0, unit = "cm")
        ),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(
            size = rel(0.8),
            margin = margin(t = 0, r = 0, b = 0.1, l = 0, unit = "cm")
        ),
        legend.key.spacing.y = unit(1, 'mm'),
        legend.key.height = unit(4, 'mm')
    )
dev.off()

# Figure 4b right
png(
    filename = file.path(DATA_OUT, "figure_4b_right.png"),
    width = 18,
    height = 10.5,
    units = "cm",
    res = 600,
    bg = "transparent",
)
ggplot(
    synthetic_ip[
        synthetic_ip$center_var == "high" & synthetic_ip$spread_var == "high" & synthetic_ip$max_shift == 24,
    ],
    aes(
        x = interval * 5,
        y = ip_change_frac_init
    )
) +
    geom_hline(
        yintercept = 0,
        colour = "black",
        linetype = 2,
        linewidth = 0.5
    ) +
    geom_smooth(aes(colour = predictor), linewidth = 0.5) +
    facet_wrap(
        facets = vars(shift_case_num),
        labeller = labeller(shift_case_num = c(
            "1" = "1: slope > 0; shift later",
            "2" = "2: slope > 0; shift earlier",
            "3" = "3: slope > 0; some shift later,\nsome shift earlier",
            "4" = "4: slope < 0; shift later",
            "5" = "5: slope < 0; shift earlier",
            "6" = "6: slope < 0; some shift later,\nsome shift earlier",
            "7" = "7: uncorrelated; shift later",
            "8" = "8: uncorrelated; shift earlier",
            "9" = "9: uncorrelated; some shift later,\nsome shift earlier"
        ))
    ) +
    labs(
        colour = "predictor of \nphenological shift",
        x = "Day of year",
        y = "Change in interaction potential\nout of max initial (decrease < 0, increase > 0)"
    ) +
    # divide x-axis into thirds for tick labels
    scale_x_continuous(
        breaks = c(1, 122, 244, 366),
        labels = as.character(c(1, 122, 244, 366)),
        limits = c(1, 366)
    ) +
    # colorblind-friendly
    scale_color_manual(
        labels = c("median day", "duration", "none"),
        values = unname(palette.colors(n = 5))[2:4]
    ) +
    guides(colour = guide_legend(override.aes = list(linewidth = 0.75))) +
    geom_rect(
        data = subset(
            synthetic_annotate_box,
            center_var == "high" & spread_var == "high"
        ),
        mapping = aes(
            xmin = q75_min_interval * 5,
            xmax = q75_max_interval * 5,
            ymin = -Inf,
            ymax = Inf,
            x = NULL,
            y = NULL
        ),
        alpha = 0.25,
        fill = "gray"
    ) +
    mytheme +
    theme(
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line = element_line(linewidth = 0.5),
        strip.text.x = element_text(
            size = rel(0.85),
            margin = margin(t = 0.1, r = 0, b = 0.1, l = 0, unit = "cm")
        ),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.title = element_text(
            size = rel(0.8),
            margin = margin(t = 0, r = 0, b = 0.1, l = 0, unit = "cm")
        ),
        legend.key.spacing.y = unit(1, 'mm'),
        legend.key.height = unit(4, 'mm')
    )
dev.off()

# Figure 5
# package interference when displaying axis labels, load only ggplot2 package
#  and synthetic_phenology_sims_out.RData
png(
    filename = file.path(DATA_OUT, "figure_5.png"),
    width = 5000,
    height = 1800,
    res = 300
)
ggplot(
    subset(
        synthetic_ip_loss_var_annual_direction_xy_meanci95_melt,
        max_shift == 24
    ),
    aes(
        x = frac_init_diff_from_highhigh_total_x,
        y = frac_init_diff_from_highhigh_total_y,
        colour = predictor,
        shape = factor(shift_case_num)
    )
) +
    facet_wrap(
        facets = vars(diff_from),
        labeller = labeller(diff_from = c(
            "lowhigh" = "loss of variation in\nmedian day",
            "highlow" = "loss of variation in\nduration",
            "lowlow" = "loss of variation in\nmedian day & duration"
        ))
    ) +
    geom_abline(
        slope = -1,
        intercept = 0,
        linetype = 2,
        alpha = 0.25
    ) +
    geom_vline(
        xintercept = 0,
        colour = "black",
        linetype = 1,
        alpha = 0.5
    ) +
    geom_hline(
        yintercept = 0,
        colour = "black",
        linetype = 1,
        alpha = 0.5
    ) +
    geom_point(
        size = 5,
        alpha = 0.825
    ) +
    # colorblind-friendly
    scale_color_manual(
        labels = c("median day", "duration", "none"),
        values = unname(palette.colors(n = 5))[2:4]
    ) +
    scale_shape_manual(
        labels = c(
            "1: slope > 0; shift later",
            "2: slope > 0; shift earlier",
            "3: slope > 0; some shift later,\nsome shift earlier",
            "4: slope < 0; shift later",
            "5: slope < 0; shift earlier",
            "6: slope < 0; some shift later,\nsome shift earlier",
            "7: uncorrelated; shift later",
            "8: uncorrelated; shift earlier",
            "9: uncorrelated; some shift later,\nsome shift earlier"
        ),
        limits = as.character(1:9),
        values = c(2, 1, 7, 5, 6, 13, 15, 16, 18)
    ) +
    labs(
        colour = "predictor of phenological shift",
        shape = "shape of shift function",
        x = "\u2190 decreases amplified          decreases dampened \u2192\n\neffect on mean cumulative annual decrease in interaction potential",
        y = "effect on mean cumulative annual increase\nin interaction potential\n\n\u2190 increases dampened     increases amplified \u2192"
    ) +
    coord_cartesian(xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5)) +
    mytheme +
    theme(
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5), colour = "black"),
        panel.spacing.x = unit(1.5, "lines"),
        strip.text.x = element_text(size = rel(1.75)),
        legend.margin = margin(b = 5, l = 8, unit = "mm"),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        legend.key.spacing.y = unit(1.5, 'mm')
    )
dev.off()

# Save to .RData files (as needed) ---------------------------------------------
save(
    base_pdfs,
    callperiods_frog,
    p_frog,
    p_frog_extremes50,
    p_frog_mid50,
    p_frog_mid96,
    file = file.path(DATA_OUT_RDATA, "anuran_phenology.RData")
)
save(
    calling_phenophases,
    synthetic_phenophases,
    file = file.path(DATA_OUT_RDATA, "phenology_sims_samples.RData")
)
save(
    calling_phenophases_var_stats,
    synthetic_phenophases_var_stats,
    file = file.path(DATA_OUT_RDATA, "phenology_sims_samples_variance.RData")
)
save(
    calling_phenophases_no_shift,
    calling_phenophases_shift,
    calling_annotate_box,
    calling_ip_no_shift,
    calling_ip_shift,
    calling_ip,
    calling_ip_change_annual_direction,
    calling_ip_change_annual_direction_xy,
    calling_ip_change_annual_direction_xy_meanci95,
    file = file.path(DATA_OUT_RDATA, "anuran_phenology_sims.RData")
)
save(
    synthetic_annotate_box,
    synthetic_phenophases_shift,
    synthetic_ip_no_shift,
    synthetic_ip_shift,
    synthetic_ip,
    synthetic_ip_change_annual_direction,
    synthetic_ip_change_annual_direction_xy,
    synthetic_ip_change_annual_direction_xy_meanci95,
    synthetic_ip_loss_var_annual_direction_xy_meanci95_melt,
    file = file.path(DATA_OUT_RDATA, "synthetic_phenology_sims.RData")
)
