
### This script pulls in PHAVR linelist file data, does a little cleaning, and creates variables

# set-up ------------------------------------------------------------------

# set working directory - change path to your own folder

if (file.exists(Sys.getenv('MCW_EPI_WD'))){
  mywd <- Sys.getenv('MCW_EPI_WD')
}

if (!exists("mywd")){
  mywd <- "C:/Users/lrein/mcw.edu/COVID19 Spatiotemporal Epidemiology - General"
}

### load R packages
packs <- c("tidyverse"
           ,"knitr"
           ,"devtools"
           ,"reshape2"
           ,"readxl"
           ,"plyr"
           ,"dplyr"
           ,"gridExtra"
           ,"lubridate"
           ,"scales"
           ,"Hmisc"
           ,"nls2"
           ,"ggrepel"
           ,"zoo"
           ,"incidence"
           ,"EpiEstim"
           ,"rgdal"
           ,"DescTools"
           ,"htmlTable"
           )
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = packs)

ls_save <- ls()

# select and import PHAVR data --------------------------------------------

# data directory
if (!exists("datdir")) datdir <- paste(mywd, "/../Beyer, Kirsten - WEDSS_Archive", sep="")

# get all WEDSS datasets from data directory, following specific file naming convention
# allfnames <- list.files(datdir, pattern= "^MKECovid19LineList_\\d{1,2}\\.\\d{1,2}\\.\\d{1,2}")
allfnames <- list.files(datdir, pattern= "^WedssCovid19MKESELineList\\d{1,2}\\.\\d{1,2}\\.\\d{1,2}")
allfnames <- allfnames[grepl(".csv$", allfnames)]

# grab extract dates from file names
# allfdates <- unlist(lapply(allfnames, function(x) strsplit(x, "MKECovid19LineList_")[[1]][2]))
allfdates <- unlist(lapply(allfnames, function(x) strsplit(x, "WedssCovid19MKESELineList")[[1]][2]))
allfdates <- gsub("\\.csv", "", allfdates)
allfdates <- unlist(lapply(allfdates, function(x) strsplit(x, "_")[[1]][1]))
allfdates <- as.Date(as.character(allfdates), format = "%m.%d.%y")

allphavr <- data.frame("fname" = allfnames
                       ,"date"  = allfdates
                       )

# grab most recent available date
if (!exists("dat_date")) dat_date <- max(allfdates)
dat_date <- as.Date(dat_date)

### only consider most recent pull
phavr_new <- subset(allphavr, date == dat_date)

# if multiple pulls were made in same day, use most recent:
if (nrow(phavr_new) > 1){
  phavr_new <- phavr_new[order(phavr_new$fname),]
  phavr_new <- tail(phavr_new, 1)
}

# grab time of last data pull
lasttime <- strsplit(as.character(phavr_new$fname), "_")[[1]][2]
lasttime <- gsub(".csv", "", lasttime)
lasttime <- gsub("\\.", ":", lasttime)

### read in data
wdat <- data.frame(read_csv(paste(datdir, phavr_new$fname, sep="/")
                              ,guess_max = 50000
                              )
                   )
# make variable names a little nicer
names(wdat) <- tolower(names(wdat))

# set last day to include in analysis, the most recent full day of data by default
if (!exists("lastdate")) lastdate <- max(phavr_new$date) - 1

# do not include anything before March 1st - there are some typos with dates from prior to 2020
alldates <- seq(as.Date("2020-03-01"), lastdate, 1)

# dates -------------------------------------------------------------------

# format dates
checkdate <- function(oldvar, newvar, n = 20, all = FALSE, dat = wdat){
  temp <- dat[!is.na(dat[,oldvar]),]
  if (all) n <- nrow(temp)
  temp[sample(1:nrow(temp), size = n), c(oldvar, newvar)]
}

wdat$incdate <- as.Date(wdat$specimendatereceived, "%m/%d/%Y")
checkdate("specimendatereceived", "incdate", 20)

wdat$coldate <- as.Date(wdat$specimendatecollected, "%m/%d/%Y")
checkdate("specimendatecollected", "coldate", 20)

wdat$epidate <- as.Date(wdat$episodedate, "%m/%d/%Y")
checkdate("episodedate", "epidate", 20)

wdat$admitdate <- as.Date(wdat$dateadmitted, "%m/%d/%Y")
checkdate("dateadmitted", "admitdate", 20)

wdat$deathdate <- as.Date(wdat$dateofdeath, "%m/%d/%Y")
checkdate("dateofdeath", "deathdate", 20)

wdat$dischargedate <- as.Date(wdat$datedischarged, "%m/%d/%Y")
checkdate("datedischarged", "dischargedate", 20)

wdat$credate <- as.Date(wdat$createdate, "%m/%d/%Y")
checkdate("createdate", "credate", 20)

wdat$onset_date <- wdat$onsetdate
wdat$onsetdate <- as.Date(wdat$onset_date, "%m/%d/%Y")
checkdate("onset_date", "onsetdate", 20)

### Specimen collection date is used to date all cases in the analysis
### if missing impute first with episode date, then record creation date

## if collection date is missing, impute with episode date:
table(is.na(wdat$coldate), wdat$resolutionstatus)
wdat$coldate[is.na(wdat$coldate)] <- wdat$epidate[is.na(wdat$coldate)]

## if collection date still missing, use createdate
table(is.na(wdat$coldate), wdat$resolutionstatus)
wdat$coldate[is.na(wdat$coldate)] <- wdat$credate[is.na(wdat$coldate)]

table(is.na(wdat$coldate), wdat$resolutionstatus)

# create event indicator variables ----------------------------------------

wdat$resolution_status <- wdat$resolutionstatus
wdat$patient_hospitalized <- wdat$patienthospitalized
wdat$patient_died_of_this_illness <- wdat$patientdiedofthisillness

# testing indicator
wdat$test01 <- 1

# case indicator
wdat$case01 <- 0
wdat$case01[wdat$resolutionstatus %in% "Confirmed"] <- 1

# hospitalized case indicator
wdat$hosp01 <- 0
wdat$hosp01[wdat$resolutionstatus %in% "Confirmed" & wdat$patienthospitalized %in% "Y"] <- 1

# deceased case indicator
wdat$died01 <- 0
wdat$died01[wdat$resolutionstatus %in% "Confirmed" & wdat$patientdiedofthisillness %in% "Y"] <- 1

wdat$case01_yn <- factor(wdat$case01, levels = c(0, 1), labels = c("No", "Yes"))
wdat$hosp01_yn <- factor(wdat$hosp01, levels = c(0, 1), labels = c("No", "Yes"))
wdat$died01_yn <- factor(wdat$died01, levels = c(0, 1), labels = c("No", "Yes"))

# review any deaths that are missing a deathdate
checkdeaths <- wdat[wdat$died01 == 1
                    ,c("clientid"
                       ,"sex"
                       ,"age"
                       ,"race"
                       ,"dateofdeath"
                       ,"deathdate"
                       ,"episodedate"
                       ,"countyfips"
                       ,"jurisdiction"
                    )]
table(checkdeaths$deathdate %in% alldates)
checkdeaths[!(checkdeaths$deathdate %in% alldates),]

# check for missing county_fips -------------------------------------------

### The cleaned 'county_fips' variable is used to subset data for milwaukee co and regional reports
### if countyfips is missing, impute with 'countyofresidence' variable 

wdat$county_fips <- wdat$countyfips
wdat$census_tract <- wdat$censustract
wdat$census_block <- wdat$censusblock

wdat$county_fips[wdat$county_fips %in% "N/A"] <- NA

## imputation by countyofresidence
wdat$county_fips[is.na(wdat$county_fips) & (wdat$countyofresidence %in% "Milwaukee" )] <- 55079
wdat$county_fips[is.na(wdat$county_fips) & (wdat$countyofresidence %in% "Ozaukee"   )] <- 55089
wdat$county_fips[is.na(wdat$county_fips) & (wdat$countyofresidence %in% "Washington")] <- 55131
wdat$county_fips[is.na(wdat$county_fips) & (wdat$countyofresidence %in% "Kenosha"   )] <- 55059
wdat$county_fips[is.na(wdat$county_fips) & (wdat$countyofresidence %in% "Walworth"  )] <- 55127
wdat$county_fips[is.na(wdat$county_fips) & (wdat$countyofresidence %in% "Racine"    )] <- 55101
wdat$county_fips[is.na(wdat$county_fips) & (wdat$countyofresidence %in% "Waukesha"  )] <- 55133

table(is.na(wdat$countyfips), wdat$resolutionstatus, useNA='ifany')
table(is.na(wdat$county_fips), wdat$resolutionstatus, useNA='ifany')

table(wdat$jurisdiction
      , is.na(wdat$county_fips)
      , wdat$resolutionstatus %in% "Confirmed" 
      , useNA='ifany'
      )

# create analytic variables -----------------------------------------------

### clean location variables - add leading zeros if needed
wdat$county_fips <- str_pad(wdat$county_fips, width=5, side="left", pad="0")
wdat$census_tract <- str_pad(wdat$census_tract, width=6, side="left", pad="0")
wdat$census_block <- str_pad(wdat$census_block, width=4, side="left", pad="0")

wdat$blockid <- paste(wdat$county_fips
                      ,wdat$census_tract
                      ,wdat$census_block
                      ,sep = ""
                      )
wdat$block10 <- substring(wdat$blockid, 1, 12)

capstr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

wdat$city <- unlist(lapply(wdat$city, function(x) capstr(tolower(x))))
table(wdat$city)
wdat$city[wdat$city %in% c("Milwuakee"
                           , "Milwaukeewi"
                           , "Milwaukke  Wi"
                           , "Milwa"
                           , "Mil"
                           , "Milw."
                           , "Milw"
                           , "Miilwaukee"
                           , "Milwa"
                           , "Milwauke"
                           , "Mke"
                           , "Miwaukee"
                           , "Miwlaukee"
                           , "Milwakee"
                           , "Mlwaukee"
                           , "Milwaukeer"
                           , "Milwaikee"
                           , "Milwauee"
                           , "Milwaueke"
                           , "Milwaauujkeee"
                           , "Mikwaukee"
                           )] <- "Milwaukee"
wdat$city[wdat$city %in% c("W Bend", "W. Bend")] <- "West Bend"
wdat$city[wdat$city %in% c("W Milwaukee", "W. Milwaukee")] <- "West Milwaukee"
wdat$city[wdat$city %in% c("S Milwaukee", "S. Milwaukee")] <- "South Milwaukee"
wdat$city[wdat$city %in% c("W Allis"
                           , "West Alis"
                           , "West  Allis"
                           , "West Allis."
                           )] <- "West Allis"
wdat$city[wdat$city %in% c("St Francis", "St. Francis")] <- "Saint Francis"
wdat$city[wdat$city %in% c("St Thomas", "St. Thomas")] <- "Saint Thomas"
wdat$city[wdat$city %in% c("St Charles", "St. Charles")] <- "Saint Charles"
wdat$city[wdat$city %in% c("Brook Field")] <- "Brookfield"
wdat$city[wdat$city %in% c("Germantown, Wi")] <- "Germantown"
wdat$city[wdat$city %in% c("Hales Corner")] <- "Hales Corners"
wdat$city[wdat$city %in% c("Hartland Wi", "Heartland Waukesh")] <- "Heartland"
wdat$city[wdat$city %in% c("Mt Pleasant")] <- "Mount Pleasant"
wdat$city[wdat$city %in% c("Menomonee Fls", "Menomonee", "Monaminee")] <- "Menomonee Falls"
wdat$city[wdat$city %in% c("Mequan", "Mequeon")] <- "Mequon" 
wdat$city[wdat$city %in% c("Keosha")] <- "Kenosha"
wdat$city[wdat$city %in% c("Wauwakeshaw", "Waukeshaw")] <- "Waukesha"
wdat$city[wdat$city %in% c("Wauwarosa"
                           , "Wauwatoa"
                           , "Wauwatosa"
                           , "Wauwatosaq"
                           , "Wauwautosa"
                           , "Wauwtosa"
                           , "Wawatosa"
                           , "Wauwautosa"
                           )] <- "Wauwatosa"
wdat$city[wdat$city %in% c("."
                           ,"/xysgt"
                           ,"0"
                           ,"Unknown"
                           ,"Not Provided"
                           ,"Unk"
                           ,"Unknoqn"
                           ,"Unkown"
                           ,"Uknown"
                           ,"Na"
                           ,"Nana"
                           ,NA
                           )] <- "Unknown"
table(wdat$city)

# combined race ethnicity variable
wdat$racethnicity <- wdat$race
wdat$racethnicity[wdat$racethnicity %in% "Black or African American"] <- "Black or AA"
wdat$racethnicity[wdat$racethnicity %in% "American Indian or Alaska Native"] <- "AIAN"
wdat$racethnicity[wdat$racethnicity %in% "Native Hawaiian or Other Pacific Islander"] <- "NHOPI"
wdat$racethnicity[wdat$racethnicity %in% c("Unknown", NA)] <- "Unknown"
wdat$racethnicity[wdat$ethnicity %in% c("Hispanic or Latino")] <- "Hispanic"
wdat$racethnicity[!is.na(wdat$racethnicity) & 
                    !(wdat$racethnicity %in% c("Black or AA"
                                               , "White"
                                               , "AIAN"
                                               , "NHOPI"
                                               , "Hispanic"
                                               , "Unknown"
                                               , "Multiple Races"
                                               , "Asian"
                                               )
                      )] <- "Other"
table("Pre" = wdat$race, "Post" = wdat$racethnicity, useNA='always')

table(wdat$sex, useNA='ifany')
wdat$gender <- wdat$sex
wdat$gender[is.na(wdat$gender) | !(wdat$gender %in% c("Male", "Female"))] <- "Other/Unknown"
table(wdat$sex, wdat$gender)

# age categories
table(wdat$agecat2 <- cut(wdat$age
                          , breaks = c(0, 18, 40, 60, 80, Inf)
                          , labels = c("< 18"
                                       , "18-39"
                                       , "40-59"
                                       , "60-79"
                                       , "80+"
                          )
                          , include.lowest = TRUE
                          , right = FALSE
)
, wdat$resolutionstatus
, useNA='ifany'
)
table(wdat$age, wdat$agecat2, useNA='ifany')

table(wdat$agecat6 <- cut(wdat$age
                          , breaks = c(0, 18, 30, 40, 60, 80, Inf)
                          , labels = c("< 18"
                                       , "18-29"
                                       , "30-39"
                                       , "40-59"
                                       , "60-79"
                                       , "80+"
                          )
                          , include.lowest = TRUE
                          , right = FALSE
)
, wdat$resolutionstatus
, useNA='ifany'
)
table(wdat$age, wdat$agecat6, useNA='ifany')

# inclusion/exclusion -----------------------------------------------------

# total number of rows:
n0 <- nrow(wdat)

# number of rows with status "Not A Case" or "Confirmed":
# wdat <- subset(wdat, resolutionstatus %in% c("Not A Case", "Confirmed"))
# n2_res <- nrow(wdat)

# number of rows with non-missing collection or episode date:
wdat <- subset(wdat, !is.na(coldate))
n2_cdate <- nrow(wdat)

# Filter out Contact Investigations per Chris Steward @ WI DHS
table(wdat$recordtype)
wdat <- subset(wdat, recordtype %in% "Disease Incident")
n3_rtype <- nrow(wdat)

# review testing types ---------------------------------------------------------

# check if we are capturing any antibody results?
wdat <- wdat[order(wdat$coldate),]

table(wdat$specimenresult)

wdat$specimenresult <- toupper(wdat$specimenresult)

abtxt <- paste(c("IGG", " AB ", "SEROLOGY", "ANTIBODY", "IGM", "SERA"), collapse="|")

table(subset(wdat, grepl(abtxt, specimenresult))$specimenresult)

table(subset(wdat, !(grepl(abtxt, specimenresult)))$specimenresult)

# Add block group variables  ----------------------------------------------

#### Merge in block group income variables from census data
bginc <- data.frame(read_csv(paste(mywd, "/Data Sources/WEDSS/census/ACS5y_SEShousing_mkebg2018_short.csv", sep="")))
names(bginc) <- tolower(names(bginc))

bginc$block10 <- as.character(bginc$id)

wdat2 <- merge(wdat, bginc, by = "block10", all.x = TRUE, all.y = FALSE)

bginc_pop <- rep(bginc$medhhinc, times = bginc$total)
quantile(bginc_pop, na.rm=T)
# 0%    25%    50%    75%   100% 
# 7872  35833  50096  68393 250001 

wdat2$medhhinc_q <- cut(wdat2$medhhinc
                        , breaks = c(-Inf, 35833, 50096, 68393, Inf)
                        , labels = c("Q1: Lowest"
                                     ,"Q2: Low-medium"
                                     ,"Q3: Medium-high"
                                     ,"Q4: Highest"
                                     )
                        )
table(wdat2$medhhinc_q, wdat2$county_fips %in% 55079, useNA='ifany')
wdat <- wdat2

table(wdat$county_fips %in% 55079, wdat$medhhinc_q, wdat$resolutionstatus %in% "Confirmed", useNA='ifany')

# Add Milwaukee Co denominators for pop based rates ----------------------------------------------

### Merge in denominators for population based rates

pop <- data.frame(read_csv(paste(mywd, "/Data Sources/census/ACS5y_county_wi_pop_groups2018.csv", sep="")))

popmke <- subset(pop, Countyfips %in% 55079)

# gender
popmke_gender <- data.frame("gender" = c("Male"
                                         , "Female"
                                         , "Other/Unknown"
                                         )
                            ,"pop_gender" = c(popmke$Male
                                              , popmke$Female
                                              , NA
                                              )
                            ,"county_fips" = 55079
                            )
sum(popmke_gender[,2], na.rm=T)

# age
popmke_agecat2 <- data.frame("agecat2" = c("< 18"
                                           , "18-39"
                                           , "40-59"
                                           , "60-79"
                                           , "80+"
                                           )
                             ,"pop_agecat2" = c(popmke$Tot18under
                                                    , popmke$Tot19_39
                                                    , popmke$Tot40_59
                                                    , popmke$Tot60_79
                                                    , popmke$Tot80over
                                                )
                             ,"county_fips" = 55079
                             )
sum(popmke_agecat2[,2], na.rm=T)

# race/eth
popmke_racethnicity <- data.frame("racethnicity" = c("Black or AA"
                                                     , "White"
                                                     , "Hispanic"
                                                     , "Asian"
                                                     , "AIAN"
                                                     , "NHOPI"
                                                     , "Multiple Races"
                                                     , "Other"
                                                     , "Unknown"
                                                     )
                                ,"pop_racethnicity" = c(popmke$NHBlack
                                                        ,popmke$NHWhite 
                                                        ,popmke$Hispanic
                                                        ,popmke$NHAsian
                                                        ,popmke$NHAIAN
                                                        ,popmke$NHNHPI
                                                        ,popmke$NHTwomore
                                                        ,popmke$NHOther
                                                        ,NA
                                                        )
                                ,"county_fips" = 55079
                                )
sum(popmke_racethnicity[,2], na.rm=T)

# income - not doing pop based rates...
incs <- names(popmke)[grepl("Inc",names(popmke))]
sum(popmke[,incs], na.rm=T)

wdat <- merge(wdat, popmke_gender      , by = c("county_fips", "gender"      ), all.x = TRUE, all.y = FALSE)
wdat <- merge(wdat, popmke_agecat2     , by = c("county_fips", "agecat2"     ), all.x = TRUE, all.y = FALSE)
wdat <- merge(wdat, popmke_racethnicity, by = c("county_fips", "racethnicity"), all.x = TRUE, all.y = FALSE)

# create variables to calculate rates per 1000 MKE Co. residents...
wdat$case01_pop_gender <- (wdat$case01/wdat$pop_gender)*1000
wdat$case01_pop_agecat2 <- (wdat$case01/wdat$pop_agecat2)*1000
wdat$case01_pop_racethnicity <- (wdat$case01/wdat$pop_racethnicity)*1000

wdat$hosp01_pop_gender <- (wdat$hosp01/wdat$pop_gender)*1000
wdat$hosp01_pop_agecat2 <- (wdat$hosp01/wdat$pop_agecat2)*1000
wdat$hosp01_pop_racethnicity <- (wdat$hosp01/wdat$pop_racethnicity)*1000

wdat$died01_pop_gender <- (wdat$died01/wdat$pop_gender)*1000
wdat$died01_pop_agecat2 <- (wdat$died01/wdat$pop_agecat2)*1000
wdat$died01_pop_racethnicity <- (wdat$died01/wdat$pop_racethnicity)*1000

# collapse wdat into a person-level dataset  ----------------------------------

idat <- wdat

### keep first confirmed incident if any, otherwise first not a case

wdat$rsorder <- 99
wdat$rsorder[wdat$resolution_status %in% "Confirmed" ] <- 1
wdat$rsorder[wdat$resolution_status %in% "Not A Case"] <- 2
wdat$rsorder[wdat$resolution_status %in% "Suspect"   ] <- 3
wdat$rsorder[wdat$resolution_status %in% "Probable"  ] <- 4

table(wdat$resolution_status)

# check how many people have more than 1 positive incident
# also check if anyone was positive more than once and had a hosp or death after the 1st incidence
wdat$deathdate[wdat$died01 == 0] <- NA

wdat2 <- wdat %>%
  dplyr::group_by(clientid) %>%
  dplyr::mutate(nincidentids = length(unique(incidentid))) %>%
  dplyr::mutate(nincidentids_pos = sum(case01)) %>%
  dplyr::mutate(died01_ever = max(died01)) %>%
  dplyr::mutate(hosp01_ever = max(hosp01)) %>%
  dplyr::mutate(deathdate_ever = min(deathdate))

wdat <- wdat2

wdat <- wdat[order(wdat$clientid, wdat$rsorder, wdat$coldate),]
wdat <- wdat[!duplicated(wdat$clientid),]

table(wdat$resolution_status)

table("First" = wdat$hosp01, "Ever" = wdat$hosp01_ever, wdat$resolution_status %in% "Confirmed")
table("First" = wdat$died01, "Ever" = wdat$died01_ever, wdat$resolution_status %in% "Confirmed")

# if someone had multiple positive incidents, make sure we capture whether that person was ever hospitalized or died after 1st dx
wdat$hosp01 <- wdat$hosp01_ever
wdat$died01 <- wdat$died01_ever
wdat$deathdate <- wdat$deathdate_ever

# remove unnecessary objects from workspace -------------------------------

wdat <- data.frame(wdat)

# save(wdat, file = paste(mywd, "/../Beyer, Kirsten - WEDSS_Archive/lrein_data/", dat_date, "_phavr_linelist.Rdata", sep=""))
# write.csv(wdat, file = paste(mywd, "/../Beyer, Kirsten - WEDSS_Archive/lrein_data/", dat_date, "_phavr_linelist.csv", sep=""))

rm(list=setdiff(ls(), c(ls_save, "wdat", "idat", "lastdate", "lasttime", "alldates")))







