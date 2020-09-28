
### This script pulls in PHAVR labs file, does a little cleaning, and creates variables

# set-up ------------------------------------------------------------------

# set working directory - may have to change path to your own folder

if(file.exists(Sys.getenv('MCW_EPI_WD'))){
  mywd <- Sys.getenv('MCW_EPI_WD')
  home <- Sys.getenv('HOME')
  datdir<- paste0(home,'/Documents/mcw.edu/Beyer, Kirsten - WEDSS_Archive/')
}

if (!exists("mywd")){
      mywd <- "C:/Users/lrein/mcw.edu/COVID19 Spatiotemporal Epidemiology - General"
}

### load R packages
packs <- c("plyr"
           ,"tidyverse"
           ,"knitr"
           ,"devtools"
           ,"reshape2"
           ,"readxl"
           ,"lubridate"
           ,"scales"
           ,"Hmisc"
           ,"zoo"
           ,"DescTools"
           ,"htmlTable"
)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(char = packs)

# save list of items in current workspace for later
ls_save <- ls()

# select and import PHAVR data --------------------------------------------

# data directory
if (!exists("datdir")) datdir <- paste(mywd, "/../Beyer, Kirsten - WEDSS_Archive", sep="")

# get all WEDSS labs datasets from data directory, following specific file naming convention
# allbnames <- list.files(datdir, pattern= "^WedssCovid19LabsMKECounty_\\d{1,2}\\.\\d{1,2}\\.\\d{1,2}")
allbnames <- list.files(datdir, pattern= "^WedssCovid19LabsMKESECounty\\d{1,2}\\.\\d{1,2}\\.\\d{1,2}")
allbnames <- allbnames[grepl(".csv$", allbnames)]

# grab extract dates from file names
# allbdates <- unlist(lapply(allbnames, function(x) strsplit(x, "WedssCovid19LabsMKECounty_")[[1]][2]))
allbdates <- unlist(lapply(allbnames, function(x) strsplit(x, "WedssCovid19LabsMKESECounty")[[1]][2]))
allbdates <- gsub("\\.csv", "", allbdates)
allbdates <- unlist(lapply(allbdates, function(x) strsplit(x, "_")[[1]][1]))
allbdates <- as.Date(as.character(allbdates), format = "%m.%d.%y")
alllab <- data.frame("fname" = allbnames
                    ,"date"  = allbdates
                    )

# grab most recent available date
if (!exists("dat_date")) dat_date <- max(allbdates)
if (dat_date > max(allbdates)) dat_date <- max(allbdates)
dat_date <- as.Date(dat_date)

### only consider most recent pull
lab_new   <- subset(alllab, date == dat_date)

# if multiple pulls were made in same day, use most recent:
if (nrow(lab_new) > 1){
  lab_new <- lab_new[order(lab_new$fname),]
  lab_new <- tail(lab_new, 1)
}

# read in data
bdat <- data.frame(read_csv(paste(datdir, lab_new$fname, sep="/")
                            ,guess_max = 25000
                            )
)

# make variable names a little nicer
names(bdat) <- tolower(names(bdat))

# collapse lab result data ------------------------------------------------

# by default, only clean data for milwaukee co. jurisdictions
if (!exists("jurs")){
  jurs <- c("Milwaukee"
            ,"Greenfield"
            ,"Cudahy" 
            ,"Oak Creek"
            ,"West Allis" 
            ,"Greendale"       
            ,"Wauwatosa"       
            ,"North Shore"     
            ,"Franklin"   
            ,"Saint Francis"
            ,"South Milwaukee"
            ,"Hales Corners" 
            
            # ,"Ozaukee County"
            # ,"Kenosha County"
            # ,"Waukesha County"
            # ,"Central Racine"
            # ,"Racine"
            # ,"Washington County"
            # ,"Walworth County"
            # 
            # ,"Waukesha County"
            )
}
bdat <- subset(bdat, jurisdiction %in% jurs)

# define regions to use for regional reports...
bdat$regions <- NA
bdat$regions[bdat$jurisdiction %in% c("Milwaukee"
                                      ,"Greenfield"
                                      ,"Cudahy" 
                                      ,"Oak Creek"
                                      ,"West Allis" 
                                      ,"Greendale"       
                                      ,"Wauwatosa"       
                                      ,"North Shore"     
                                      ,"Franklin"   
                                      ,"Saint Francis"
                                      ,"South Milwaukee"
                                      ,"Hales Corners" )] <- "Milwaukee"
bdat$regions[bdat$jurisdiction %in% c("Racine")] <- "City of Racine"
bdat$regions[bdat$jurisdiction %in% c("Kenosha County")] <- "Kenosha"
bdat$regions[bdat$jurisdiction %in% c("Walworth County")] <- "Walworth"
bdat$regions[bdat$jurisdiction %in% c("Washington County")] <- "Washington"
bdat$regions[bdat$jurisdiction %in% c("Ozaukee County")] <- "Ozaukee"  

## correct lab dataset dates:
bdat$coldate <- as.Date(bdat$speccollecteddate, "%m/%d/%Y")

bdat$incdate <- as.Date(bdat$resultdate, "%m/%d/%Y")

bdat$epidate <- as.Date(bdat$episodedate, "%m/%d/%Y")

bdat$credate <- as.Date(bdat$createdate, "%m/%d/%Y")

### Specimen collection date is used to date all cases in the analysis
### if missing impute first with episode date, then record creation date
bdat$coldate[is.na(bdat$coldate)] <- bdat$epidate[is.na(bdat$coldate)]
bdat$coldate[is.na(bdat$coldate)] <- bdat$credate[is.na(bdat$coldate)]

## look at staging status:
table(bdat$importstatus, useNA='always')
table(bdat$importstatus, bdat$result, useNA='always')

## look at recordtypes:
table(bdat$recordtype, useNA='always')
table(bdat$recordtype, bdat$result, useNA='always')

## remove any antibody test codes:
abcodes <- c('41460-7'
            ,'94505-5'
            ,'94507-1'
            ,'94562-6'
            ,'94563-4'
            ,'94564-2'
            ,'94762-2'
            ,'PLT2335'
            ,'94547-7'
            ,'94508-9'
            ,'64564-2'
            ,'94761-4'
            ,'94661-6'
            )
bdat$abtest <- 0
bdat$abtest[bdat$testcode %in% abcodes] <- 1

table(bdat$abtest, bdat$serology, useNA='always') # what is the 'serology' field?

bdat <- subset(bdat, !(testcode %in% abcodes))

### de-duplicate same collection day tests

# rows per incident:
table(table(bdat$incidentid))

bdat$facilityname <- toupper(bdat$facilityname)

bdatc <- bdat %>%
  dplyr::group_by(regions, jurisdiction, incidentid, coldate) %>%
  dplyr::summarise(labs_nresults = sum(result %in% c("Negative", "Positive"))
                  ,labs_anypos = ifelse(sum(result %in% "Positive") >= 1, 1, 0)
                  ,nfacs = sum(!is.na(unique(facilityname)), na.rm=T)
                  ,facility = paste(unique(facilityname[!is.na(facilityname)]), collapse = "; ")
                  )
bdatc <- subset(bdatc, labs_nresults >= 1)

# unique dates per incident:
table(table(bdatc$incidentid))

# generate test result status for each unique incidence testing day:
ldat <- bdatc

ldat$labs_status <- NA
ldat$labs_status[ldat$labs_anypos >= 1] <- 1
ldat$labs_status[ldat$labs_anypos == 0] <- 0

table(ldat$labs_status, useNA='ifany')

# remove unnecessary objects from workspace -------------------------------

# save(ldat, file = paste(mywd, "/../Beyer, Kirsten - WEDSS_Archive/lrein_data/", dat_date, "_phavr_labs.Rdata", sep=""))
# write.csv(ldat, file = paste(mywd, "/../Beyer, Kirsten - WEDSS_Archive/lrein_data/", dat_date, "_phavr_labs.csv", sep=""))

ldat <- data.frame(ldat)
bdat <- data.frame(bdat)

rm(list=setdiff(ls(), c(ls_save, "ldat", "bdat")))








