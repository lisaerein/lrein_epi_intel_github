# lrein_epi_intel_github

2020-07-22: Added data prep scripts

This repository includes R scripts used to prepare WEDSS data for Milwaukee County Covid-19 Epidemiology Intel Team analyses.

We download WEDSS data daily (via PHAVR), starting around 1:30pm each day. 

The PHAVR reports we use are the "LineList" file for cases and "Labs" file for individual tests.

------ LineList report --------------------------

Our PHAVR LineList report is an incident level dataset. 
- Note individuals (clientid) can have multiple disease incidents (incidentid) and each indicident may be a composite of multiple tests.

We use the LineList file for counts, analysis, etc. of confirmed cases.
- We do not use this dataset to count negative tests (we use the Labs file for this).

JURIDICTIONS: Our report includes the following jurisdictions in Milwaukee Co. and surrounding regions:
- Milwaukee
- Greenfield
- Cudahy
- Oak Creek
- West Allis 
- Greendale
- Wauwatosa
- North Shore
- Franklin
- Saint Francis
- South Milwaukee
- Hales Corners
- Ozaukee County
- Kenosha County
- Waukesha County
- Racine
- Washington County
- Walworth County
- Central Racine (not currently used in any analysis)
- Waukesha County (not currently used in any analysis)

DATES: The date we use for all analyses is the specimen collection date. 
- If unavailable (~ 5% missing), this is imputed first with the episode date then creation date if still missing.

LOCATIONS: For county specific reports, we subset the data using the 'countyfips' field. 
- If 'countyfips' is unavailable, we impute using 'countyofresidence' information. 
- Note that jurisdiction can differ from city/county fields due to health departments helping one another with follow-up, especially in the Milwaukee Co area. 
- We have found that county and city fields are better aligned with a person's reported residence. 
- City of Racine and City of Milwaukee cases are placed within the county by block ID from arcGIS

GENERAL DATA PREP STEPS:

1) Select and import the most recent LineList csv from our data directory.
2) Format and impute dates as needed.
3) Create indicators for events:
	- Cases: 'resolutionstatus' == "Confirmed"
	- Hospitalized cases: 'resolutionstatus' == "Confirmed" + 'patienthospitalized' == "True"
	- Deceased cases: 'resolutionstatus' == "Confirmed" + 'patient_died_of_this_illness' == "True"
		Note: 
		- Cases and hospitalizations are placed in time by the specimen collection date. 
		- Deaths are placed in time by the deathdate. 
		- Missing deathdates are not shown in our plots or counted in death totals in our report. 
		  Recent deaths tend to have a delay in entering deathdate by a few days.
4) Clean location and demographic variables used in analysis.
5) Exclude any incidents with a 'recordtype' of "Contact Investigation" keeping only "Disease Incident"
	- This was a recommendation from WI DHS
6) Merge in population information for Milwaukee County report to calculate rates.
7) Collapse data to get a person ('clientid') level dataset.
	- Right now our analyses are of individual people (for cases) and tests (for testing analyses)
	- If a person has more than 1 positive disease incidence we are treating that as a single case. Very few of these exist currently.
	- If a person has more than 1 positive disease incidence, we use the date of the first positive incident. 
		- If they are ever hospitalized or died (even in a subsequent positive incident) that is still counted for that person.

------ Labs report ------------------------------

Our PHAVR Labs report is a reported test result level dataset. 
- A single test can have the same result reported to WEDSS by multiple entities in the Labs file
- We de-duplicate the data by collapsing by 'incidentid' and date, as recommended by WI DHS (detailed steps outlined below)
	
JURIDICTIONS: Our report includes the following jurisdictions in Milwaukee Co. and surrounding regions:
- Milwaukee
- Greenfield
- Cudahy
- Oak Creek
- West Allis 
- Greendale
- Wauwatosa
- North Shore
- Franklin
- Saint Francis
- South Milwaukee
- Hales Corners
- Ozaukee County
- Kenosha County
- Waukesha County
- Racine
- Washington County
- Walworth County
- Central Racine (not currently used in any analysis)
- Waukesha County (not currently used in any analysis)

DATES: The date we use for all analyses is the specimen collection date. 
- If unavailable, this is imputed first with the episode date then creation date if still missing.
	
LOCATIONS: For county specific reports, we subset the data by JURISDICTION.
- Note this differs from the LineList process. 
	- Recent test results in 'staging' may not yet be geocoded or have reliable county of residence information. 
	- We have done things this way to better align with the Milwaukee Co. OEM dashboard.
		
GENERAL DATA PREP STEPS:

1) Select and import the most recent Labs csv from our data directory.
2) Format and impute dates as needed.
3) Review the 'testcode' field to check for Antibody tests and exclude any of these results (none last time I checked).
4) Collapse data by 'incidentid' and date ('coldate' in our dataset).
	- If any positive results are reported on that day count as a positive test.
	- If all are negative count as a negative test.
