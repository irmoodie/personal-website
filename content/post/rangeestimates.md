---
title: "Estimating Species Ranges using R, GBIF and GADM"
date: 2021-01-19T15:43:10Z
Description: ""
Tags: []
Categories: [blogpost, guide, R]
DisableComments: false
---

While working on a phylogenetics project at Uppsala University, I ended up getting a bit carried away with what was supposed to be a tiny part of the overall project. However, I found it super interesting and felt it worthwhile to share more widely.

#### Some context

During the course *Evolutionary Patterns* the final project consisted of constructing a phylogeny to answer a specific chosen question. My group was interested in assessing if two large moss orders were monophyletic or not, and if some understudied "weird" taxa were supported as belonging to either order. Part of our initial interest stemmed from the biogeography of these species, and this is where the interest in gathering range estimate data came from. However, we dropped this aspect of the project as it became clear it would require much more time than we had to do it justice.

When I say range here, I am referring to the geographic distribution of a species. When I say range size, I am generally referring to the summed area of the regions where a species is found.

#### The problem

Finding a reliable and, importantly, quantified estimate of the range of understudied moss species in the literature proved difficult. To overcome this, I put together a script that would use occurrence data from the Global Biodiversity Information Facility ([GBIF](https://www.gbif.org/)) and shapefile data from the Database of Global Administrative Areas ([GADM](https://gadm.org/)) to produce an estimate of range size. I need to state here I am not the first to do this by any means, and much of my code has been adapted from other sources (see citations), I just wish to demonstrate how easy something like this is to do, and the power of having access to open public databases.

#### The solution

The pipeline looks a little something like this:

```
Input: 
		A list of taxa names for the species of interest

Process: 	
		1. Using the species names, find the usage keys for each by querying GBIF
		2. Using the usage keys, submit a download request for all entries for each key
		3. Using the location data associated with each entry, create a list of all regions that a species is present in
		4. Using shapefiles from GADM, calculate the areas for these regions and sum for each species
			
Output:
		A table with species range estimates for each taxa input
```

I chose to implement this in R, but it could be easily replicated using python as many of the packages used have their python counterpart (e.g.`rgbif` = `pygbif`). I also use a combination of `base` R and `tidyverse`, so my apologies to the purists out there. Okay let's start.

I have a file with the names of each of our study taxa (`species.csv`), formatted as a single column with the heading `species` and each species name written in the format `Homo_sapiens`. I begin by clearing the environment and loading the required packages:

```R
rm(list=ls()) # clear environment

library(rgbif) # for pulling data from GBIF
library(taxize) # useful for dealing with taxa names
library(magrittr) # for tee operator
library(tidyverse) # load last to avoid conflicts
```

First we need to get the usage keys for each one of our taxa. We do this by reading in the file `species.csv` with `readr::read_csv()`.  Using the pipe operator `%>%` I "pipe" it into `dplyr::mutate()` so that we can remove the underscore and replace it with a space, as this is the format GBIF will be expecting for species names. Next I use `pull()` to create a list of species names, and then use `rgbif::get_gbifid_()` to get the GBIF backbone taxon IDs (the strings of numbers that uniquely identifies a species and all its synonym names in the GBIF database) via the GBIF API. If we stopped at this step, we'd have have a list of dataframes (one for each species we searched), each with 24 columns of data (including taxonomic data, usage keys, metadata on the search, etc). Next we add in the species name as a column, as otherwise when we use `dplyr::bind_rows()` we would lack  a column that would exactly match the names we provided to GBIF. This is going to make any further work with this dataset much easier. Using a "Tee Pipe" `%T>% `  and `readr::write_tsv()` , we will save our current dataset to a file, as it saves us having to query the GBIF database in the future when we need this information. The tee pipe is very hand for things like this, think of it like a little side action, a "but also could you quickly do this", and is not acknowledged by the rest of the code after it is done.  Finally we do some filtering using `dplyr::filter()` to make sure our dataset is what we want, and then again use `pull()` to create a list of the usage keys, which are what we will use to tell GBIF what species we want the entries of in the next step.

```R
taxon_keys <- read_csv("species.csv") %>% # read in the file containing species names
  mutate(species = str_replace(species, "_", " ")) %>% # change "Homo_sapiens" to "Homo sapiens"
  pull("species") %>% # create a list of species names
  get_gbifid_(method="backbone") %>% # get the taxon IDs
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add col for name
  bind_rows() %T>%
  write_tsv(path = "all_matches.tsv") %>% # output the full file for records
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # ensure only exact and accepted matches are taken to avoid duplicates
  filter(kingdom == "Plantae") %>% # just in case
  pull(usagekey) # makes a list of the usage keys for submitting request
```

The next step is very simple. We use the function `rgbif::occ_download()` to query the GBIF database using our `taxon_keys`. The download request will then appear in the downloads section of your GBIF account.

```R
occ_download(
  pred_in("taxonKey", taxon_keys), # use taxon keys to id species of interest
  format = "SIMPLE_CSV", # specifies the format we want in the download, and the simple CSV contains enough for our purposes
  user=*******,pwd=********,email=*********** # you need to put in your own user info
)
```

Once the request has been completed, and we have our downloaded dataset. For my 32 moss taxa the dataset ended up being approx. 260MB, so be aware they can get big quickly. Next, we move onto estimating our species ranges from the entries. I had separated this part into a separate script so I could reuse the query script for other projects. We start again by clearing the environment and loading in the required packages. Here I make use of the package `raster` for working with shapefiles, as well as a package called `countrycode`, which proved to be a very handy tool for converting two letter code country codes (SE) to three letter country codes (SWE), and vice versa.

```R
rm(list=ls()) # clear environment

library(raster) # for dealing with shapefiles
library(countrycode) # easy way to convert iso2c to iso3c and viceversa
library(tidyverse) 
```

Next we read in our downloaded dataset that I had saved as `gbif.csv` using `readr::read_tsv()` (the downloaded file was tab separated, not comma), and then filter out the entries that were useless for our purposes, as well as removing those that likely do not represent the natural distribution of the species (keeping only those that have coordinates associated with the entry).

```R
data <- read_tsv("gbif.csv") %>% # import dataset downloaded from GBIF
  filter(occurrenceStatus == "PRESENT", # remove absence entries
         decimalLatitude != "NA", # remove entries with missing lat/long
         decimalLongitude != "NA") 
```

At this stage we need to download our required shapefiles from GADM. In my case, I was fairly sure we would need shapefiles for all countries, but if we wanted to check which countries are present in the dataset, simply use `dplyr::group_by()` and `count()` to get a list of all countries present, as below.

```R
data %>%
	group_by(countryCode) %>%
  	count()
```

We will download the whole world shapefiles that are separated by level from [here](https://gadm.org/download_world.html). We are really only interested in level 0, which is at the level of the country. As to not bias our estimates towards species in very large countries we also need to download the level 1 (state/province) shapefiles for each large country present in the dataset. For me, that was Australia, Brazil, Russia, China, USA, India and Canada. Once we have these files, I extracted them into the R working directory. Let's begin with the large countries. The goal here is to import the shapefile for the large country, get a list of all provinces that we need to measure for the dataset, measure those provinces and sum the totals area of all provinces a species is found in. We do this for each country in the list by making use of a `for` loop. We start by creating a list of all the ISO3C codes for each of our large countries (`largecountries`), and an empty data frame that we will add the range estimates to (`rangeDataset`). Next we start our loop by writing for each country in the list `largecountries` please do the following:

1. Import the required shapefile from the folder that matches the countries ISO3C code (e.g. for CHN it will be found in `gadm36_CHN_shp/gadm36_CHN_1.shp`) and store it as `shp`.
2. Filter our dataset (`data`) to only include entries that are for the large country (As GBIF uses ISO2C, we use the function `countrycode::countrycode()` to translate the ISO3C code we are using) and for entries that have a matching region name in the GADM shapefile `shp` (you might want to do some manual checking to ensure the names match up, especially when you get into regions like Russia). Group these entries by `species `and `stateProvince` to get all unique combinations, and save it as `countrydata`. This gives us a data frame with three columns, `species`, `stateProvince` and `count`.
3. Create a vector called `prov` that lists all the unique provinces found in `countrydata`, and create an empty vector called `area` that has the same length as `prov`.
4.  For each province in `prov` , calculate the area of the province in km^2 using `raster::area()` and add it to the vector `area`. 
5. Combine the two vectors `prov` and `area` into a data frame, which I have called `area`. This now gives us two columns, one with the province name, and the other with its area. We also rename the column `prov` to `stateProvince` so we can easily join it in the next step.
6. Using `dplyr::left_join()` add a column of areas for each province in `countrydata`.
7. Calculate the range of each species by summing the areas of provinces where each species is found, and call this data frame `range`.
8. Using `rbind()`, add `range` to our data frame `rangeDataset`.
9. Remove the variables used in the loop.

The loop then goes to the next large country in the list and does the same, and we are left with the range estimates for our species within these countries.

```R
largecountries <- c("AUS", "USA", "CAN", "IND", "RUS", "BRA", "CHN") # list large countries by iso3c code

rangeDataset <- data.frame(species = character(),
                           range = double()) # create empty dataset to enter range data into

for (i in seq_along(largecountries)) { # for each country in the list do the following:
  
  shp <- shapefile(paste0("gadm36_", largecountries[i], "_shp/gadm36_", largecountries[i], "_1.shp", sep = "")) # import country shape file that contains province/states
  
  countrydata <- data %>%
    filter(countryCode == countrycode(largecountries[i], origin = "iso3c", destination = "iso2c"), 
           stateProvince %in% shp$NAME_1) %>% # filter to data with correct prov names only
    group_by(species, stateProvince) %>% # group by unique combinations of species and provs
    count() # add a count (not actually used anymore)
  
  prov <- unique(countrydata$stateProvince) # list provs in dataset
  area <- vector(length = length(prov)) # make empty vector of same length as prov vector
  
  for (i in 1:length(prov)){ # for the number of provs do the following:
    area[i] <- raster::area(shp[which(shp$NAME_1 == paste0(prov[i])),])/1000000 # measure area of prov
  }
  
  area <- data.frame(prov, area) # match up areas with provinces
  area <- rename(area, stateProvince = prov) # rename variable for joining
  countrydata <- left_join(countrydata, area, by = "stateProvince") # join areas to gbif data
  
  range <- countrydata %>% 
    group_by(species) %>%
    summarise(range = sum(area)) # add up all areas of countires in which a species is located
  
  rangeDataset <- rbind(rangeDataset, range) # add ranges to rangeDataset
  
  rm(area, countrydata, prov, shp, range) # remove temp variables before next loop
  
}

rm(largecountries, i)
```

Next, the rest of the world. We'll use a very similar approach as above, but this time we don't need to separate the countries into regions, and we can use the same shapefile throughout. 

``` R
shp <- shapefile("gadm36_levels_shp/gadm36_0.shp") #import world shape file

worlddata <- data %>% # filter out countires already done
  filter(countryCode != "AU",
         countryCode != "CA" ,
         countryCode != "IN" ,
         countryCode != "RU" ,
         countryCode != "US",
         countryCode != "BR",
         countryCode != "CN") %>%
  group_by(species, countryCode) %>%
  count()

prov <- unique(worlddata$countryCode, na.rm=TRUE) # get list of countryCodes in dataset
prov <- countrycode(prov, origin = "iso2c", destination = "iso3c") # convert 2 letter codes to 3 letter, as used in the shapefile
prov <- na.omit(prov) # remove NAs (such as ZZ, cannot estimate unknown country)
area <- vector(length = length(prov)) # make empty vector for loop

for (i in 1:length(prov)){ # almost same loop as above
  area[i] <- raster::area(shp[which(shp$GID_0==paste0(prov[i])),])/1000000
}

prov <- countrycode(prov, origin = "iso3c", destination = "iso2c") # convert codes back
area <- data.frame(prov, area) # match up areas with provinces
area <- rename(area, countryCode = prov) # rename varaible for join
worlddata <- left_join(worlddata, area, by = "countryCode") # join areas to gbif data
worlddata <- worlddata %>%
  filter(area != "NA") # remove NAs

range <- worlddata %>%
  group_by(species) %>%
  summarise(range = sum(area)) # add up all areas of countires in which a species is located

rangeDataset <- rbind(rangeDataset, range) # add ranges to rangeDataset

rangeDataset <- rangeDataset %>%
  group_by(species) %>%
  summarise(range = sum(range)) # summarise range dataset to list of species and ranges

rm(area, worlddata, prov, shp, range, i) # remove temp variables

write_csv(rangeDataset, file = "rangedataset.csv") # save dataset to a file
```

And that's it! We've estimated the range sizes of all of our species. I think this method can be used for many group of species as a quick and easy way to estimate the range sizes, and possibly population sizes? This maybe has some added assumptions, but I would expect the two to correlate quite closely when looking within a group of closely related taxa.

Quickly as a last step, we can plot our dataset using `ggplot2`.

```R
ggplot(data = rangeDataset, aes(y = range/1000000, x = species)) +
  geom_bar(stat = "identity") +
  ylab(label = c(expression("Range size in millions of" ~km^2))) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  coord_flip() +
  theme_classic() +
  theme(axis.title.y=element_blank())
```

![Range Estimates](/images/rangeestimates.png)

---

#### Citations

For code and troubleshooting related to GBIF see:

[Downloading occurrences from a long list of species in R and Python - Waller and Grosjean (2019)](https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/)

Original idea (and some code) for estimating species ranges came from this paper by Blow et al. See Supplementary materials and the associated GitHub repo:

[Blow R, Willink B, Svensson EI. 2020. A molecular phylogeny of forktail damselflies (genus Ischnura) reveals a dynamic macroevolutionary history of female colour polymorphisms. bioRxiv 2020.06.06.137828.](https://www.biorxiv.org/content/10.1101/2020.06.06.137828v1.full) 

[GitHub: rachelblow/IschnuraPolymorphismEvolutionaryHistory](https://github.com/rachelblow/IschnuraPolymorphismEvolutionaryHistory)