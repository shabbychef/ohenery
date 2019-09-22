# /usr/bin/r
#
# Copyright 2019-2019 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav 
#
# This file is part of ohenery.
#
# ohenery is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ohenery is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with ohenery.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2019.09.04
# Copyright: Steven E. Pav, 2019
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

suppressMessages(library(docopt))       # we need docopt (>= 0.3) as on CRAN

doc <- "Usage: gen_race_data.r [-v] [-L <LOC>] [-m <MINDATE>] [-M <MAXDATE>] [-O <OUTFILE>]

-L LOC --loc=LOC                 Location of upstream data [default: ~/github/varo/data-raw/ts/summary]
-O OUTFILE --outfile=OUTFILE     Location of output data [default: race_data.csv]
-m MINDATE --mindate=MINDATE     Minimum date to include, in YYYYMMDD [default: 20190304]
-M MAXDATE --maxdate=MAXDATE     Minimum date to include, in YYYYMMDD [default: 20190324]
-v --verbose                     Be more verbose
-h --help                        show this help text"

opt <- docopt(doc)

suppressMessages({
	library(readr)
	library(dplyr)
	library(tidyr)
	library(magrittr)
})

# briscode lookup
brs_data <- jsonlite::fromJSON('~/github/varo/data-raw/ts/raw/CurrentRace_20181212T221956.json')
brs_lookups <- brs_data$CurrentRace %>%
  select(BrisCode,DisplayName) %>%
  rename(Track=DisplayName)

dlist <- dir(path=opt$loc,pattern='ts_sum.+\\.csv',full.names=TRUE)
dateus <- as.Date(gsub('^.+ts_sum_(\\d{8}).csv','\\1',dlist),format='%Y%m%d')
isok   <- (dateus >= as.Date(opt$mindate,format='%Y%m%d')) & (dateus <= as.Date(opt$maxdate,format='%Y%m%d')) 
dlist  <- dlist[isok]

ctype <- readr::cols(
 .default=col_double(),
  date = col_date(format = ""),
  Track = col_character(),
  Type = col_character(),
  RaceNum = col_integer(),
  ProgramNumber = col_integer(),
  finish = col_integer(),
  PostPosition = col_integer(),
  HorseName = col_character(),
  Jockey = col_character(),
  Trainer = col_character(),
  Medication = col_character(),
  ApprenticeWeight = col_integer(),
  EquipmentIndicator = col_integer(),
  Sire = col_character(),
  Dam = col_character(),
  Age = col_integer(),
  Sex = col_character(),
  MTO = col_character(),
  Weight_lbs = col_integer(),
  PostTime = col_character(),
  CorrectedPostTime = col_datetime(format = ""),
  SurfaceText = col_character(),
  Purse = col_character(),
  Turf = col_character(),
  Dirt = col_character(),
  Country = col_character(),
  AveragePaceE1 = col_integer(),
  AveragePaceE2 = col_integer(),
  AveragePaceLp = col_integer(),
  AverageSpeed = col_integer(),
  MaxClaimPrice = col_integer(),
  AltPurse = col_integer(),
  HorseId = col_integer(),
  JockeyWeight = col_integer(),
  SpeedPoints = col_integer(),
  TrainerId = col_integer(),
  JockeyId = col_integer()
)

race_data <- lapply(dlist,readr::read_csv,col_type=ctype) %>%
	bind_rows() %>%
	select(-TrainerId,-JockeyId,-SpeedPoints,-HorseId,
				 -matches('^Average(Pace(E1|E2|Lp)|Speed)$'),
				 -ApprenticeWeight,-EquipmentIndicator,
				 -Purse,  # give them numbers?
				 -MedWeight) %>%  # redundant that one
	rename(BrisCode=Track) %>%
	mutate(BrisCode=tolower(BrisCode)) %>%
	arrange(date,BrisCode,RaceNum,finish,PostPosition) %>%
	rename(Date=date,Finish=finish) %>%
	mutate(EventId=group_indices(.,Date,BrisCode,RaceNum)) %>%
	mutate(HorseId=group_indices(.,HorseName,Sex,Sire,Dam)) %>%
	#mutate(JockeyId=group_indices(.,Jockey)) %>%
	#mutate(TrainerId=group_indices(.,Trainer)) %>%
	left_join(brs_lookups,by='BrisCode') %>%
	mutate(TrackId=group_indices(.,BrisCode)) %>%
	select(EventId,
				 #Date,
				 #Country,BrisCode,Track,
				 TrackId,
				 Type,RaceNum,
				 #PostTime,
				 CorrectedPostTime,
				 Yards,SurfaceText,
				 #Turf,Dirt,
				 #MaxClaimPrice,AltPurse,
				 HorseName,HorseId,
				 #Sire,Dam,
				 Age,Sex,Weight_lbs, #ProgramNumber,
				 PostPosition,
				 Medication,   #MTO,  # dunno what this is.
				 MorningLine,WN_pool,PL_pool,SH_pool,
				 #WN_payoff,PL_payoff,SH_payoff,
				 #Jockey,JockeyId,JockeyWeight,Trainer,TrainerId
				 Finish) %>%
	#rename(PurseDollars=AltPurse) %>%
	mutate(Sex=toupper(Sex)) 

race_data %>% readr::write_csv(opt$outfile)
usethis::use_data(race_data,overwrite=TRUE)

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
